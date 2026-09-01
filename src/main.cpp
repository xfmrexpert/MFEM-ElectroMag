// Copyright (c) 2026 T. C. Raymond
// SPDX-License-Identifier: MIT

#include <filesystem>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include "mfem.hpp"
#include "build_info.hpp"
#include "config/input_parser.hpp"
#include "solvers/physics_solver.hpp"
#include "solvers/solver_factory.hpp"
#include "config/config_validator.hpp"
#include "io/status_reporter.hpp"

namespace {

void PrintUsage(const char* prog) {
    std::cerr <<
        "Usage: " << prog << " <config.json> [options]\n"
        "Options:\n"
        "  --results-path <directory> Override Gmsh results output directory\n"
        "  --export-refine <N>        Refinement factor for export mesh (default = solve order)\n"
        "  --export-vector-space <L2|H1>  Reserved; L2 is currently the only supported choice\n"
        "  --verbosity <0|1|2>        0=status/timing, 1=solver output, 2=diagnostics\n"
        "  --machine-readable         Emit flushed JSON Lines progress on stdout\n"
        "  --version                  Print version/build information and exit\n"
        "  -h, --help                 Show this help\n";
}

} // namespace

int main(int argc, char *argv[]) {
    // Convert MFEM internal errors (MFEM_ASSERT / MFEM_VERIFY / MFEM_ABORT)
    // into catchable mfem::ErrorException objects so a malformed mesh does
    // not call std::abort() from deep inside the Mesh constructor.
    // Requires MFEM built with MFEM_USE_EXCEPTIONS=ON (see top-level
    // CMakeLists.txt).
#ifdef MFEM_USE_EXCEPTIONS
    mfem::set_error_action(mfem::MFEM_ERROR_THROW);
#endif
    StatusReporter& reporter = StatusReporter::Global();
    bool machine_readable_requested = false;
    for (int i = 1; i < argc; ++i) {
        if (std::string(argv[i]) == "--machine-readable") {
            machine_readable_requested = true;
            break;
        }
    }
    if (machine_readable_requested) {
        reporter.SetFormat(StatusReporter::Format::JsonLines);
    }

    try {
        std::string config_file = "config.json";
        std::string cli_results_path;
        int  cli_export_refine = 0;        // 0 = unset
        std::string cli_vector_space;
        int cli_verbosity = 0;
        bool verbosity_explicit = false;

        for (int i = 1; i < argc; ++i) {
            std::string a = argv[i];
            auto need_value = [&](const std::string& opt) -> std::string {
                if (i + 1 >= argc) {
                    throw std::runtime_error(opt + " requires a value");
                }
                return std::string(argv[++i]);
            };

            if (a == "-h" || a == "--help") {
                PrintUsage(argv[0]);
                return 0;
            } else if (a == "--version") {
                std::cout << build_info::Describe() << '\n';
                return 0;
            } else if (a == "--results-path") {
                cli_results_path = need_value(a);
            } else if (a == "--export-refine") {
                cli_export_refine = std::stoi(need_value(a));
                if (cli_export_refine < 1) {
                    throw std::runtime_error("--export-refine must be >= 1");
                }
            } else if (a == "--export-vector-space") {
                cli_vector_space = need_value(a);
                if (cli_vector_space != "L2" && cli_vector_space != "H1") {
                    throw std::runtime_error("--export-vector-space must be L2 or H1");
                }
                if (cli_vector_space == "H1") {
                    reporter.Warning("--export-vector-space H1 is not yet implemented; using L2.");
                }
            } else if (a == "--verbosity") {
                cli_verbosity = std::stoi(need_value(a));
                StatusReporter::VerbosityFromInt(cli_verbosity);
                verbosity_explicit = true;
            } else if (a == "--machine-readable") {
                machine_readable_requested = true;
            } else if (!a.empty() && a[0] == '-') {
                throw std::runtime_error("Unknown option: " + a);
            } else {
                config_file = a;
            }
        }

        if (machine_readable_requested && !verbosity_explicit) {
            cli_verbosity = 1;
        }
        reporter.SetVerbosity(StatusReporter::VerbosityFromInt(cli_verbosity));
        if (machine_readable_requested) {
            mfem::out.SetStream(reporter.SolverOutput());
        }

        // Report build identity up front so a run can always be traced back to
        // a specific binary; a stale executable is otherwise only detectable by
        // the confusing downstream errors it produces. Routed through the
        // reporter so --machine-readable still emits well-formed JSON Lines.
        reporter.Status(build_info::Describe());

        // Shared Infrastructure
        std::error_code config_ec;
        auto config_abs = std::filesystem::weakly_canonical(config_file, config_ec);
        reporter.Diagnostic("Config file: "
            + (config_ec ? std::filesystem::absolute(config_file).string()
                         : config_abs.string()));

        std::unique_ptr<InputParser> parser_ptr;
        {
            auto operation = reporter.Start("configuration loading");
            parser_ptr = std::make_unique<InputParser>(config_file);
        }
        InputParser& parser = *parser_ptr;

        // Inject CLI overrides into the shared json so downstream parsing picks
        // them up (solvers re-run InputParser internally inside Setup()).
        // const_cast is safe: parser owns the json when constructed from a path.
        json& mutable_json = const_cast<json&>(parser.config);
        if (!mutable_json.contains("simulation") || !mutable_json["simulation"].is_object()) {
            mutable_json["simulation"] = json::object();
        }
        if (!cli_results_path.empty()) {
            mutable_json["simulation"]["results_path"] = cli_results_path;
        }
        if (cli_export_refine > 0) {
            mutable_json["simulation"]["export_refine"] = cli_export_refine;
        }

        // Validate Configuration (schema and basic semantics before decoding)
        ConfigValidator validator;
        {
            auto operation = reporter.Start("configuration validation");
            validator.ValidateOrThrow(parser.config);
        }

        ProblemConfig config = parser.GetProblemConfig();

        // Load Mesh
        std::error_code mesh_ec;
        auto mesh_abs = std::filesystem::weakly_canonical(config.MeshPath, mesh_ec);
        reporter.Diagnostic("Mesh file: "
            + (mesh_ec ? std::filesystem::absolute(config.MeshPath).string()
                       : mesh_abs.string()));

        // MFEM reports a missing/unreadable file through the same error path
        // as a malformed one, so check it up front to keep the two failures
        // distinguishable for the user.
        std::error_code exists_ec;
        if (!std::filesystem::is_regular_file(config.MeshPath, exists_ec)) {
            throw std::runtime_error(
                "Mesh file not found or not readable: '" + config.MeshPath + "'");
        }

        // Load without auto-fix so we can diagnose bad meshes uniformly in
        // Debug and Release (MFEM_ASSERT is a no-op in Release, which would
        // otherwise let an inconsistent mesh through silently).
        //
        // NOTE: the mfem::Mesh(filename, ...) constructor internally calls
        // Load() -> Finalize() -> CheckBdrElementOrientation(). For meshes
        // with orphan boundary elements (boundary edges whose endpoints are
        // not shared by any 2D element), CheckBdrElementOrientation indexes
        // faces_info[-1] and triggers MFEM_ASSERT. Requires MFEM built with
        // MFEM_USE_EXCEPTIONS=ON (enforced below) so that becomes a catchable
        // mfem::ErrorException instead of an std::abort().
#ifndef MFEM_USE_EXCEPTIONS
#error "MFEM must be built with MFEM_USE_EXCEPTIONS=ON; see top-level CMakeLists.txt"
#endif
        std::unique_ptr<mfem::Mesh> mesh_ptr;
        {
            auto operation = reporter.Start("mesh loading");
            try {
                mesh_ptr = std::make_unique<mfem::Mesh>(
                    config.MeshPath, /*generate_edges=*/1, /*refine=*/0,
                    /*fix_orientation=*/false);
            }
            catch (const mfem::ErrorException& e) {
                throw std::runtime_error(
                    "MFEM could not load mesh '" + config.MeshPath + "'. Check "
                    "that the file format is supported and that 2D elements use "
                    "counter-clockwise winding with boundary lines whose nodes "
                    "lie on element edges.\nMFEM detail: " + std::string(e.what()));
            }

            // Orientation is reported (not repaired) so a bad mesh fails loudly
            // at its source rather than being silently patched every run.
            const int bad_el  = mesh_ptr->CheckElementOrientation(/*fix_it=*/false);
            const int bad_bdr = mesh_ptr->CheckBdrElementOrientation(/*fix_it=*/false);
            if (bad_el != 0 || bad_bdr != 0) {
                throw std::runtime_error(
                    "Invalid mesh '" + config.MeshPath + "': "
                    + std::to_string(bad_el) + " mis-oriented element(s), "
                    + std::to_string(bad_bdr)
                    + " mis-oriented or orphan boundary element(s). "
                    "Regenerate the mesh with counter-clockwise 2D element "
                    "winding and boundary lines whose nodes lie on element edges.");
            }
        }

        mfem::Mesh& mesh = *mesh_ptr;

        // Validate Configuration (with mesh for attribute checking)
        {
            auto operation = reporter.Start("configuration validation");
            validator.ValidateOrThrow(parser.config, &mesh);
        }

        // Factory Logic - Create Solver
        std::unique_ptr<PhysicsSolver> solver;
        {
            auto operation = reporter.Start("solver creation");
            solver = SolverFactory::Instance().Create(mesh, config);
        }

        // Execution
        {
            auto operation = reporter.Start("solver setup");
            solver->Setup();
        }
        {
            auto operation = reporter.Start("solution process");
            solver->Run();
        }
        {
            auto operation = reporter.Start("analysis result writing");
            solver->SaveAnalysis();
        }

        return 0;
    }
    catch (const std::exception& e) {
        reporter.Error(e.what());
        return 1;
    }
    catch (...) {
        reporter.Error("Unknown exception occurred");
        return 2;
    }
}