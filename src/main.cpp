// Copyright (c) 2026 T. C. Raymond
// SPDX-License-Identifier: MIT

#include <filesystem>
#include <iostream>
#include <memory>
#include <string>
#include "input_parser.hpp"
#include "physics_solver.hpp"
#include "solver_factory.hpp"
#include "config_validator.hpp"

namespace {

void PrintUsage(const char* prog) {
    std::cerr <<
        "Usage: " << prog << " <config.json> [options]\n"
        "Options:\n"
        "  --results-file <path>      Override Gmsh MSH 2.2 results output path\n"
        "  --export-refine <N>        Refinement factor for export mesh (default = solve order)\n"
        "  --export-vector-space <L2|H1>  Reserved; L2 is currently the only supported choice\n"
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

    try {
        std::string config_file = "config.json";
        std::string cli_results_file;
        int  cli_export_refine = 0;        // 0 = unset
        std::string cli_vector_space;

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
            } else if (a == "--results-file") {
                cli_results_file = need_value(a);
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
                    std::cerr << "Warning: --export-vector-space H1 is not yet "
                                 "implemented; using L2." << std::endl;
                }
            } else if (!a.empty() && a[0] == '-') {
                throw std::runtime_error("Unknown option: " + a);
            } else {
                config_file = a;
            }
        }

        // 1. Shared Infrastructure
        std::error_code config_ec;
        auto config_abs = std::filesystem::weakly_canonical(config_file, config_ec);
        std::cout << "Config file: "
                  << (config_ec ? std::filesystem::absolute(config_file).string()
                                : config_abs.string())
                  << std::endl;

        InputParser parser(config_file);

        // Inject CLI overrides into the shared json so downstream parsing picks
        // them up (solvers re-run InputParser internally inside Setup()).
        // const_cast is safe: parser owns the json when constructed from a path.
        json& mutable_json = const_cast<json&>(parser.config);
        if (!mutable_json.contains("simulation") || !mutable_json["simulation"].is_object()) {
            mutable_json["simulation"] = json::object();
        }
        if (!cli_results_file.empty()) {
            mutable_json["simulation"]["results_file"] = cli_results_file;
        }
        if (cli_export_refine > 0) {
            mutable_json["simulation"]["export_refine"] = cli_export_refine;
        }

        ProblemConfig config = parser.GetProblemConfig();

        // 2. Validate Configuration (basic validation before mesh loading)
        ConfigValidator validator;
        validator.ValidateOrThrow(parser.config);

        // 3. Load Mesh
        std::error_code mesh_ec;
        auto mesh_abs = std::filesystem::weakly_canonical(config.MeshPath, mesh_ec);
        std::cout << "Mesh file: "
                  << (mesh_ec ? std::filesystem::absolute(config.MeshPath).string()
                              : mesh_abs.string())
                  << std::endl;

        // Load without auto-fix so we can diagnose bad meshes uniformly in
        // Debug and Release (MFEM_ASSERT is a no-op in Release, which would
        // otherwise let an inconsistent mesh through silently).
        //
        // NOTE: the mfem::Mesh(filename, ...) constructor internally calls
        // Load() -> Finalize() -> CheckBdrElementOrientation(). For meshes
        // with orphan boundary elements (boundary edges whose endpoints are
        // not shared by any 2D element), CheckBdrElementOrientation indexes
        // faces_info[-1] and triggers MFEM_ASSERT. With MFEM_USE_EXCEPTIONS
        // that becomes a catchable mfem::ErrorException; without it, it
        // would std::abort(). We catch it here to produce an actionable
        // diagnostic instead of crashing.
        std::unique_ptr<mfem::Mesh> mesh_ptr;
        try {
            mesh_ptr = std::make_unique<mfem::Mesh>(
                config.MeshPath, /*generate_edges=*/1, /*refine=*/0,
                /*fix_orientation=*/false);
        }
#ifdef MFEM_USE_EXCEPTIONS
        catch (const mfem::ErrorException& e) {
            throw std::runtime_error(
                "Invalid mesh '" + config.MeshPath + "': MFEM rejected it "
                "during load. This usually means orphan boundary elements "
                "(boundary edges whose endpoints are not shared by any 2D "
                "element) or mis-oriented elements. Regenerate the mesh "
                "with counter-clockwise 2D element winding and boundary "
                "lines whose nodes lie on element edges.\nMFEM detail: "
                + std::string(e.what()));
        }
#endif
        mfem::Mesh& mesh = *mesh_ptr;

        int bad_el  = mesh.CheckElementOrientation(false);
        int bad_bdr = mesh.CheckBdrElementOrientation(false);
        if (bad_el != 0 || bad_bdr != 0) {
            throw std::runtime_error(
                "Invalid mesh '" + config.MeshPath + "': "
                + std::to_string(bad_el) + " mis-oriented element(s), "
                + std::to_string(bad_bdr)
                + " mis-oriented or orphan boundary element(s). "
                "Regenerate the mesh with counter-clockwise 2D element "
                "winding and boundary lines whose nodes lie on element edges.");
        }

        mesh.Finalize(/*refine=*/true, /*fix_orientation=*/true);

        // 4. Validate Configuration (with mesh for attribute checking)
        validator.ValidateOrThrow(parser.config, &mesh);

        // 5. Factory Logic - Create Solver
        auto solver = SolverFactory::Instance().Create(config.PhysicsType, mesh, parser.config);

        // 6. Execution
        solver->Setup();
        solver->Run();
        // Save is now called inside of Run()
        solver->SaveAnalysis();

        return 0;
    }
    catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    catch (...) {
        std::cerr << "Error: Unknown exception occurred" << std::endl;
        return 2;
    }
}