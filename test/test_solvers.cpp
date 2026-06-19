// Copyright (c) 2026 T. C. Raymond
// SPDX-License-Identifier: MIT

#include <catch2/catch_test_macros.hpp>
#include "../src/input_parser.hpp"
#include "../src/electrostatic_solver.hpp"
#include "../src/magnetostatic_solver.hpp"
#include "../src/magnetoquasistatic_solver.hpp"
#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <array>
#include <cmath>
#include <cstddef>

namespace fs = std::filesystem;

// Helper function to create a minimal test mesh file
void CreateTestMesh(const std::string& filename) {
    std::ofstream mesh_file(filename);
    // Create a simple 2D axisymmetric mesh (single triangle element)
    mesh_file << "MFEM mesh v1.0\n\n";
    mesh_file << "dimension\n2\n\n";
    mesh_file << "elements\n1\n";
    mesh_file << "1 2 0 1 2\n\n";  // Triangle with attribute 1
    mesh_file << "boundary\n3\n";
    mesh_file << "1 1 0 1\n";       // Boundary segment with attribute 1
    mesh_file << "1 1 1 2\n";       // Boundary segment with attribute 1
    mesh_file << "1 1 2 0\n\n";     // Boundary segment with attribute 1
    mesh_file << "vertices\n3\n2\n";
    mesh_file << "0.0 0.0\n";
    mesh_file << "1.0 0.0\n";
    mesh_file << "0.5 1.0\n";
    mesh_file.close();
}

TEST_CASE("ElectrostaticSolver can be constructed", "[solvers]") {
    std::string mesh_file = "test_electrostatic.mesh";
    CreateTestMesh(mesh_file);

    json test_config = {
        {"simulation", {
            {"physics", "electrostatics"},
            {"mesh", mesh_file},
            {"order", 1},
            {"geometry", "axisymmetric"}
        }},
        {"materials", json::array({
            {
                {"name", "dielectric"},
                {"attributes", {1}},
                {"properties", {{"epsilon_r", 2.0}}}
            }
        })},
        {"boundaries", json::array({
            {
                {"name", "ground"},
                {"attributes", {1}},
                {"type", "Dirichlet"},
                {"value", 0.0}
            }
        })}
    };

    InputParser parser(test_config);
    mfem::Mesh mesh(mesh_file.c_str(), 1, 1);

    REQUIRE_NOTHROW(ElectrostaticSolver(mesh, parser.config));

    // Cleanup
    fs::remove(mesh_file);
}

TEST_CASE("MagnetostaticSolver can be constructed", "[solvers]") {
    std::string mesh_file = "test_magnetostatic.mesh";
    CreateTestMesh(mesh_file);

    json test_config = {
        {"simulation", {
            {"physics", "magnetostatics"},
            {"mesh", mesh_file},
            {"order", 1},
            {"geometry", "axisymmetric"}
        }},
        {"materials", json::array({
            {
                {"name", "iron"},
                {"attributes", {1}},
                {"properties", {{"mu_r", 1000.0}}}
            }
        })},
        {"sources", json::array({
            {
                {"name", "coil"},
                {"attributes", {1}},
                {"type", "CurrentDensity"},
                {"value", 1000.0}
            }
        })},
        {"boundaries", json::array({
            {
                {"name", "far_field"},
                {"attributes", {1}},
                {"type", "Dirichlet"},
                {"value", 0.0}
            }
        })}
    };

    InputParser parser(test_config);
    mfem::Mesh mesh(mesh_file.c_str(), 1, 1);

    REQUIRE_NOTHROW(MagnetostaticSolver(mesh, parser.config));

    // Cleanup
    fs::remove(mesh_file);
}

TEST_CASE("MagnetoquasistaticSolver can be constructed", "[solvers]") {
    std::string mesh_file = "test_mqs.mesh";
    CreateTestMesh(mesh_file);

    json test_config = {
        {"simulation", {
            {"physics", "magnetoquasistatics"},
            {"mesh", mesh_file},
            {"order", 1},
            {"geometry", "axisymmetric"},
            {"frequency", 60.0}
        }},
        {"materials", json::array({
            {
                {"name", "conductor"},
                {"attributes", {1}},
                {"properties", {
                    {"mu_r", 1.0},
                    {"sigma", 5.8e7}
                }}
            }
        })},
        {"sources", json::array({
            {
                {"name", "coil"},
                {"attributes", {1}},
                {"type", "CurrentDensity"},
                {"value", 1000.0}
            }
        })},
        {"boundaries", json::array({
            {
                {"name", "far_field"},
                {"attributes", {1}},
                {"type", "Dirichlet"},
                {"value", 0.0}
            }
        })}
    };

    InputParser parser(test_config);
    mfem::Mesh mesh(mesh_file.c_str(), 1, 1);

    REQUIRE_NOTHROW(MagnetoquasistaticSolver(mesh, parser.config));

    // Cleanup
    fs::remove(mesh_file);
}

TEST_CASE("Solver factory logic works correctly", "[solvers]") {
    std::string mesh_file = "test_factory.mesh";
    CreateTestMesh(mesh_file);

    SECTION("electrostatics") {
        json config = {
            {"simulation", {
                {"physics", "electrostatics"},
                {"mesh", mesh_file},
                {"order", 1}
            }},
            {"materials", json::array()},
            {"boundaries", json::array()}
        };

        InputParser parser(config);
        mfem::Mesh mesh(mesh_file.c_str(), 1, 1);

        std::string type = parser.config["simulation"]["physics"];
        REQUIRE(type == "electrostatics");
    }

    SECTION("magnetostatics") {
        json config = {
            {"simulation", {
                {"physics", "magnetostatics"},
                {"mesh", mesh_file},
                {"order", 1}
            }},
            {"materials", json::array()},
            {"sources", json::array()},
            {"boundaries", json::array()}
        };

        InputParser parser(config);
        mfem::Mesh mesh(mesh_file.c_str(), 1, 1);

        std::string type = parser.config["simulation"]["physics"];
        REQUIRE(type == "magnetostatics");
    }

    SECTION("magnetoquasistatics") {
        json config = {
            {"simulation", {
                {"physics", "magnetoquasistatics"},
                {"mesh", mesh_file},
                {"order", 1},
                {"frequency", 60.0}
            }},
            {"materials", json::array()},
            {"sources", json::array()},
            {"boundaries", json::array()}
        };

        InputParser parser(config);
        mfem::Mesh mesh(mesh_file.c_str(), 1, 1);

        std::string type = parser.config["simulation"]["physics"];
        REQUIRE(type == "magnetoquasistatics");
    }

    // Cleanup
    fs::remove(mesh_file);
}

// ===========================================================================
// Adaptive Mesh Refinement (AMR) regression tests
// ===========================================================================

namespace {

// Build a 2D axisymmetric "coaxial capacitor" triangle mesh spanning
// r in [r_inner, r_outer], z in [0, height], and write it as an MFEM v1.0 mesh.
// Boundary attributes are assigned by geometry (independent of MFEM's internal
// edge numbering): 1 = inner conductor (r=r_inner), 2 = outer conductor
// (r=r_outer), 3 = top/bottom symmetry edges (z=0 or z=height). This exercises
// the axisymmetric ZienkiewiczZhu flux/energy path, and the analytic field
// E_r ~ 1/r concentrates error near the inner conductor so AMR has something to
// resolve.
void CreateCoaxMesh(const std::string& filename,
                    double r_inner, double r_outer, double height,
                    int nr, int nz) {
    const int nvr = nr + 1;
    const int nvz = nz + 1;
    auto vid = [nvr](int i, int j) { return j * nvr + i; }; // i along r, j along z

    std::vector<std::array<double, 2>> verts;
    verts.reserve(static_cast<size_t>(nvr) * nvz);
    for (int j = 0; j < nvz; ++j) {
        const double z = height * static_cast<double>(j) / nz;
        for (int i = 0; i < nvr; ++i) {
            const double r = r_inner + (r_outer - r_inner) * static_cast<double>(i) / nr;
            verts.push_back({ r, z });
        }
    }

    std::ofstream m(filename);
    m << "MFEM mesh v1.0\n\n";
    m << "dimension\n2\n\n";

    // Two triangles per cell, all domain attribute 1.
    m << "elements\n" << (2 * nr * nz) << "\n";
    for (int j = 0; j < nz; ++j) {
        for (int i = 0; i < nr; ++i) {
            const int v00 = vid(i, j);
            const int v10 = vid(i + 1, j);
            const int v11 = vid(i + 1, j + 1);
            const int v01 = vid(i, j + 1);
            m << "1 2 " << v00 << " " << v10 << " " << v11 << "\n";
            m << "1 2 " << v00 << " " << v11 << " " << v01 << "\n";
        }
    }
    m << "\n";

    // Boundary edges with geometry-based attributes.
    std::vector<std::array<int, 3>> bdr; // {attr, va, vb}
    for (int j = 0; j < nz; ++j) {       // vertical edges at r=r_inner / r=r_outer
        bdr.push_back({ 1, vid(0, j),  vid(0, j + 1) });
        bdr.push_back({ 2, vid(nr, j), vid(nr, j + 1) });
    }
    for (int i = 0; i < nr; ++i) {       // horizontal edges at z=0 / z=height
        bdr.push_back({ 3, vid(i, 0),  vid(i + 1, 0) });
        bdr.push_back({ 3, vid(i, nz), vid(i + 1, nz) });
    }

    m << "boundary\n" << bdr.size() << "\n";
    for (const auto& b : bdr) {
        m << b[0] << " 1 " << b[1] << " " << b[2] << "\n";
    }
    m << "\n";

    m << "vertices\n" << verts.size() << "\n2\n";
    for (const auto& v : verts) {
        m << v[0] << " " << v[1] << "\n";
    }
    m.close();
}

// Common AMR-enabled coaxial-capacitor config (Field analysis, one energized
// scenario). Outputs are off by default; callers enable Gmsh output as needed.
json MakeCoaxAmrConfig(const std::string& mesh_file, int max_iterations) {
    return json{
        {"simulation", {
            {"physics", "electrostatics"},
            {"mesh", mesh_file},
            {"order", 1},
            {"geometry_type", "axisymmetric"},
            {"analysis_type", "field"},
            {"amr", {
                {"enabled", true},
                {"max_iterations", max_iterations},
                {"error_fraction", 0.7},
                {"max_dofs", 0},          // disabled; iteration cap governs
                {"error_tolerance", 0.0}  // disabled; iteration cap governs
            }}
        }},
        {"entity_groups", json::array({
            {{"name", "Dielectric"}, {"dim", 2}, {"attribute_ids", {1}}},
            {{"name", "Inner"},      {"dim", 1}, {"attribute_ids", {1}}},
            {{"name", "Outer"},      {"dim", 1}, {"attribute_ids", {2}}}
        })},
        {"regions", json::array({
            {{"name", "Dielectric"}, {"entity_group", "Dielectric"}, {"material", 1}}
        })},
        {"materials", json::array({
            {{"name", "Vacuum"}, {"properties", {{"epsilon_r", 1.0}}}}
        })},
        {"terminals", json::array({
            {{"name", "Inner"}, {"excitation", "voltage"}, {"entity_group", "Inner"}},
            {{"name", "Outer"}, {"excitation", "voltage"}, {"entity_group", "Outer"}}
        })},
        {"scenarios", json::array({
            {{"name", "energized"}, {"excitations", json::array({
                {{"terminal", "Inner"}, {"value", 1.0}},
                {{"terminal", "Outer"}, {"value", 0.0}}
            })}}
        })}
    };
}

// Extract the raw text of a named MSH section, e.g. "$Nodes" ... "$EndNodes".
// Returns the lines strictly between the markers (exclusive), joined by '\n'.
std::string ExtractMshSection(const std::string& path, const std::string& tag) {
    std::ifstream f(path);
    REQUIRE(f.is_open());
    const std::string begin = "$" + tag;
    const std::string end = "$End" + tag;
    std::string line;
    bool in = false;
    std::ostringstream out;
    while (std::getline(f, line)) {
        if (!in) {
            if (line.rfind(begin, 0) == 0) { in = true; }
        } else {
            if (line.rfind(end, 0) == 0) { break; }
            out << line << '\n';
        }
    }
    return out.str();
}

} // namespace

TEST_CASE("AMR refines an axisymmetric coax and stays conforming", "[solvers][amr]") {
    const std::string mesh_file = "test_amr_coax.mesh";
    CreateCoaxMesh(mesh_file, /*r_inner=*/1.0, /*r_outer=*/4.0,
                   /*height=*/1.0, /*nr=*/6, /*nz=*/4);

    json config = MakeCoaxAmrConfig(mesh_file, /*max_iterations=*/4);

    InputParser parser(config);
    mfem::Mesh mesh(mesh_file.c_str(), 1, 1);

    ElectrostaticSolver solver(mesh, parser.config);
    solver.Setup();
    REQUIRE_NOTHROW(solver.Run());

    const auto& history = solver.GetAmrHistory();
    REQUIRE(history.size() >= 2);

    // (1) Conforming throughout: a simplex mesh refined conformingly never
    //     allocates an ncmesh; RefineConforming() also throws otherwise.
    REQUIRE(mesh.ncmesh == nullptr);
    REQUIRE_FALSE(mesh.Nonconforming());

    // (2) Work grows monotonically: each refined iteration adds true DOFs.
    for (std::size_t k = 1; k < history.size(); ++k) {
        REQUIRE(history[k].true_dofs > history[k - 1].true_dofs);
    }

    // (3) The recovery-based global error estimate trends down as we refine.
    //     Compare first vs last to avoid over-constraining intermediate steps.
    REQUIRE(history.back().global_error < history.front().global_error);

    // (4) Peak |E| stays finite/positive and physically plausible. With an
    //     analytic field E_r = V0 / (r ln(r_o/r_i)) the maximum is at the inner
    //     surface: 1 / (1 * ln(4)). Linear FEM recovers the peak from below and
    //     approaches it under refinement, so require the final peak to sit in a
    //     band around analytic - tight enough to catch a wrong axisymmetric
    //     measure, NaN, or sign error, loose enough not to depend on the exact
    //     (Dorfler + bisection) refinement depth.
    const double analytic_peak = 1.0 / (1.0 * std::log(4.0));
    const double last_peak = history.back().peak_absE;
    REQUIRE(std::isfinite(last_peak));
    REQUIRE(last_peak > 0.5 * analytic_peak);
    REQUIRE(last_peak < 1.2 * analytic_peak);

    fs::remove(mesh_file);
}

TEST_CASE("AMR with multiple scenarios writes a shared conforming mesh", "[solvers][amr]") {
    const std::string mesh_file = "test_amr_shared.mesh";
    CreateCoaxMesh(mesh_file, /*r_inner=*/1.0, /*r_outer=*/4.0,
                   /*height=*/1.0, /*nr=*/6, /*nz=*/4);

    // Unique temp output dir so the two results files land next to the mesh.
    const fs::path tmp_dir = fs::temp_directory_path() / "mfem_amr_shared_test";
    fs::create_directories(tmp_dir);
    const fs::path mesh_in_tmp = tmp_dir / mesh_file;
    fs::copy_file(mesh_file, mesh_in_tmp, fs::copy_options::overwrite_existing);

    json config = MakeCoaxAmrConfig(mesh_in_tmp.string(), /*max_iterations=*/3);
    config["simulation"]["output_gmsh"] = true;
    // Two scenarios; AMR builds ONE shared mesh and writes both results on it.
    config["scenarios"] = json::array({
        json{{"name", "driveA"}, {"excitations", json::array({
            json{{"terminal", "Inner"}, {"value", 1.0}},
            json{{"terminal", "Outer"}, {"value", 0.0}}
        })}},
        json{{"name", "driveB"}, {"excitations", json::array({
            json{{"terminal", "Inner"}, {"value", 0.0}},
            json{{"terminal", "Outer"}, {"value", 1.0}}
        })}}
    });

    InputParser parser(config);
    mfem::Mesh mesh(mesh_in_tmp.string().c_str(), 1, 1);

    ElectrostaticSolver solver(mesh, parser.config);
    solver.Setup();
    REQUIRE_NOTHROW(solver.Run());

    REQUIRE(mesh.ncmesh == nullptr);

    const fs::path resA = tmp_dir / "driveA.results.msh";
    const fs::path resB = tmp_dir / "driveB.results.msh";
    REQUIRE(fs::exists(resA));
    REQUIRE(fs::exists(resB));

    // Spec: every <scenario>.results.msh shares an identical mesh tiling. The
    // $Nodes and $Elements sections must be byte-identical across scenarios.
    REQUIRE(ExtractMshSection(resA.string(), "Nodes") ==
            ExtractMshSection(resB.string(), "Nodes"));
    REQUIRE(ExtractMshSection(resA.string(), "Elements") ==
            ExtractMshSection(resB.string(), "Elements"));

    fs::remove(mesh_file);
    std::error_code ec;
    fs::remove_all(tmp_dir, ec);
}
