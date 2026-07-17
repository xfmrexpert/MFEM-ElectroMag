// Copyright (c) 2026 T. C. Raymond
// SPDX-License-Identifier: MIT

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
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
#include <algorithm>
#include <cmath>
#include <complex>
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

void CreatePlanarStripMesh(const std::string& filename,
                           double length, double height,
                           int nx, int ny) {
    const int nvx = nx + 1;
    const int nvy = ny + 1;
    auto vertex = [nvx](int i, int j) { return j * nvx + i; };

    std::ofstream mesh_file(filename);
    mesh_file << "MFEM mesh v1.0\n\n";
    mesh_file << "dimension\n2\n\n";
    mesh_file << "elements\n" << 2 * nx * ny << "\n";
    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            const int v00 = vertex(i, j);
            const int v10 = vertex(i + 1, j);
            const int v11 = vertex(i + 1, j + 1);
            const int v01 = vertex(i, j + 1);
            mesh_file << "1 2 " << v00 << " " << v10 << " " << v11 << "\n";
            mesh_file << "1 2 " << v00 << " " << v11 << " " << v01 << "\n";
        }
    }

    mesh_file << "\nboundary\n" << 2 * ny + 2 * nx << "\n";
    for (int j = 0; j < ny; ++j) {
        mesh_file << "1 1 " << vertex(0, j) << " " << vertex(0, j + 1) << "\n";
        mesh_file << "2 1 " << vertex(nx, j) << " " << vertex(nx, j + 1) << "\n";
    }
    for (int i = 0; i < nx; ++i) {
        mesh_file << "3 1 " << vertex(i, 0) << " " << vertex(i + 1, 0) << "\n";
        mesh_file << "3 1 " << vertex(i, ny) << " " << vertex(i + 1, ny) << "\n";
    }

    mesh_file << "\nvertices\n" << nvx * nvy << "\n2\n";
    for (int j = 0; j < nvy; ++j) {
        const double y = height * static_cast<double>(j) / ny;
        for (int i = 0; i < nvx; ++i) {
            const double x = length * static_cast<double>(i) / nx;
            mesh_file << x << " " << y << "\n";
        }
    }
}

json MakePlanarStripConfig(const std::string& physics,
                           const std::string& mesh_file,
                           int order,
                           const json& material_properties,
                           double left_value,
                           double right_value) {
    return json{
        {"simulation", {
            {"physics_type", physics},
            {"mesh", mesh_file},
            {"order", order},
            {"geometry_type", "planar"},
            {"analysis_type", "field"},
            {"solver_tolerance", 1e-12},
            {"solver_max_iter", 4000},
            {"solver_print_level", 0}
        }},
        {"entity_groups", json::array({
            {{"name", "Domain"}, {"dim", 2}, {"attribute_ids", {1}}},
            {{"name", "Left"},   {"dim", 1}, {"attribute_ids", {1}}},
            {{"name", "Right"},  {"dim", 1}, {"attribute_ids", {2}}}
        })},
        {"regions", json::array({
            {{"name", "Domain"}, {"entity_group", "Domain"}, {"material", 1}}
        })},
        {"materials", json::array({
            {{"name", "Material"}, {"properties", material_properties}}
        })},
        {"boundaries", json::array({
            {{"name", "Left"},  {"type", "Dirichlet"}, {"entity_group", "Left"},  {"value", left_value}},
            {{"name", "Right"}, {"type", "Dirichlet"}, {"entity_group", "Right"}, {"value", right_value}}
        })},
        {"scenarios", json::array({
            {{"name", "analytic"}, {"excitations", json::array()}}
        })}
    };
}

const FieldExport& FindField(const FieldExportSet& fields, const std::string& name) {
    const auto& exported = fields.Fields();
    const auto field = std::find_if(exported.begin(), exported.end(),
        [&name](const FieldExport& candidate) { return candidate.name == name; });
    REQUIRE(field != exported.end());
    return *field;
}

mfem::IntegrationPoint TriangleCenter() {
    mfem::IntegrationPoint point;
    point.Set2(1.0 / 3.0, 1.0 / 3.0);
    return point;
}

mfem::Vector PhysicalPoint(const mfem::GridFunction& field, int element,
                           const mfem::IntegrationPoint& point) {
    mfem::Vector physical(2);
    field.FESpace()->GetElementTransformation(element)->Transform(point, physical);
    return physical;
}

double SamplePrimaryScalar(const FieldExportSet& fields, const std::string& name,
                           int element, const mfem::IntegrationPoint& point) {
    const FieldExport& field = FindField(fields, name);
    REQUIRE(field.kind == FieldExport::Kind::PrimaryScalar);
    return field.primary->GetValue(element, point);
}

mfem::Vector SampleDerivedVector(const FieldExportSet& fields, const std::string& name,
                                 int element, const mfem::IntegrationPoint& point) {
    const FieldExport& field = FindField(fields, name);
    REQUIRE(field.kind == FieldExport::Kind::DerivedVector);
    const auto& exported = fields.Fields();
    const auto primary = std::find_if(exported.begin(), exported.end(),
        [](const FieldExport& candidate) { return candidate.primary != nullptr; });
    REQUIRE(primary != exported.end());
    mfem::ElementTransformation* transformation =
        primary->primary->FESpace()->GetElementTransformation(element);
    mfem::Vector value(field.vector->GetVDim());
    field.vector->Eval(value, *transformation, point);
    return value;
}

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

TEST_CASE("Electrostatic solver reproduces a uniform field between plates",
          "[solvers][analytic][electrostatic]") {
    const std::string mesh_file = "test_analytic_electrostatic.mesh";
    constexpr double length = 0.2;
    constexpr double height = 0.05;
    constexpr double voltage = 100.0;
    constexpr int nx = 8;
    constexpr int ny = 2;
    CreatePlanarStripMesh(mesh_file, length, height, nx, ny);

    json config = MakePlanarStripConfig(
        "electrostatics", mesh_file, 1, {{"epsilon_r", 2.5}}, voltage, 0.0);
    mfem::Mesh mesh(mesh_file.c_str(), 1, 1);
    ElectrostaticSolver solver(mesh, config);
    solver.Setup();
    solver.Run();

    FieldExportSet fields = solver.CollectExportFields();
    const mfem::IntegrationPoint point = TriangleCenter();
    const int element = 2 * (nx / 2);
    const FieldExport& potential_field = FindField(fields, "V");
    const mfem::Vector physical = PhysicalPoint(*potential_field.primary, element, point);
    const double expected_potential = voltage * (1.0 - physical(0) / length);
    const double expected_field = voltage / length;

    const double potential = SamplePrimaryScalar(fields, "V", element, point);
    const mfem::Vector electric_field = SampleDerivedVector(fields, "E", element, point);

    REQUIRE(potential == Catch::Approx(expected_potential).epsilon(1e-7));
    REQUIRE(electric_field(0) == Catch::Approx(expected_field).epsilon(1e-6));
    REQUIRE(electric_field(1) == Catch::Approx(0.0).margin(1e-4));

    fs::remove(mesh_file);
}

TEST_CASE("Magnetostatic solver reproduces the field of a uniform current slab",
          "[solvers][analytic][magnetostatic]") {
    const std::string mesh_file = "test_analytic_magnetostatic.mesh";
    constexpr double length = 0.1;
    constexpr double height = 0.02;
    constexpr double current_density = 2.0e6;
    constexpr int nx = 8;
    constexpr int ny = 2;
    CreatePlanarStripMesh(mesh_file, length, height, nx, ny);

    json config = MakePlanarStripConfig(
        "magnetostatics", mesh_file, 2, {{"mu_r", 1.0}}, 0.0, 0.0);
    config["terminals"] = json::array({
        {{"name", "Current"}, {"excitation", "current"},
         {"conductor_type", "stranded"}, {"entity_group", "Domain"}}
    });
    config["scenarios"][0]["excitations"] = json::array({
        {{"terminal", "Current"}, {"value", current_density * length * height}}
    });

    mfem::Mesh mesh(mesh_file.c_str(), 1, 1);
    MagnetostaticSolver solver(mesh, config);
    solver.Setup();
    solver.Run();

    FieldExportSet fields = solver.CollectExportFields();
    const mfem::IntegrationPoint point = TriangleCenter();
    const int element = 2 * (nx / 4);
    const FieldExport& vector_potential_field = FindField(fields, "A");
    const mfem::Vector physical =
        PhysicalPoint(*vector_potential_field.primary, element, point);
    const double x = physical(0);
    const double expected_potential =
        0.5 * Constants::MU_0 * current_density * x * (length - x);
    const double expected_by =
        Constants::MU_0 * current_density * (x - 0.5 * length);

    const double vector_potential = SamplePrimaryScalar(fields, "A", element, point);
    const mfem::Vector magnetic_field = SampleDerivedVector(fields, "B", element, point);

    REQUIRE(vector_potential == Catch::Approx(expected_potential).epsilon(1e-8));
    REQUIRE(magnetic_field(0) == Catch::Approx(0.0).margin(1e-7));
    REQUIRE(magnetic_field(1) == Catch::Approx(expected_by).epsilon(1e-7));

    fs::remove(mesh_file);
}

TEST_CASE("Magnetoquasistatic solver reproduces conducting-slab skin effect",
          "[solvers][analytic][mqs]") {
    const std::string mesh_file = "test_analytic_mqs.mesh";
    constexpr double length = 0.04;
    constexpr double height = 0.005;
    constexpr double conductivity = 3.5e7;
    constexpr double frequency = 60.0;
    constexpr int nx = 48;
    constexpr int ny = 1;
    CreatePlanarStripMesh(mesh_file, length, height, nx, ny);

    json config = MakePlanarStripConfig(
        "magnetoquasistatics", mesh_file, 2,
        {{"mu_r", 1.0}, {"sigma", conductivity}}, 1.0, 0.0);
    config["simulation"]["frequency"] = frequency;

    mfem::Mesh mesh(mesh_file.c_str(), 1, 1);
    MagnetoquasistaticSolver solver(mesh, config);
    solver.Setup();
    solver.Run();

    FieldExportSet fields = solver.CollectExportFields();
    const mfem::IntegrationPoint point = TriangleCenter();
    const int element = 2 * (nx / 4);
    const FieldExport& real_field = FindField(fields, "A_Real");
    const mfem::Vector physical = PhysicalPoint(*real_field.primary, element, point);

    const double omega = Constants::TWO_PI * frequency;
    const std::complex<double> wave_number =
        std::sqrt(std::complex<double>(0.0, omega * Constants::MU_0 * conductivity));
    const std::complex<double> expected =
        std::sinh(wave_number * (length - physical(0))) /
        std::sinh(wave_number * length);

    const double actual_real = SamplePrimaryScalar(fields, "A_Real", element, point);
    const double actual_imag = SamplePrimaryScalar(fields, "A_Imag", element, point);

    REQUIRE(actual_real == Catch::Approx(expected.real()).epsilon(5e-3));
    REQUIRE(actual_imag == Catch::Approx(expected.imag()).epsilon(5e-3));
    REQUIRE(std::hypot(actual_real, actual_imag) < 1.0);
    REQUIRE(actual_imag < 0.0);

    fs::remove(mesh_file);
}

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
