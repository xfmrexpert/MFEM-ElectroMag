// Copyright (c) 2026 T. C. Raymond
// SPDX-License-Identifier: MIT

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include "../src/input_parser.hpp"
#include "../src/electrostatic_solver.hpp"
#include "../src/magnetostatic_solver.hpp"
#include "../src/magnetoquasistatic_solver.hpp"
#include "../src/solver_factory.hpp"
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
#include <functional>

namespace fs = std::filesystem;

ProblemConfig DecodeConfig(const json& config);

namespace {
void CreatePlanarStripMesh(const std::string& filename,
                           double length, double height, int nx, int ny);
json MakePlanarStripConfig(const std::string& physics,
                           const std::string& mesh_file,
                           int order,
                           const json& material_properties,
                           double left_value,
                           double right_value);
const FieldExport& FindField(const FieldExportSet& fields,
                             const std::string& name);
mfem::IntegrationPoint TriangleCenter();
mfem::Vector PhysicalPoint(const mfem::GridFunction& field, int element,
                           const mfem::IntegrationPoint& point);
double SamplePrimaryScalar(const FieldExportSet& fields,
                           const std::string& name, int element,
                           const mfem::IntegrationPoint& point);
mfem::Vector SampleDerivedVector(const FieldExportSet& fields,
                                 const std::string& name, int element,
                                 const mfem::IntegrationPoint& point);
} // namespace

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

TEST_CASE("Electrostatic solver applies natural-flux Neumann data",
          "[solvers][analytic][electrostatic][neumann]") {
    const std::string mesh_file = "test_electrostatic_neumann.mesh";
    constexpr double length = 0.2;
    constexpr double height = 0.05;
    constexpr double relative_permittivity = 2.5;
    constexpr double gradient = 3.0;
    constexpr int nx = 8;
    constexpr int ny = 2;
    CreatePlanarStripMesh(mesh_file, length, height, nx, ny);

    const double permittivity = Constants::EPSILON_0 * relative_permittivity;
    json config = MakePlanarStripConfig(
        "electrostatics", mesh_file, 1,
        {{"epsilon_r", relative_permittivity}}, 0.0, 0.0);
    config["boundaries"][1]["type"] = "Neumann";
    config["boundaries"][1]["value"] = permittivity * gradient;

    mfem::Mesh mesh(mesh_file.c_str(), 1, 1);
    ElectrostaticSolver solver(mesh, DecodeConfig(config));
    solver.Setup();
    solver.Run();

    const FieldExportSet fields = solver.CollectExportFields();
    const mfem::IntegrationPoint point = TriangleCenter();
    const int element = 2 * (nx / 2);
    const FieldExport& potential_field = FindField(fields, "V");
    const mfem::Vector physical =
        PhysicalPoint(*potential_field.primary, element, point);

    REQUIRE(SamplePrimaryScalar(fields, "V", element, point) ==
        Catch::Approx(gradient * physical(0)).epsilon(1.0e-7));
    const mfem::Vector electric_field =
        SampleDerivedVector(fields, "E", element, point);
    REQUIRE(electric_field(0) == Catch::Approx(-gradient).epsilon(1.0e-7));
    REQUIRE(electric_field(1) == Catch::Approx(0.0).margin(1.0e-6));

    fs::remove(mesh_file);
}

TEST_CASE("Boundary closures and voltage terminals have distinct ownership",
          "[solvers][boundaries][overlap]") {
    const std::string mesh_file = "test_boundary_ownership.mesh";
    CreatePlanarStripMesh(mesh_file, 0.2, 0.05, 2, 1);

    auto terminal_config = [&]() {
        json config = MakePlanarStripConfig(
            "electrostatics", mesh_file, 1, {{"epsilon_r", 1.0}}, 0.0, 0.0);
        config["terminals"] = json::array({
            {{"name", "LeftTerminal"}, {"excitation_type", "voltage"},
             {"entity_group", "Left"}}
        });
        config["boundaries"].erase(config["boundaries"].begin());
        return config;
    };

    SECTION("same boundary attribute is rejected") {
        json config = terminal_config();
        config["boundaries"].push_back(
            {{"name", "LeftFlux"}, {"type", "Neumann"},
             {"entity_group", "Left"}, {"value", 1.0}});

        mfem::Mesh mesh(mesh_file.c_str(), 1, 1);
        ElectrostaticSolver solver(mesh, DecodeConfig(config));
        REQUIRE_THROWS_WITH(solver.Setup(),
            Catch::Matchers::ContainsSubstring("one physical role"));
    }

    SECTION("different attributes may meet at a corner") {
        json config = terminal_config();
        config["entity_groups"].push_back(
            {{"name", "Horizontal"}, {"dim", 1}, {"attribute_ids", {3}}});
        config["boundaries"].push_back(
            {{"name", "HorizontalFlux"}, {"type", "Neumann"},
             {"entity_group", "Horizontal"}, {"value", 0.0}});

        mfem::Mesh mesh(mesh_file.c_str(), 1, 1);
        ElectrostaticSolver solver(mesh, DecodeConfig(config));
        REQUIRE_NOTHROW(solver.Setup());
    }

    fs::remove(mesh_file);
}

TEST_CASE("Axisymmetric magnetic solvers enforce zero A_phi on the axis",
          "[solvers][axisymmetric][axis][boundaries]") {
    const std::string mesh_file = "test_magnetic_axis_boundary.mesh";
    CreatePlanarStripMesh(mesh_file, 0.2, 0.05, 2, 1);

    auto magnetic_config = [&](const std::string& physics, double axis_value) {
        const json material = physics == "magnetoquasistatics"
            ? json{{"mu_r", 1.0}, {"sigma", 0.0}}
            : json{{"mu_r", 1.0}};
        json config = MakePlanarStripConfig(
            physics, mesh_file, 1, material, axis_value, 0.0);
        config["simulation"]["geometry_type"] = "axisymmetric";
        if (physics == "magnetoquasistatics") {
            config["scenarios"][0]["frequency"] = 60.0;
        }
        return config;
    };

    auto add_horizontal_dirichlet = [](json& config, double value) {
        config["entity_groups"].push_back(
            {{"name", "Horizontal"}, {"dim", 1}, {"attribute_ids", {3}}});
        config["boundaries"].push_back(
            {{"name", "Horizontal"}, {"type", "Dirichlet"},
             {"entity_group", "Horizontal"}, {"value", value}});
    };

    SECTION("magnetostatics rejects a nonzero axis value") {
        mfem::Mesh mesh(mesh_file.c_str(), 1, 1);
        MagnetostaticSolver solver(
            mesh, DecodeConfig(magnetic_config("magnetostatics", 1.0)));
        REQUIRE_THROWS_WITH(solver.Setup(),
            Catch::Matchers::ContainsSubstring("Axis regularity requires A_phi = 0"));
    }

    SECTION("magnetoquasistatics rejects a nonzero axis value") {
        mfem::Mesh mesh(mesh_file.c_str(), 1, 1);
        MagnetoquasistaticSolver solver(
            mesh, DecodeConfig(magnetic_config("magnetoquasistatics", 1.0)));
        REQUIRE_THROWS_WITH(solver.Setup(),
            Catch::Matchers::ContainsSubstring("Axis regularity requires A_phi = 0"));
    }

    SECTION("an explicit zero axis value is accepted") {
        mfem::Mesh mesh(mesh_file.c_str(), 1, 1);
        MagnetostaticSolver solver(
            mesh, DecodeConfig(magnetic_config("magnetostatics", 0.0)));
        REQUIRE_NOTHROW(solver.Setup());
    }

    SECTION("magnetostatics rejects a nonzero value at an axis corner") {
        json config = magnetic_config("magnetostatics", 0.0);
        add_horizontal_dirichlet(config, 1.0);
        mfem::Mesh mesh(mesh_file.c_str(), 1, 1);
        MagnetostaticSolver solver(mesh, DecodeConfig(config));
        REQUIRE_THROWS_WITH(solver.Setup(),
            Catch::Matchers::ContainsSubstring("on the magnetic symmetry axis"));
    }

    SECTION("magnetoquasistatics rejects a nonzero value at an axis corner") {
        json config = magnetic_config("magnetoquasistatics", 0.0);
        add_horizontal_dirichlet(config, 1.0);
        mfem::Mesh mesh(mesh_file.c_str(), 1, 1);
        MagnetoquasistaticSolver solver(mesh, DecodeConfig(config));
        REQUIRE_THROWS_WITH(solver.Setup(),
            Catch::Matchers::ContainsSubstring("on the magnetic symmetry axis"));
    }

    SECTION("a zero value on a boundary adjacent to the axis is accepted") {
        json config = magnetic_config("magnetostatics", 0.0);
        add_horizontal_dirichlet(config, 0.0);
        mfem::Mesh mesh(mesh_file.c_str(), 1, 1);
        MagnetostaticSolver solver(mesh, DecodeConfig(config));
        REQUIRE_NOTHROW(solver.Setup());
    }

    SECTION("planar problems do not apply magnetic axis regularity") {
        json config = magnetic_config("magnetostatics", 1.0);
        config["simulation"]["geometry_type"] = "planar";
        mfem::Mesh mesh(mesh_file.c_str(), 1, 1);
        MagnetostaticSolver solver(mesh, DecodeConfig(config));
        REQUIRE_NOTHROW(solver.Setup());
    }

    SECTION("annular domains may use a nonzero inner-boundary value") {
        mfem::Mesh mesh(mesh_file.c_str(), 1, 1);
        for (int v = 0; v < mesh.GetNV(); ++v) {
            mesh.GetVertex(v)[0] += 1.0;
        }
        MagnetostaticSolver solver(
            mesh, DecodeConfig(magnetic_config("magnetostatics", 1.0)));
        REQUIRE_NOTHROW(solver.Setup());
    }

    fs::remove(mesh_file);
}

ProblemConfig DecodeConfig(const json& config) {
    return InputParser(config).GetProblemConfig();
}

TEST_CASE("Magnetoquasistatic Neumann data loads only the real field",
          "[solvers][analytic][mqs][neumann]") {
    const std::string mesh_file = "test_mqs_neumann.mesh";
    constexpr double length = 0.1;
    constexpr double height = 0.02;
    constexpr double gradient = 0.75;
    constexpr int nx = 8;
    constexpr int ny = 2;
    CreatePlanarStripMesh(mesh_file, length, height, nx, ny);

    const double reluctivity = 1.0 / Constants::MU_0;
    json config = MakePlanarStripConfig(
        "magnetoquasistatics", mesh_file, 1,
        {{"mu_r", 1.0}, {"sigma", 0.0}}, 0.0, 0.0);
    config["boundaries"][1]["type"] = "Neumann";
    config["boundaries"][1]["value"] = reluctivity * gradient;
    config["scenarios"][0]["frequency"] = 60.0;

    mfem::Mesh mesh(mesh_file.c_str(), 1, 1);
    MagnetoquasistaticSolver solver(mesh, DecodeConfig(config));
    solver.Setup();
    solver.Run();

    const FieldExportSet fields = solver.CollectExportFields();
    const mfem::IntegrationPoint point = TriangleCenter();
    const int element = 2 * (nx / 2);
    const FieldExport& real_field = FindField(fields, "A_Real");
    const mfem::Vector physical = PhysicalPoint(*real_field.primary, element, point);

    REQUIRE(SamplePrimaryScalar(fields, "A_Real", element, point) ==
        Catch::Approx(gradient * physical(0)).epsilon(1.0e-8));
    REQUIRE(SamplePrimaryScalar(fields, "A_Imag", element, point) ==
        Catch::Approx(0.0).margin(1.0e-10));

    fs::remove(mesh_file);
}

TEST_CASE("Solvers reject reserved Robin boundary conditions during setup",
          "[solvers][boundaries][robin]") {
    const std::string mesh_file = "test_robin_rejection.mesh";
    CreatePlanarStripMesh(mesh_file, 0.1, 0.02, 2, 1);

    auto robin_config = [&](const std::string& physics, const json& material) {
        json config = MakePlanarStripConfig(
            physics, mesh_file, 1, material, 0.0, 0.0);
        config["boundaries"][1]["type"] = "Robin";
        config["boundaries"][1]["robin_coefficient"] = 1.0;
        if (physics == "magnetoquasistatics") {
            config["scenarios"][0]["frequency"] = 60.0;
        }
        return config;
    };

    SECTION("electrostatics") {
        mfem::Mesh mesh(mesh_file.c_str(), 1, 1);
        ElectrostaticSolver solver(mesh, DecodeConfig(
            robin_config("electrostatics", {{"epsilon_r", 1.0}})));
        REQUIRE_THROWS_WITH(solver.Setup(),
            Catch::Matchers::ContainsSubstring("reserved but not implemented"));
    }

    SECTION("magnetostatics") {
        mfem::Mesh mesh(mesh_file.c_str(), 1, 1);
        MagnetostaticSolver solver(mesh, DecodeConfig(
            robin_config("magnetostatics", {{"mu_r", 1.0}})));
        REQUIRE_THROWS_WITH(solver.Setup(),
            Catch::Matchers::ContainsSubstring("reserved but not implemented"));
    }

    SECTION("magnetoquasistatics") {
        mfem::Mesh mesh(mesh_file.c_str(), 1, 1);
        MagnetoquasistaticSolver solver(mesh, DecodeConfig(
            robin_config("magnetoquasistatics",
                         {{"mu_r", 1.0}, {"sigma", 0.0}})));
        REQUIRE_THROWS_WITH(solver.Setup(),
            Catch::Matchers::ContainsSubstring("reserved but not implemented"));
    }

    fs::remove(mesh_file);
}

class ScenarioOrderProbe : public PhysicsSolver {
public:
    ScenarioOrderProbe(mfem::Mesh& mesh, const ProblemConfig& config)
        : PhysicsSolver(mesh, config) {}

    std::vector<std::string> ScenarioNames() const {
        std::vector<std::string> names;
        for (const auto& [name, scenario] : BuildSolveScenarios()) {
            names.push_back(name);
        }
        return names;
    }

    void Setup() override {}
    void SaveAnalysis() override {}
    FieldExportSet CollectExportFields() const override { return {}; }

protected:
    void BuildOperators() override {}
    void RunOnCurrentMesh() override {}
    void EstimateCurrentSolutionError(mfem::Vector& errors) override {
        errors.SetSize(mesh.GetNE());
        errors = 0.0;
    }
    double ComputePeakFieldMagnitude() const override { return 0.0; }
};

class AmrLifecycleProbe : public PhysicsSolver {
public:
    AmrLifecycleProbe(mfem::Mesh& mesh, const ProblemConfig& config)
        : PhysicsSolver(mesh, config) {}

    void Setup() override {
        fec = std::make_unique<mfem::H1_FECollection>(1, mesh.Dimension());
        BuildOperators();
    }

    void SaveAnalysis() override {}
    FieldExportSet CollectExportFields() const override { return {}; }

    int OperatorBuilds() const { return operator_builds; }
    int ErrorEstimates() const { return error_estimates; }
    int SolvePasses() const { return solve_passes; }
    const std::vector<amr::AmrIterationInfo>& AmrHistory() const {
        return GetAmrHistory();
    }

protected:
    void BuildOperators() override {
        ++operator_builds;
        fespace = std::make_unique<mfem::FiniteElementSpace>(&mesh, fec.get());
    }

    // One scenario per pass, folding its error indicator in as production
    // solvers do.
    void RunOnCurrentMesh() override {
        ++solve_passes;
        AccumulateScenarioError();
    }

    void EstimateCurrentSolutionError(mfem::Vector& errors) override {
        ++error_estimates;
        errors.SetSize(mesh.GetNE());
        errors = 1.0;
    }

    double ComputePeakFieldMagnitude() const override { return 2.5; }

private:
    int operator_builds = 0;
    int error_estimates = 0;
    int solve_passes = 0;
};

class ElectrostaticAmrProbe : public ElectrostaticSolver {
public:
    using ElectrostaticSolver::ElectrostaticSolver;

    const std::vector<amr::AmrIterationInfo>& AmrHistory() const {
        return GetAmrHistory();
    }

    mfem::Vector LocalErrors() {
        mfem::Vector errors;
        EstimateCurrentSolutionError(errors);
        return errors;
    }
};

class MagnetostaticAmrProbe : public MagnetostaticSolver {
public:
    using MagnetostaticSolver::MagnetostaticSolver;

    mfem::Vector LocalErrors() {
        mfem::Vector errors;
        EstimateCurrentSolutionError(errors);
        return errors;
    }
};

TEST_CASE("ElectrostaticSolver can be constructed", "[solvers]") {
    std::string mesh_file = "test_electrostatic.mesh";
    CreateTestMesh(mesh_file);

    json test_config = {
        {"simulation", {
            {"physics_type", "electrostatics"},
            {"mesh", mesh_file},
            {"order", 1},
            {"geometry_type", "axisymmetric"}
        }},
        {"entity_groups", json::array({
            {{"name", "Domain"}, {"dim", 2}, {"attribute_ids", {1}}},
            {{"name", "Ground"}, {"dim", 1}, {"attribute_ids", {1}}}
        })},
        {"regions", json::array({
            {{"entity_group", "Domain"}, {"material", "dielectric"}}
        })},
        {"materials", json::array({
            {
                {"name", "dielectric"},
                {"properties", {{"epsilon_r", 2.0}}}
            }
        })},
        {"boundaries", json::array({
            {
                {"name", "ground"},
                {"entity_group", "Ground"},
                {"type", "Dirichlet"},
                {"value", 0.0}
            }
        })}
    };

    InputParser parser(test_config);
    ProblemConfig config = parser.GetProblemConfig();
    mfem::Mesh mesh(mesh_file.c_str(), 1, 1);

    REQUIRE_NOTHROW(ElectrostaticSolver(mesh, config));

    // Cleanup
    fs::remove(mesh_file);
}

TEST_CASE("PhysicsSolver owns the adaptive lifecycle", "[solvers][amr]") {
    const std::string mesh_file = "test_amr_lifecycle.mesh";
    CreateTestMesh(mesh_file);
    mfem::Mesh mesh(mesh_file.c_str(), 1, 1);

    ProblemConfig config;
    config.Amr.Enabled = true;
    config.Amr.MaxIterations = 3;
    config.Amr.MaxDofs = 0;
    config.Amr.ErrorFraction = 1.0;

    AmrLifecycleProbe probe(mesh, config);
    probe.Setup();
    probe.Run();

    REQUIRE(probe.OperatorBuilds() == 3);
    REQUIRE(probe.ErrorEstimates() == 3);
    // The AMR loop now runs exactly one solve pass per iteration; the last one
    // is on the converged mesh, so there is no extra pass afterwards.
    REQUIRE(probe.SolvePasses() == 3);
    REQUIRE(probe.AmrHistory().size() == 3);
    REQUIRE(probe.AmrHistory().back().peak_field_magnitude == 2.5);
    REQUIRE(mesh.GetNE() > 1);
    REQUIRE_FALSE(mesh.Nonconforming());

    fs::remove(mesh_file);
}

TEST_CASE("Field scenarios preserve declaration order", "[solvers][scenarios]") {
    const std::string mesh_file = "test_scenario_order.mesh";
    CreateTestMesh(mesh_file);
    mfem::Mesh mesh(mesh_file.c_str(), 1, 1);

    ProblemConfig config;
    config.AnalysisType = AnalysisType::Field;
    config.Scenarios.emplace_back("Zulu", Scenario{});
    config.Scenarios.emplace_back("Alpha", Scenario{});
    config.Scenarios.emplace_back("Middle", Scenario{});

    ScenarioOrderProbe probe(mesh, config);
    REQUIRE((probe.ScenarioNames() == std::vector<std::string>{"Zulu", "Alpha", "Middle"}));

    fs::remove(mesh_file);
}

TEST_CASE("MagnetostaticSolver can be constructed", "[solvers]") {
    std::string mesh_file = "test_magnetostatic.mesh";
    CreateTestMesh(mesh_file);

    json test_config = {
        {"simulation", {
            {"physics_type", "magnetostatics"},
            {"mesh", mesh_file},
            {"order", 1},
            {"geometry_type", "axisymmetric"}
        }},
        {"entity_groups", json::array({
            {{"name", "Domain"}, {"dim", 2}, {"attribute_ids", {1}}},
            {{"name", "FarField"}, {"dim", 1}, {"attribute_ids", {1}}}
        })},
        {"regions", json::array({
            {{"entity_group", "Domain"}, {"material", "iron"}}
        })},
        {"materials", json::array({
            {
                {"name", "iron"},
                {"properties", {{"mu_r", 1000.0}}}
            }
        })},
        {"boundaries", json::array({
            {
                {"name", "far_field"},
                {"entity_group", "FarField"},
                {"type", "Dirichlet"},
                {"value", 0.0}
            }
        })}
    };

    InputParser parser(test_config);
    ProblemConfig config = parser.GetProblemConfig();
    mfem::Mesh mesh(mesh_file.c_str(), 1, 1);

    REQUIRE_NOTHROW(MagnetostaticSolver(mesh, config));

    // Cleanup
    fs::remove(mesh_file);
}

TEST_CASE("MagnetoquasistaticSolver can be constructed", "[solvers]") {
    std::string mesh_file = "test_mqs.mesh";
    CreateTestMesh(mesh_file);

    json test_config = {
        {"simulation", {
            {"physics_type", "magnetoquasistatics"},
            {"mesh", mesh_file},
            {"order", 1},
            {"geometry_type", "axisymmetric"}
        }},
        {"entity_groups", json::array({
            {{"name", "Domain"}, {"dim", 2}, {"attribute_ids", {1}}},
            {{"name", "FarField"}, {"dim", 1}, {"attribute_ids", {1}}}
        })},
        {"regions", json::array({
            {{"entity_group", "Domain"}, {"material", "conductor"}}
        })},
        {"materials", json::array({
            {
                {"name", "conductor"},
                {"properties", {
                    {"mu_r", 1.0},
                    {"sigma", 5.8e7}
                }}
            }
        })},
        {"boundaries", json::array({
            {
                {"name", "far_field"},
                {"entity_group", "FarField"},
                {"type", "Dirichlet"},
                {"value", 0.0}
            }
        })},
        {"scenarios", json::array({
            {{"name", "default"}, {"frequency", 60.0}, {"excitations", json::array()}}
        })}
    };

    InputParser parser(test_config);
    ProblemConfig config = parser.GetProblemConfig();
    mfem::Mesh mesh(mesh_file.c_str(), 1, 1);

    REQUIRE_NOTHROW(MagnetoquasistaticSolver(mesh, config));

    // Cleanup
    fs::remove(mesh_file);
}

TEST_CASE("Solver factory logic works correctly", "[solvers]") {
    std::string mesh_file = "test_factory.mesh";
    CreateTestMesh(mesh_file);
    mfem::Mesh mesh(mesh_file.c_str(), 1, 1);

    SECTION("electrostatics") {
        ProblemConfig config;
        config.PhysicsType = PhysicsType::Electrostatics;
        auto solver = SolverFactory::Instance().Create(mesh, config);
        REQUIRE(dynamic_cast<ElectrostaticSolver*>(solver.get()) != nullptr);
    }

    SECTION("magnetostatics") {
        ProblemConfig config;
        config.PhysicsType = PhysicsType::Magnetostatics;
        auto solver = SolverFactory::Instance().Create(mesh, config);
        REQUIRE(dynamic_cast<MagnetostaticSolver*>(solver.get()) != nullptr);
    }

    SECTION("magnetoquasistatics") {
        ProblemConfig config;
        config.PhysicsType = PhysicsType::Magnetoquasistatics;
        auto solver = SolverFactory::Instance().Create(mesh, config);
        REQUIRE(dynamic_cast<MagnetoquasistaticSolver*>(solver.get()) != nullptr);
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

void CreateLayeredStripMesh(const std::string& filename,
                            double length, double height,
                            int nx, int ny, int interface_column) {
    REQUIRE(interface_column > 0);
    REQUIRE(interface_column < nx);

    const int nvx = nx + 1;
    const int nvy = ny + 1;
    auto vertex = [nvx](int i, int j) { return j * nvx + i; };

    std::ofstream mesh_file(filename);
    mesh_file << "MFEM mesh v1.0\n\n";
    mesh_file << "dimension\n2\n\n";
    mesh_file << "elements\n" << 2 * nx * ny << "\n";
    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            const int attribute = i < interface_column ? 1 : 2;
            const int v00 = vertex(i, j);
            const int v10 = vertex(i + 1, j);
            const int v11 = vertex(i + 1, j + 1);
            const int v01 = vertex(i, j + 1);
            mesh_file << attribute << " 2 " << v00 << " " << v10 << " " << v11 << "\n";
            mesh_file << attribute << " 2 " << v00 << " " << v11 << " " << v01 << "\n";
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
            {{"name", "Domain"}, {"entity_group", "Domain"}, {"material", "Material"}}
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

double IntegrateVectorMagnitudeSquared(const FieldExportSet& fields,
                                       const std::string& name) {
    const FieldExport& field = FindField(fields, name);
    REQUIRE(field.kind == FieldExport::Kind::DerivedVector);
    const auto& exported = fields.Fields();
    const auto primary = std::find_if(exported.begin(), exported.end(),
        [](const FieldExport& candidate) { return candidate.primary != nullptr; });
    REQUIRE(primary != exported.end());

    const mfem::FiniteElementSpace* space = primary->primary->FESpace();
    double integral = 0.0;
    for (int element = 0; element < space->GetNE(); ++element) {
        const mfem::FiniteElement* finite_element = space->GetFE(element);
        mfem::ElementTransformation* transformation =
            space->GetElementTransformation(element);
        const mfem::IntegrationRule& rule = mfem::IntRules.Get(
            finite_element->GetGeomType(), 2 * finite_element->GetOrder() + 2);
        for (int q = 0; q < rule.GetNPoints(); ++q) {
            const mfem::IntegrationPoint& point = rule.IntPoint(q);
            transformation->SetIntPoint(&point);
            mfem::Vector value(field.vector->GetVDim());
            field.vector->Eval(value, *transformation, point);
            integral += point.weight * transformation->Weight() * (value * value);
        }
    }
    return integral;
}

double RelativeComplexL2Error(
    const FieldExportSet& fields,
    const std::function<std::complex<double>(const mfem::Vector&)>& exact) {
    const FieldExport& real_field = FindField(fields, "A_Real");
    const FieldExport& imag_field = FindField(fields, "A_Imag");
    REQUIRE(real_field.kind == FieldExport::Kind::PrimaryScalar);
    REQUIRE(imag_field.kind == FieldExport::Kind::PrimaryScalar);

    const mfem::FiniteElementSpace* space = real_field.primary->FESpace();
    double error_squared = 0.0;
    double exact_squared = 0.0;
    for (int element = 0; element < space->GetNE(); ++element) {
        const mfem::FiniteElement* finite_element = space->GetFE(element);
        mfem::ElementTransformation* transformation =
            space->GetElementTransformation(element);
        const mfem::IntegrationRule& rule = mfem::IntRules.Get(
            finite_element->GetGeomType(), 2 * finite_element->GetOrder() + 4);
        for (int q = 0; q < rule.GetNPoints(); ++q) {
            const mfem::IntegrationPoint& point = rule.IntPoint(q);
            transformation->SetIntPoint(&point);
            mfem::Vector physical(2);
            transformation->Transform(point, physical);
            const std::complex<double> expected = exact(physical);
            const std::complex<double> actual(
                real_field.primary->GetValue(element, point),
                imag_field.primary->GetValue(element, point));
            const double weight = point.weight * transformation->Weight();
            error_squared += weight * std::norm(actual - expected);
            exact_squared += weight * std::norm(expected);
        }
    }
    return std::sqrt(error_squared / exact_squared);
}

struct CsvMatrix {
    std::vector<std::string> labels;
    std::vector<std::vector<double>> values;
};

CsvMatrix ReadCsvMatrix(const std::string& filename) {
    std::ifstream input(filename);
    REQUIRE(input.is_open());

    CsvMatrix matrix;
    std::string line;
    REQUIRE(static_cast<bool>(std::getline(input, line)));
    std::istringstream header(line);
    std::string cell;
    REQUIRE(static_cast<bool>(std::getline(header, cell, ',')));
    while (std::getline(header, cell, ',')) {
        matrix.labels.push_back(cell);
    }

    while (std::getline(input, line)) {
        std::istringstream row(line);
        REQUIRE(static_cast<bool>(std::getline(row, cell, ',')));
        std::vector<double> values;
        while (std::getline(row, cell, ',')) {
            values.push_back(std::stod(cell));
        }
        REQUIRE(values.size() == matrix.labels.size());
        matrix.values.push_back(std::move(values));
    }
    REQUIRE(matrix.values.size() == matrix.labels.size());
    return matrix;
}

void RequireMatricesEqual(const CsvMatrix& actual, const CsvMatrix& expected,
                          double tolerance = 1e-10) {
    REQUIRE(actual.labels == expected.labels);
    REQUIRE(actual.values.size() == expected.values.size());
    for (std::size_t row = 0; row < actual.values.size(); ++row) {
        REQUIRE(actual.values[row].size() == expected.values[row].size());
        for (std::size_t column = 0; column < actual.values[row].size(); ++column) {
            REQUIRE(actual.values[row][column] ==
                Catch::Approx(expected.values[row][column]).epsilon(tolerance));
        }
    }
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
                    int nr, int nz, int interface_column = 0) {
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

    // Two triangles per cell. With a positive interface column, radial cells to
    // its right use domain attribute 2; otherwise the whole domain is attribute 1.
    m << "elements\n" << (2 * nr * nz) << "\n";
    for (int j = 0; j < nz; ++j) {
        for (int i = 0; i < nr; ++i) {
            const int attribute = interface_column > 0 && i >= interface_column ? 2 : 1;
            const int v00 = vid(i, j);
            const int v10 = vid(i + 1, j);
            const int v11 = vid(i + 1, j + 1);
            const int v01 = vid(i, j + 1);
            m << attribute << " 2 " << v00 << " " << v10 << " " << v11 << "\n";
            m << attribute << " 2 " << v00 << " " << v11 << " " << v01 << "\n";
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

// Build a 2D axisymmetric shield-and-turns mesh: an aluminum shield band
// nearest the axis, then two rectangular conductor turns further out in r,
// each separated by air. All conductors are held off the symmetry axis so the
// massive-port conductance integral sigma/(2*pi*r) stays finite.
//
// Radial layout (r increasing), with the domain attribute in parentheses:
//   air (1) | shield (2) | air (1) | turn A (3) | air (1) | turn B (4) | air (1)
//
// Boundary attributes: 1 = inner (r=r_min), 2 = outer (r=r_max),
// 3 = top/bottom (z=0 or z=height).
void CreateShieldedTurnsMesh(const std::string& filename,
                             double r_min, double r_max, double height,
                             int nz = 4, int cells_per_band = 2) {
    // Radial band edges as fractions of the span, one entry per band, each
    // paired with the domain attribute carried by that band.
    struct Band { int cells; int attribute; };
    const std::vector<Band> bands{
        { cells_per_band, 1 },  // air gap next to the axis side
        { cells_per_band, 2 },  // aluminum shield
        { cells_per_band, 1 },  // air between shield and turns
        { cells_per_band, 3 },  // turn A
        { cells_per_band, 1 },  // air between turns
        { cells_per_band, 4 },  // turn B
        { cells_per_band, 1 }   // air out to the far boundary
    };

    int nr = 0;
    for (const Band& band : bands) nr += band.cells;
    const int nvr = nr + 1;
    const int nvz = nz + 1;
    auto vid = [nvr](int i, int j) { return j * nvr + i; };

    // Per-radial-cell attribute, expanded from the band table.
    std::vector<int> cell_attribute;
    cell_attribute.reserve(nr);
    for (const Band& band : bands) {
        for (int c = 0; c < band.cells; ++c) cell_attribute.push_back(band.attribute);
    }

    std::ofstream m(filename);
    m << "MFEM mesh v1.0\n\n";
    m << "dimension\n2\n\n";

    m << "elements\n" << (2 * nr * nz) << "\n";
    for (int j = 0; j < nz; ++j) {
        for (int i = 0; i < nr; ++i) {
            const int attribute = cell_attribute[i];
            const int v00 = vid(i, j);
            const int v10 = vid(i + 1, j);
            const int v11 = vid(i + 1, j + 1);
            const int v01 = vid(i, j + 1);
            m << attribute << " 2 " << v00 << " " << v10 << " " << v11 << "\n";
            m << attribute << " 2 " << v00 << " " << v11 << " " << v01 << "\n";
        }
    }
    m << "\n";

    std::vector<std::array<int, 3>> bdr; // {attr, va, vb}
    for (int j = 0; j < nz; ++j) {
        bdr.push_back({ 1, vid(0, j),  vid(0, j + 1) });
        bdr.push_back({ 2, vid(nr, j), vid(nr, j + 1) });
    }
    for (int i = 0; i < nr; ++i) {
        bdr.push_back({ 3, vid(i, 0),  vid(i + 1, 0) });
        bdr.push_back({ 3, vid(i, nz), vid(i + 1, nz) });
    }

    m << "boundary\n" << bdr.size() << "\n";
    for (const auto& b : bdr) {
        m << b[0] << " 1 " << b[1] << " " << b[2] << "\n";
    }
    m << "\n";

    m << "vertices\n" << (nvr * nvz) << "\n2\n";
    for (int j = 0; j < nvz; ++j) {
        const double z = height * static_cast<double>(j) / nz;
        for (int i = 0; i < nvr; ++i) {
            const double r = r_min + (r_max - r_min) * static_cast<double>(i) / nr;
            m << r << " " << z << "\n";
        }
    }
    m.close();
}

// Coupling-matrix config for the shield-and-turns mesh at a single frequency
// point. Region index 1 is the shield, so callers switch its material or add a
// current constraint to select between the shielding cases under test.
json MakeShieldedTurnsConfig(const std::string& mesh_file, double frequency) {
    return json{
        {"simulation", {
            {"physics_type", "magnetoquasistatics"},
            {"mesh", mesh_file},
            {"order", 1},
            {"geometry_type", "axisymmetric"},
            {"analysis_type", "coupling_matrix"},
            {"solver_tolerance", 1e-12},
            {"solver_max_iter", 4000},
            {"solver_print_level", 0}
        }},
        {"entity_groups", json::array({
            {{"name", "AirDomain"},    {"dim", 2}, {"attribute_ids", {1}}},
            {{"name", "ShieldDomain"}, {"dim", 2}, {"attribute_ids", {2}}},
            {{"name", "TurnADomain"},  {"dim", 2}, {"attribute_ids", {3}}},
            {{"name", "TurnBDomain"},  {"dim", 2}, {"attribute_ids", {4}}},
            {{"name", "Outer"},        {"dim", 1}, {"attribute_ids", {2}}}
        })},
        {"regions", json::array({
            {{"name", "Air"},    {"entity_group", "AirDomain"},    {"material", "Air"}},
            {{"name", "Shield"}, {"entity_group", "ShieldDomain"}, {"material", "Aluminum"}},
            {{"name", "TurnA"},  {"entity_group", "TurnADomain"},  {"material", "Copper"}},
            {{"name", "TurnB"},  {"entity_group", "TurnBDomain"},  {"material", "Copper"}}
        })},
        {"materials", json::array({
            {{"name", "Air"},      {"properties", {{"mu_r", 1.0}, {"sigma", 0.0}}}},
            {{"name", "Aluminum"}, {"properties", {{"mu_r", 1.0}, {"sigma", 3.5e7}}}},
            {{"name", "Copper"},   {"properties", {{"mu_r", 1.0}, {"sigma", 5.8e7}}}}
        })},
        {"boundaries", json::array({
            {{"name", "Outer"}, {"type", "Dirichlet"}, {"entity_group", "Outer"}, {"value", 0.0}}
        })},
        {"terminals", json::array({
            {{"name", "TurnA"}, {"excitation_type", "current"},
             {"conductor_type", "massive"}, {"entity_group", "TurnADomain"}},
            {{"name", "TurnB"}, {"excitation_type", "current"},
             {"conductor_type", "massive"}, {"entity_group", "TurnBDomain"}}
        })},
        {"scenarios", json::array({
            {{"name", "point"},
             {"frequency", {{"scale", "linear"}, {"start", frequency},
                            {"stop", frequency}, {"points", 1}}},
             {"excitations", json::array()}}
        })}
    };
}

// Common AMR-enabled coaxial-capacitor config (Field analysis, one energized
// scenario). Outputs are off by default; callers enable Gmsh output as needed.
json MakeCoaxAmrConfig(const std::string& mesh_file, int max_iterations) {
    return json{
        {"simulation", {
            {"physics_type", "electrostatics"},
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
            {{"name", "Dielectric"}, {"entity_group", "Dielectric"}, {"material", "Vacuum"}}
        })},
        {"materials", json::array({
            {{"name", "Vacuum"}, {"properties", {{"epsilon_r", 1.0}}}}
        })},
        {"terminals", json::array({
            {{"name", "Inner"}, {"excitation_type", "voltage"}, {"entity_group", "Inner"}},
            {{"name", "Outer"}, {"excitation_type", "voltage"}, {"entity_group", "Outer"}}
        })},
        {"scenarios", json::array({
            {{"name", "energized"}, {"excitations", json::array({
                {{"terminal", "Inner"}, {"value", 1.0}},
                {{"terminal", "Outer"}, {"value", 0.0}}
            })}}
        })}
    };
}

mfem::Vector AxisymmetricDiffusionErrors(
    mfem::Mesh& mesh, mfem::GridFunction& solution,
    mfem::Coefficient& coefficient, bool with_coefficient) {
    mfem::H1_FECollection flux_fec(1, mesh.Dimension());
    mfem::FiniteElementSpace flux_fes(
        &mesh, &flux_fec, mesh.SpaceDimension());
    AxisymmetricDiffusionIntegrator integrator(coefficient);
    mfem::ZienkiewiczZhuEstimator estimator(
        integrator, solution, flux_fes);
    estimator.SetWithCoeff(with_coefficient);
    estimator.SetFluxAveraging(1);
    return estimator.GetLocalErrors();
}

mfem::Vector AxisymmetricCurlCurlErrors(
    mfem::Mesh& mesh, mfem::GridFunction& solution,
    mfem::Coefficient& coefficient, bool with_coefficient) {
    mfem::H1_FECollection flux_fec(1, mesh.Dimension());
    mfem::FiniteElementSpace flux_fes(
        &mesh, &flux_fec, mesh.SpaceDimension());
    AxisymmetricCurlCurlIntegrator integrator(
        coefficient, axisym::InspectMesh(mesh).tolerance);
    mfem::ZienkiewiczZhuEstimator estimator(
        integrator, solution, flux_fes);
    estimator.SetWithCoeff(with_coefficient);
    estimator.SetFluxAveraging(1);
    return estimator.GetLocalErrors();
}

void RequireErrorsEqual(const mfem::Vector& actual,
                        const mfem::Vector& expected) {
    REQUIRE(actual.Size() == expected.Size());
    for (int i = 0; i < actual.Size(); ++i) {
        REQUIRE(actual(i) == Catch::Approx(expected(i)).epsilon(1e-12));
    }
}

bool ErrorsDiffer(const mfem::Vector& first, const mfem::Vector& second) {
    if (first.Size() != second.Size()) return true;
    for (int i = 0; i < first.Size(); ++i) {
        const double scale = std::max(std::abs(first(i)), std::abs(second(i)));
        if (std::abs(first(i) - second(i)) > 1e-8 * scale) return true;
    }
    return false;
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
    ElectrostaticSolver solver(mesh, DecodeConfig(config));
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

TEST_CASE("Electrostatic field energy gives analytic capacitance and conserves flux",
          "[solvers][analytic][electrostatic][energy][conservation]") {
    const std::string mesh_file = "test_electrostatic_energy.mesh";
    constexpr double length = 0.2;
    constexpr double height = 0.05;
    constexpr double voltage = 100.0;
    constexpr double relative_permittivity = 2.5;
    constexpr int nx = 8;
    constexpr int ny = 3;
    CreatePlanarStripMesh(mesh_file, length, height, nx, ny);

    json config = MakePlanarStripConfig(
        "electrostatics", mesh_file, 1,
        {{"epsilon_r", relative_permittivity}}, voltage, 0.0);
    mfem::Mesh mesh(mesh_file.c_str(), 1, 1);
    ElectrostaticSolver solver(mesh, DecodeConfig(config));
    solver.Setup();
    solver.Run();

    FieldExportSet fields = solver.CollectExportFields();
    const double permittivity = Constants::EPSILON_0 * relative_permittivity;
    const double field_integral = IntegrateVectorMagnitudeSquared(fields, "E");
    const double electric_energy = 0.5 * permittivity * field_integral;
    const double computed_capacitance = 2.0 * electric_energy / (voltage * voltage);
    const double analytic_capacitance = permittivity * height / length;

    REQUIRE(electric_energy == Catch::Approx(
        0.5 * analytic_capacitance * voltage * voltage).epsilon(1e-8));
    REQUIRE(computed_capacitance ==
        Catch::Approx(analytic_capacitance).epsilon(1e-8));

    const mfem::IntegrationPoint center = TriangleCenter();
    const double boundary_segment_length = height / ny;
    double outward_flux = 0.0;
    double absolute_flux = 0.0;
    for (int j = 0; j < ny; ++j) {
        const int left_element = 2 * (j * nx) + 1;
        const int right_element = 2 * (j * nx + nx - 1);
        const mfem::Vector left_field =
            SampleDerivedVector(fields, "E", left_element, center);
        const mfem::Vector right_field =
            SampleDerivedVector(fields, "E", right_element, center);
        const double left_flux = -permittivity * left_field(0) * boundary_segment_length;
        const double right_flux = permittivity * right_field(0) * boundary_segment_length;
        outward_flux += left_flux + right_flux;
        absolute_flux += std::abs(left_flux) + std::abs(right_flux);
    }
    REQUIRE(std::abs(outward_flux) < 1e-7 * absolute_flux);

    fs::remove(mesh_file);
}

TEST_CASE("Electrostatic solver satisfies dielectric interface conditions",
          "[solvers][analytic][electrostatic][materials][interface]") {
    const std::string mesh_file = "test_dielectric_interface.mesh";
    constexpr double length = 0.2;
    constexpr double height = 0.05;
    constexpr double voltage = 10.0;
    constexpr double epsilon_r_left = 2.0;
    constexpr double epsilon_r_right = 5.0;
    constexpr int nx = 8;
    constexpr int ny = 2;
    constexpr int interface_column = nx / 2;
    CreateLayeredStripMesh(
        mesh_file, length, height, nx, ny, interface_column);

    json config = {
        {"simulation", {
            {"physics_type", "electrostatics"},
            {"mesh", mesh_file},
            {"order", 1},
            {"geometry_type", "planar"},
            {"analysis_type", "field"},
            {"solver_tolerance", 1e-12},
            {"solver_max_iter", 4000},
            {"solver_print_level", 0}
        }},
        {"entity_groups", json::array({
            {{"name", "LeftDielectric"},  {"dim", 2}, {"attribute_ids", {1}}},
            {{"name", "RightDielectric"}, {"dim", 2}, {"attribute_ids", {2}}},
            {{"name", "Left"},  {"dim", 1}, {"attribute_ids", {1}}},
            {{"name", "Right"}, {"dim", 1}, {"attribute_ids", {2}}}
        })},
        {"regions", json::array({
            {{"name", "LeftDielectric"}, {"entity_group", "LeftDielectric"}, {"material", "LowPermittivity"}},
            {{"name", "RightDielectric"}, {"entity_group", "RightDielectric"}, {"material", "HighPermittivity"}}
        })},
        {"materials", json::array({
            {{"name", "LowPermittivity"}, {"properties", {{"epsilon_r", epsilon_r_left}}}},
            {{"name", "HighPermittivity"}, {"properties", {{"epsilon_r", epsilon_r_right}}}}
        })},
        {"boundaries", json::array({
            {{"name", "Left"},  {"type", "Dirichlet"}, {"entity_group", "Left"},  {"value", voltage}},
            {{"name", "Right"}, {"type", "Dirichlet"}, {"entity_group", "Right"}, {"value", 0.0}}
        })},
        {"scenarios", json::array({
            {{"name", "interface"}, {"excitations", json::array()}}
        })}
    };

    mfem::Mesh mesh(mesh_file.c_str(), 1, 1);
    ElectrostaticSolver solver(mesh, DecodeConfig(config));
    solver.Setup();
    solver.Run();

    FieldExportSet fields = solver.CollectExportFields();
    const mfem::IntegrationPoint center = TriangleCenter();
    const int left_element = 2 * (interface_column - 1);
    const int right_element = 2 * interface_column;
    const mfem::Vector left_field =
        SampleDerivedVector(fields, "E", left_element, center);
    const mfem::Vector right_field =
        SampleDerivedVector(fields, "E", right_element, center);
    const double epsilon_left = Constants::EPSILON_0 * epsilon_r_left;
    const double epsilon_right = Constants::EPSILON_0 * epsilon_r_right;
    const double left_width = length * interface_column / nx;
    const double right_width = length - left_width;
    const double analytic_displacement =
        voltage / (left_width / epsilon_left + right_width / epsilon_right);

    REQUIRE(epsilon_left * left_field(0) ==
        Catch::Approx(analytic_displacement).epsilon(1e-6));
    REQUIRE(epsilon_right * right_field(0) ==
        Catch::Approx(analytic_displacement).epsilon(1e-6));
    REQUIRE(std::abs(left_field(1)) < 1e-6 * std::abs(left_field(0)));
    REQUIRE(std::abs(right_field(1)) < 1e-6 * std::abs(right_field(0)));

    fs::remove(mesh_file);
}

TEST_CASE("Electrostatic capacitance matrix is analytic and reciprocal",
          "[solvers][analytic][electrostatic][coupling][reciprocity]") {
    const std::string mesh_file = "test_capacitance_reciprocity.mesh";
    const std::string matrix_file = "capacitance_matrix.csv";
    constexpr double length = 0.2;
    constexpr double height = 0.05;
    constexpr double relative_permittivity = 2.5;
    CreatePlanarStripMesh(mesh_file, length, height, 8, 2);

    json config = MakePlanarStripConfig(
        "electrostatics", mesh_file, 1,
        {{"epsilon_r", relative_permittivity}}, 0.0, 0.0);
    config["simulation"]["analysis_type"] = "coupling_matrix";
    config["boundaries"] = json::array();
    config["terminals"] = json::array({
        {{"name", "Left"}, {"excitation_type", "voltage"}, {"entity_group", "Left"}},
        {{"name", "Right"}, {"excitation_type", "voltage"}, {"entity_group", "Right"}}
    });

    mfem::Mesh mesh(mesh_file.c_str(), 1, 1);
    ElectrostaticSolver solver(mesh, DecodeConfig(config));
    solver.Setup();
    solver.Run();
    solver.SaveAnalysis();

    const CsvMatrix matrix = ReadCsvMatrix(matrix_file);
    REQUIRE(matrix.labels == std::vector<std::string>{"Left", "Right"});
    const double analytic_capacitance =
        Constants::EPSILON_0 * relative_permittivity * height / length;
    REQUIRE(matrix.values[0][1] == Catch::Approx(matrix.values[1][0]).epsilon(1e-7));
    REQUIRE(matrix.values[0][0] == Catch::Approx(analytic_capacitance).epsilon(1e-8));
    REQUIRE(matrix.values[1][1] == Catch::Approx(analytic_capacitance).epsilon(1e-8));
    REQUIRE(matrix.values[0][1] == Catch::Approx(-analytic_capacitance).epsilon(1e-8));
    REQUIRE(matrix.values[1][0] == Catch::Approx(-analytic_capacitance).epsilon(1e-8));
    REQUIRE(matrix.values[0][0] + matrix.values[0][1] ==
        Catch::Approx(0.0).margin(1e-7 * analytic_capacitance));
    REQUIRE(matrix.values[1][0] + matrix.values[1][1] ==
        Catch::Approx(0.0).margin(1e-7 * analytic_capacitance));

    fs::remove(matrix_file);
    fs::remove(mesh_file);
}

TEST_CASE("Electrostatic coupling ignores fixed Neumann background",
          "[solvers][electrostatic][coupling][m2]") {
    const std::string mesh_file = "test_es_coupling_background.mesh";
    const std::string matrix_file = "capacitance_matrix.csv";
    CreatePlanarStripMesh(mesh_file, 0.2, 0.05, 4, 2);

    json config = MakePlanarStripConfig(
        "electrostatics", mesh_file, 1, {{"epsilon_r", 2.5}}, 0.0, 0.0);
    config["simulation"]["analysis_type"] = "coupling_matrix";
    config["boundaries"] = json::array();
    config["terminals"] = json::array({
        {{"name", "Left"}, {"excitation_type", "voltage"}, {"entity_group", "Left"}},
        {{"name", "Right"}, {"excitation_type", "voltage"}, {"entity_group", "Right"}}
    });

    auto solve = [&]() {
        mfem::Mesh mesh(mesh_file.c_str(), 1, 1);
        ElectrostaticSolver solver(mesh, DecodeConfig(config));
        solver.Setup();
        solver.Run();
        solver.SaveAnalysis();
        return ReadCsvMatrix(matrix_file);
    };

    const CsvMatrix baseline = solve();
    config["entity_groups"].push_back(
        {{"name", "Horizontal"}, {"dim", 1}, {"attribute_ids", {3}}});
    config["boundaries"].push_back(
        {{"name", "BackgroundFlux"}, {"type", "Neumann"},
         {"entity_group", "Horizontal"}, {"value", 3.0}});
    const CsvMatrix with_background = solve();

    RequireMatricesEqual(with_background, baseline);
    fs::remove(matrix_file);
    fs::remove(mesh_file);
}

// Guards the axisymmetric 2*pi normalization end to end: the stiffness
// integrator omits the global 2*pi and GatherChargeColumn restores it exactly
// once. Applying it twice (or zero times) shifts every entry by a factor of
// 2*pi, which this analytic comparison catches.
TEST_CASE("Axisymmetric capacitance matches the analytic coaxial value",
          "[solvers][analytic][electrostatic][coupling][axisymmetric]") {
    const std::string mesh_file = "test_coax_capacitance.mesh";
    const std::string matrix_file = "capacitance_matrix.csv";
    constexpr double r_inner = 0.01;
    constexpr double r_outer = 0.03;
    constexpr double height = 0.05;
    CreateCoaxMesh(mesh_file, r_inner, r_outer, height, 64, 1);

    json config = MakeCoaxAmrConfig(mesh_file, 1);
    config["simulation"]["amr"]["enabled"] = false;
    config["simulation"]["analysis_type"] = "coupling_matrix";

    mfem::Mesh mesh(mesh_file.c_str(), 1, 1);
    ElectrostaticSolver solver(mesh, DecodeConfig(config));
    solver.Setup();
    solver.Run();
    solver.SaveAnalysis();

    const CsvMatrix matrix = ReadCsvMatrix(matrix_file);
    REQUIRE(matrix.labels == std::vector<std::string>{"Inner", "Outer"});

    // Coaxial annulus of axial extent `height` closed by symmetry planes:
    //   C = 2*pi*eps*height / ln(r_outer / r_inner)
    const double analytic_capacitance = Constants::TWO_PI * Constants::EPSILON_0 *
        height / std::log(r_outer / r_inner);

    REQUIRE(matrix.values[0][0] == Catch::Approx(analytic_capacitance).epsilon(0.01));
    REQUIRE(matrix.values[1][1] == Catch::Approx(analytic_capacitance).epsilon(0.01));
    REQUIRE(matrix.values[0][1] == Catch::Approx(-analytic_capacitance).epsilon(0.01));
    // Tolerance is set by the CSV's 6 significant digits, not by the solve.
    REQUIRE(matrix.values[0][1] == Catch::Approx(matrix.values[1][0]).epsilon(1e-5));

    fs::remove(matrix_file);
    fs::remove(mesh_file);
}

// The inductance matrix must be built from each MEASURED terminal's winding
// functional, not from the driving scenario's source. Reusing the drive makes
// every row of a column identical, which the asymmetry checks below reject.
TEST_CASE("Magnetostatic inductance matrix is reciprocal and distinguishes rows",
          "[solvers][magnetostatic][coupling][reciprocity]") {
    const std::string mesh_file = "test_magnetostatic_two_coil.mesh";
    const std::string matrix_file = "inductance_matrix.csv";
    CreateLayeredStripMesh(mesh_file, 0.2, 0.05, 4, 2, 2);

    json config = MakePlanarStripConfig(
        "magnetostatics", mesh_file, 1, {{"mu_r", 1.0}}, 0.0, 0.0);
    config["simulation"]["analysis_type"] = "coupling_matrix";
    config["entity_groups"].push_back(
        {{"name", "CoilA"}, {"dim", 2}, {"attribute_ids", {1}}});
    config["entity_groups"].push_back(
        {{"name", "CoilB"}, {"dim", 2}, {"attribute_ids", {2}}});
    config["regions"] = json::array({
        {{"name", "CoilA"}, {"entity_group", "CoilA"}, {"material", "Material"}},
        {{"name", "CoilB"}, {"entity_group", "CoilB"}, {"material", "Material"}}
    });
    config["terminals"] = json::array({
        {{"name", "CoilA"}, {"excitation_type", "current"}, {"entity_group", "CoilA"}},
        {{"name", "CoilB"}, {"excitation_type", "current"}, {"entity_group", "CoilB"}}
    });

    mfem::Mesh mesh(mesh_file.c_str(), 1, 1);
    MagnetostaticSolver solver(mesh, DecodeConfig(config));
    solver.Setup();
    solver.Run();
    solver.SaveAnalysis();

    const CsvMatrix matrix = ReadCsvMatrix(matrix_file);
    REQUIRE(matrix.labels == std::vector<std::string>{"CoilA", "CoilB"});
    REQUIRE(matrix.values[0][0] > 0.0);
    REQUIRE(matrix.values[1][1] > 0.0);
    REQUIRE(matrix.values[0][1] == Catch::Approx(matrix.values[1][0]).epsilon(1e-7));
    REQUIRE(matrix.values[1][0] < matrix.values[0][0]);
    REQUIRE(matrix.values[0][1] < matrix.values[1][1]);

    fs::remove(matrix_file);
    fs::remove(mesh_file);
}

TEST_CASE("Magnetostatic coupling ignores fixed Neumann background",
          "[solvers][magnetostatic][coupling][m2]") {
    const std::string mesh_file = "test_ms_coupling_background.mesh";
    const std::string matrix_file = "inductance_matrix.csv";
    CreateLayeredStripMesh(mesh_file, 0.2, 0.05, 4, 2, 2);

    json config = MakePlanarStripConfig(
        "magnetostatics", mesh_file, 1, {{"mu_r", 1.0}}, 0.0, 0.0);
    config["simulation"]["analysis_type"] = "coupling_matrix";
    config["entity_groups"].push_back(
        {{"name", "CoilA"}, {"dim", 2}, {"attribute_ids", {1}}});
    config["entity_groups"].push_back(
        {{"name", "CoilB"}, {"dim", 2}, {"attribute_ids", {2}}});
    config["regions"] = json::array({
        {{"name", "CoilA"}, {"entity_group", "CoilA"}, {"material", "Material"}},
        {{"name", "CoilB"}, {"entity_group", "CoilB"}, {"material", "Material"}}
    });
    config["terminals"] = json::array({
        {{"name", "CoilA"}, {"excitation_type", "current"}, {"entity_group", "CoilA"}},
        {{"name", "CoilB"}, {"excitation_type", "current"}, {"entity_group", "CoilB"}}
    });

    auto solve = [&]() {
        mfem::Mesh mesh(mesh_file.c_str(), 1, 1);
        MagnetostaticSolver solver(mesh, DecodeConfig(config));
        solver.Setup();
        solver.Run();
        return ReadCsvMatrix(matrix_file);
    };

    const CsvMatrix baseline = solve();
    config["entity_groups"].push_back(
        {{"name", "Horizontal"}, {"dim", 1}, {"attribute_ids", {3}}});
    config["boundaries"].push_back(
        {{"name", "BackgroundFlux"}, {"type", "Neumann"},
         {"entity_group", "Horizontal"}, {"value", 2.0}});
    const CsvMatrix with_background = solve();

    RequireMatricesEqual(with_background, baseline);
    fs::remove(matrix_file);
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
        {{"name", "Current"}, {"excitation_type", "current"},
         {"conductor_type", "stranded"}, {"entity_group", "Domain"}}
    });
    config["scenarios"][0]["excitations"] = json::array({
        {{"terminal", "Current"}, {"value", current_density * length * height}}
    });

    mfem::Mesh mesh(mesh_file.c_str(), 1, 1);
    MagnetostaticSolver solver(mesh, DecodeConfig(config));
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

TEST_CASE("Magnetoquasistatic frequency sweep updates conducting-slab skin effect",
          "[solvers][analytic][mqs]") {
    const std::string mesh_file = "test_analytic_mqs.mesh";
    constexpr double length = 0.04;
    constexpr double height = 0.005;
    constexpr double conductivity = 3.5e7;
    constexpr double low_frequency = 60.0;
    constexpr double frequency = 600.0;
    constexpr int nx = 48;
    constexpr int ny = 1;
    CreatePlanarStripMesh(mesh_file, length, height, nx, ny);

    json config = MakePlanarStripConfig(
        "magnetoquasistatics", mesh_file, 2,
        {{"mu_r", 1.0}, {"sigma", conductivity}}, 1.0, 0.0);
    config["scenarios"] = json::array({
        {{"name", "skin"},
         {"frequency", {{"scale", "linear"}, {"start", low_frequency},
                        {"stop", frequency}, {"points", 2}}},
         {"excitations", json::array()}}
    });
    config["simulation"]["amr"] = {{"enabled", true}, {"max_iterations", 1}};

    mfem::Mesh mesh(mesh_file.c_str(), 1, 1);
    MagnetoquasistaticSolver solver(mesh, DecodeConfig(config));
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
    const std::complex<double> low_wave_number =
        std::sqrt(std::complex<double>(
            0.0, Constants::TWO_PI * low_frequency * Constants::MU_0 * conductivity));
    const std::complex<double> low_expected =
        std::sinh(low_wave_number * (length - physical(0))) /
        std::sinh(low_wave_number * length);

    const double actual_real = SamplePrimaryScalar(fields, "A_Real", element, point);
    const double actual_imag = SamplePrimaryScalar(fields, "A_Imag", element, point);

    REQUIRE(actual_real == Catch::Approx(expected.real()).epsilon(5e-3));
    REQUIRE(actual_imag == Catch::Approx(expected.imag()).epsilon(5e-3));
    REQUIRE(std::abs(std::complex<double>(actual_real, actual_imag) - expected) <
            std::abs(std::complex<double>(actual_real, actual_imag) - low_expected));
    REQUIRE(std::hypot(actual_real, actual_imag) < 1.0);
    REQUIRE(actual_imag < 0.0);

    fs::remove(mesh_file);
}

TEST_CASE("Magnetoquasistatic skin-effect solution converges under mesh refinement",
          "[solvers][analytic][mqs][convergence]") {
    constexpr double length = 0.04;
    constexpr double height = 0.005;
    constexpr double conductivity = 3.5e7;
    constexpr double frequency = 60.0;
    const double omega = Constants::TWO_PI * frequency;
    const std::complex<double> wave_number =
        std::sqrt(std::complex<double>(0.0, omega * Constants::MU_0 * conductivity));
    const auto exact = [=](const mfem::Vector& point) {
        return std::sinh(wave_number * (length - point(0))) /
               std::sinh(wave_number * length);
    };

    std::vector<double> errors;
    for (const int nx : {8, 16, 32}) {
        const std::string mesh_file =
            "test_mqs_convergence_" + std::to_string(nx) + ".mesh";
        CreatePlanarStripMesh(mesh_file, length, height, nx, 1);

        json config = MakePlanarStripConfig(
            "magnetoquasistatics", mesh_file, 1,
            {{"mu_r", 1.0}, {"sigma", conductivity}}, 1.0, 0.0);
        config["scenarios"][0]["frequency"] = frequency;

        mfem::Mesh mesh(mesh_file.c_str(), 1, 1);
        MagnetoquasistaticSolver solver(mesh, DecodeConfig(config));
        solver.Setup();
        solver.Run();
        errors.push_back(RelativeComplexL2Error(solver.CollectExportFields(), exact));

        fs::remove(mesh_file);
    }

    REQUIRE(errors[1] < errors[0]);
    REQUIRE(errors[2] < errors[1]);
    REQUIRE(errors[0] / errors[1] > 2.5);
    REQUIRE(errors[1] / errors[2] > 2.5);
}

TEST_CASE("MQS coupling supports mixed massive and stranded conductors",
          "[solvers][mqs][coupling][reciprocity]") {
    const std::string mesh_file = "test_mqs_mixed_conductors.mesh";
    const std::string low_inductance_file = "inductance_matrix_low_100Hz.csv";
    const std::string low_resistance_file = "resistance_matrix_low_100Hz.csv";
    const std::string high_inductance_file = "inductance_matrix_high_1000Hz.csv";
    const std::string high_resistance_file = "resistance_matrix_high_1000Hz.csv";
    CreateLayeredStripMesh(mesh_file, 0.2, 0.05, 2, 1, 1);

    json config = MakePlanarStripConfig(
        "magnetoquasistatics", mesh_file, 1, {{"mu_r", 1.0}}, 0.0, 0.0);
    config["simulation"]["analysis_type"] = "coupling_matrix";
    config["scenarios"] = json::array({
        {{"name", "low"}, {"frequency", 100.0}, {"excitations", json::array()}},
        {{"name", "high"}, {"frequency", 1000.0}, {"excitations", json::array()}}
    });
    config["entity_groups"].push_back(
        {{"name", "StrandedDomain"}, {"dim", 2}, {"attribute_ids", {1}}});
    config["entity_groups"].push_back(
        {{"name", "MassiveDomain"}, {"dim", 2}, {"attribute_ids", {2}}});
    config["materials"] = json::array({
        {{"name", "StrandedMaterial"},
         {"properties", {{"mu_r", 1.0}, {"sigma", 0.0}}}},
        {{"name", "MassiveMaterial"},
         {"properties", {{"mu_r", 1.0}, {"sigma", 1.0e6}}}}
    });
    config["regions"] = json::array({
        {{"name", "StrandedRegion"}, {"entity_group", "StrandedDomain"},
         {"material", "StrandedMaterial"}},
        {{"name", "MassiveRegion"}, {"entity_group", "MassiveDomain"},
         {"material", "MassiveMaterial"}}
    });
    config["terminals"] = json::array({
        {{"name", "Massive"}, {"excitation_type", "current"},
         {"conductor_type", "massive"}, {"entity_group", "MassiveDomain"}},
        {{"name", "Stranded"}, {"excitation_type", "current"},
         {"conductor_type", "stranded"}, {"entity_group", "StrandedDomain"}}
    });

    mfem::Mesh mesh(mesh_file.c_str(), 1, 1);
    MagnetoquasistaticSolver solver(mesh, DecodeConfig(config));
    solver.Setup();
    solver.Run();
    solver.SaveAnalysis();

    REQUIRE(fs::exists(low_inductance_file));
    REQUIRE(fs::exists(low_resistance_file));
    REQUIRE(fs::exists(high_inductance_file));
    REQUIRE(fs::exists(high_resistance_file));
    const CsvMatrix low_inductance = ReadCsvMatrix(low_inductance_file);
    const CsvMatrix low_resistance = ReadCsvMatrix(low_resistance_file);
    const CsvMatrix inductance = ReadCsvMatrix(high_inductance_file);
    const CsvMatrix resistance = ReadCsvMatrix(high_resistance_file);
    const std::vector<std::string> expected_labels{"Massive", "Stranded"};
    REQUIRE(inductance.labels == expected_labels);
    REQUIRE(resistance.labels == expected_labels);
    REQUIRE(inductance.values[0][1] ==
        Catch::Approx(inductance.values[1][0]).epsilon(1e-7));
    REQUIRE(resistance.values[0][1] ==
        Catch::Approx(resistance.values[1][0]).epsilon(1e-7));
    REQUIRE(inductance.values[0][0] > 0.0);
    REQUIRE(inductance.values[1][1] > 0.0);
    REQUIRE(resistance.values[0][0] > 0.0);
    REQUIRE(std::abs(inductance.values[0][0] - low_inductance.values[0][0]) > 1e-12);
    REQUIRE(std::abs(resistance.values[0][0] - low_resistance.values[0][0]) > 1e-12);

    fs::remove(low_resistance_file);
    fs::remove(low_inductance_file);
    fs::remove(high_resistance_file);
    fs::remove(high_inductance_file);
    fs::remove(mesh_file);
}

TEST_CASE("MQS coupling ignores fixed Neumann background",
          "[solvers][mqs][coupling][m2]") {
    const std::string mesh_file = "test_mqs_coupling_background.mesh";
    const std::string inductance_file = "inductance_matrix_point_1000Hz.csv";
    const std::string resistance_file = "resistance_matrix_point_1000Hz.csv";
    CreateLayeredStripMesh(mesh_file, 0.2, 0.05, 2, 1, 1);

    json config = MakePlanarStripConfig(
        "magnetoquasistatics", mesh_file, 1,
        {{"mu_r", 1.0}, {"sigma", 1.0e6}}, 0.0, 0.0);
    config["simulation"]["analysis_type"] = "coupling_matrix";
    config["scenarios"] = json::array({
        {{"name", "point"}, {"frequency", 1000.0}, {"excitations", json::array()}}
    });
    config["entity_groups"].push_back(
        {{"name", "DriveDomain"}, {"dim", 2}, {"attribute_ids", {1}}});
    config["regions"] = json::array({
        {{"name", "Domain"}, {"entity_group", "Domain"}, {"material", "Material"}}
    });
    config["terminals"] = json::array({
        {{"name", "Drive"}, {"excitation_type", "current"},
         {"conductor_type", "massive"}, {"entity_group", "DriveDomain"}}
    });

    auto solve = [&]() {
        mfem::Mesh mesh(mesh_file.c_str(), 1, 1);
        MagnetoquasistaticSolver solver(mesh, DecodeConfig(config));
        solver.Setup();
        solver.Run();
        solver.SaveAnalysis();
        return std::pair{
            ReadCsvMatrix(inductance_file), ReadCsvMatrix(resistance_file)};
    };

    const auto baseline = solve();
    config["entity_groups"].push_back(
        {{"name", "Horizontal"}, {"dim", 1}, {"attribute_ids", {3}}});
    config["boundaries"].push_back(
        {{"name", "BackgroundFlux"}, {"type", "Neumann"},
         {"entity_group", "Horizontal"}, {"value", 2.0}});
    const auto with_background = solve();

    RequireMatricesEqual(with_background.first, baseline.first);
    RequireMatricesEqual(with_background.second, baseline.second);
    fs::remove(resistance_file);
    fs::remove(inductance_file);
    fs::remove(mesh_file);
}

TEST_CASE("MQS heterogeneous massive-port conductance is piecewise and order independent",
          "[solvers][mqs][coupling][conductance]") {
    const std::string mesh_file = "test_mqs_heterogeneous_port.mesh";
    const std::string inductance_file =
        "inductance_matrix_heterogeneous_60Hz.csv";
    const std::string resistance_file =
        "resistance_matrix_heterogeneous_60Hz.csv";
    constexpr double length = 0.2;
    constexpr double height = 0.05;
    constexpr double sigma_left = 2.0;
    constexpr double sigma_right = 5.0;
    CreateLayeredStripMesh(mesh_file, length, height, 2, 1, 1);

    auto resistance_for = [&](const json& port_attributes) {
        json config = MakePlanarStripConfig(
            "magnetoquasistatics", mesh_file, 1,
            {{"mu_r", 1.0}, {"sigma", sigma_left}}, 0.0, 0.0);
        config["simulation"]["analysis_type"] = "coupling_matrix";
        config["scenarios"] = json::array({
            {{"name", "heterogeneous"}, {"frequency", 60.0},
             {"excitations", json::array()}}
        });
        config["entity_groups"].push_back(
            {{"name", "RightDomain"}, {"dim", 2}, {"attribute_ids", {2}}});
        config["entity_groups"].push_back(
            {{"name", "PortDomain"}, {"dim", 2},
             {"attribute_ids", port_attributes}});
        config["entity_groups"].push_back(
            {{"name", "Horizontal"}, {"dim", 1}, {"attribute_ids", {3}}});
        config["materials"].push_back(
            {{"name", "RightMaterial"},
             {"properties", {{"mu_r", 1.0}, {"sigma", sigma_right}}}});
        config["regions"].push_back(
            {{"name", "RightRegion"}, {"entity_group", "RightDomain"},
             {"material", "RightMaterial"}});
        config["boundaries"].push_back(
            {{"name", "Horizontal"}, {"type", "Dirichlet"},
             {"entity_group", "Horizontal"}, {"value", 0.0}});
        config["terminals"] = json::array({
            {{"name", "Port"}, {"excitation_type", "current"},
             {"conductor_type", "massive"}, {"entity_group", "PortDomain"}}
        });

        mfem::Mesh mesh(mesh_file.c_str(), 1, 1);
        MagnetoquasistaticSolver solver(mesh, DecodeConfig(config));
        solver.Setup();
        solver.Run();
        solver.SaveAnalysis();
        const CsvMatrix resistance = ReadCsvMatrix(resistance_file);
        REQUIRE(resistance.labels == std::vector<std::string>{"Port"});
        REQUIRE(resistance.values.size() == 1);
        REQUIRE(resistance.values[0].size() == 1);
        return resistance.values[0][0];
    };

    const double left_area = 0.5 * length * height;
    const double right_area = 0.5 * length * height;
    const double expected_resistance =
        1.0 / (sigma_left * left_area + sigma_right * right_area);
    const double forward = resistance_for(json::array({1, 2}));
    const double reversed = resistance_for(json::array({2, 1}));

    REQUIRE(forward == Catch::Approx(expected_resistance).epsilon(1e-6));
    REQUIRE(reversed == Catch::Approx(expected_resistance).epsilon(1e-6));
    REQUIRE(reversed == Catch::Approx(forward).epsilon(1e-12));

    fs::remove(resistance_file);
    fs::remove(inductance_file);
    fs::remove(mesh_file);
}

TEST_CASE("MQS coupling keeps open-current regions passive and off-matrix",
          "[solvers][mqs][coupling][regions]") {
    const std::string mesh_file = "test_mqs_open_current_region.mesh";
    const std::string inductance_file = "inductance_matrix_point_f1_1000Hz.csv";
    const std::string resistance_file = "resistance_matrix_point_f1_1000Hz.csv";
    const std::string legacy_inductance_file =
        "inductance_matrix_point_f1_1000Hz_1000Hz.csv";
    const std::string legacy_resistance_file =
        "resistance_matrix_point_f1_1000Hz_1000Hz.csv";
    fs::remove(legacy_inductance_file);
    fs::remove(legacy_resistance_file);
    CreateLayeredStripMesh(mesh_file, 0.2, 0.05, 4, 1, 2);

    json config = MakePlanarStripConfig(
        "magnetoquasistatics", mesh_file, 1,
        {{"mu_r", 1.0}, {"sigma", 1.0e6}}, 0.0, 0.0);
    config["simulation"]["analysis_type"] = "coupling_matrix";
    config["scenarios"] = json::array({
        {{"name", "point"},
         {"frequency", {{"scale", "linear"}, {"start", 1000.0},
                          {"stop", 1000.0}, {"points", 1}}},
         {"excitations", json::array()}}
    });
    config["entity_groups"].push_back(
        {{"name", "DriveDomain"}, {"dim", 2}, {"attribute_ids", {1}}});
    config["entity_groups"].push_back(
        {{"name", "PassiveDomain"}, {"dim", 2}, {"attribute_ids", {2}}});
    config["regions"] = json::array({
        {{"name", "DriveRegion"}, {"entity_group", "DriveDomain"}, {"material", "Material"}},
        {{"name", "PassiveRegion"}, {"entity_group", "PassiveDomain"}, {"material", "Material"}}
    });
    config["terminals"] = json::array({
        {{"name", "Drive"}, {"excitation_type", "current"},
         {"conductor_type", "massive"}, {"entity_group", "DriveDomain"}}
    });

    mfem::Mesh baseline_mesh(mesh_file.c_str(), 1, 1);
    MagnetoquasistaticSolver baseline_solver(baseline_mesh, DecodeConfig(config));
    baseline_solver.Setup();
    baseline_solver.Run();
    baseline_solver.SaveAnalysis();
    const CsvMatrix baseline_inductance = ReadCsvMatrix(inductance_file);
    const CsvMatrix baseline_resistance = ReadCsvMatrix(resistance_file);

    config["regions"][1]["current_constraint"] = "open";
    mfem::Mesh constrained_mesh(mesh_file.c_str(), 1, 1);
    MagnetoquasistaticSolver constrained_solver(constrained_mesh, DecodeConfig(config));
    constrained_solver.Setup();
    constrained_solver.Run();
    constrained_solver.SaveAnalysis();
    const CsvMatrix constrained_inductance = ReadCsvMatrix(inductance_file);
    const CsvMatrix constrained_resistance = ReadCsvMatrix(resistance_file);

    const std::vector<std::string> expected_labels{"Drive"};
    REQUIRE(constrained_inductance.labels == expected_labels);
    REQUIRE(constrained_resistance.labels == expected_labels);
    REQUIRE(constrained_inductance.values.size() == 1);
    REQUIRE(constrained_resistance.values.size() == 1);
    REQUIRE_FALSE(fs::exists(legacy_inductance_file));
    REQUIRE_FALSE(fs::exists(legacy_resistance_file));
    const double inductance_change = std::abs(
        constrained_inductance.values[0][0] - baseline_inductance.values[0][0]);
    const double resistance_change = std::abs(
        constrained_resistance.values[0][0] - baseline_resistance.values[0][0]);
    REQUIRE(inductance_change + resistance_change > 1e-12);

    fs::remove(resistance_file);
    fs::remove(inductance_file);
    fs::remove(mesh_file);
}

// Axisymmetric companion to the planar open-current test: an aluminum shield
// modeled as a passive open-current region, with two driven massive turns
// outside it. Covers the axisymmetric massive-port path (G_dc = integral
// sigma/(2*pi*r)) for a region that is not a terminal, which the planar test
// cannot exercise.
//
// Mesh resolution note: the shield must be resolved against the skin depth
// delta = sqrt(2/(omega*mu*sigma)). At 1 kHz in aluminum delta is about 2.7 mm,
// and the radial cell size here is about 2.7 mm, so delta/cell is roughly 1.
// Coarser meshes drive the extracted port resistance negative (unphysical
// negative loss), so do not reduce cells_per_band without also lowering the
// frequency.
TEST_CASE("MQS axisymmetric shield stays passive and loads the turns",
          "[solvers][mqs][coupling][regions][axisymmetric]") {
    const std::string mesh_file = "test_mqs_axisym_shield.mesh";
    const std::string inductance_file = "inductance_matrix_point_f1_1000Hz.csv";
    const std::string resistance_file = "resistance_matrix_point_f1_1000Hz.csv";
    CreateShieldedTurnsMesh(mesh_file, /*r_min=*/0.05, /*r_max=*/0.20,
                            /*height=*/0.04, /*nz=*/8, /*cells_per_band=*/8);

    json config = MakeShieldedTurnsConfig(mesh_file, 1000.0);

    // Reference 1: no shield at all. The shield band is air, so nothing is
    // induced there and the turns see their unshielded coupling.
    json air_config = config;
    air_config["regions"][1]["material"] = "Air";
    mfem::Mesh air_mesh(mesh_file.c_str(), 1, 1);
    MagnetoquasistaticSolver air_solver(air_mesh, DecodeConfig(air_config));
    air_solver.Setup();
    air_solver.Run();
    air_solver.SaveAnalysis();
    const CsvMatrix air_inductance = ReadCsvMatrix(inductance_file);

    // Reference 2: a conducting shield with no current constraint. Without a
    // port unknown there is no net-current equation, so the region behaves as a
    // short-circuited closed loop and shields most strongly.
    mfem::Mesh shorted_mesh(mesh_file.c_str(), 1, 1);
    MagnetoquasistaticSolver shorted_solver(shorted_mesh, DecodeConfig(config));
    shorted_solver.Setup();
    shorted_solver.Run();
    shorted_solver.SaveAnalysis();
    const CsvMatrix shorted_inductance = ReadCsvMatrix(inductance_file);

    // Case under test: the shield becomes an open-current region, gaining a
    // port unknown whose net current is pinned to zero. This is the physical
    // model for a shield that is not a closed turn.
    config["regions"][1]["current_constraint"] = "open";
    mfem::Mesh shielded_mesh(mesh_file.c_str(), 1, 1);
    MagnetoquasistaticSolver shielded_solver(shielded_mesh, DecodeConfig(config));
    shielded_solver.Setup();
    shielded_solver.Run();
    shielded_solver.SaveAnalysis();
    const CsvMatrix shielded_inductance = ReadCsvMatrix(inductance_file);
    const CsvMatrix shielded_resistance = ReadCsvMatrix(resistance_file);

    // The shield is passive: it never becomes a matrix row/column, so the
    // matrices stay 2x2 over the two driven turns.
    const std::vector<std::string> expected_labels{"TurnA", "TurnB"};
    REQUIRE(shielded_inductance.labels == expected_labels);
    REQUIRE(shielded_resistance.labels == expected_labels);
    REQUIRE(shielded_inductance.values.size() == 2);
    REQUIRE(shielded_resistance.values.size() == 2);

    // Reciprocity: mutual terms must match in both matrices.
    REQUIRE(shielded_inductance.values[0][1] ==
        Catch::Approx(shielded_inductance.values[1][0]).epsilon(1e-8));
    REQUIRE(shielded_resistance.values[0][1] ==
        Catch::Approx(shielded_resistance.values[1][0]).epsilon(1e-8));

    // Self terms stay physical: positive inductance and positive loss.
    for (int turn = 0; turn < 2; ++turn) {
        REQUIRE(shielded_inductance.values[turn][turn] > 0.0);
        REQUIRE(shielded_resistance.values[turn][turn] > 0.0);
    }

    // The shield is not inert: constraining net current to zero still permits
    // eddy currents to circulate within the cross-section, so the turns must
    // see different coupling than with no shield present.
    double change_vs_air = 0.0;
    for (int row = 0; row < 2; ++row) {
        for (int column = 0; column < 2; ++column) {
            change_vs_air += std::abs(shielded_inductance.values[row][column] -
                                      air_inductance.values[row][column]);
        }
    }
    REQUIRE(change_vs_air > 1e-14);

    // Shielding strength is ordered by how much net current the shield may
    // carry. A short-circuited closed loop develops a large opposing net
    // current and shields strongly. Pinning the net current to zero removes
    // that bulk opposition entirely, leaving only redistribution within the
    // cross-section, so the turns behave nearly as if the shield were absent.
    for (int turn = 0; turn < 2; ++turn) {
        const double shorted = shorted_inductance.values[turn][turn];
        const double open = shielded_inductance.values[turn][turn];
        const double air = air_inductance.values[turn][turn];

        // The closed loop shields substantially, by a margin well outside
        // discretization noise. The exact factor depends on how far the turn
        // sits from the shield, so only the direction is asserted.
        REQUIRE(shorted < 0.8 * air);
        // The zero-net-current shield does not: it stays within a few percent
        // of the unshielded value, and clearly apart from the shorted case.
        REQUIRE(std::abs(open - air) < 0.05 * air);
        REQUIRE(shorted < 0.9 * open);
    }

    fs::remove(resistance_file);
    fs::remove(inductance_file);
    fs::remove(mesh_file);
}

// Cross-check on the shielding mechanism: a short-circuited closed shield
// excludes flux more strongly as frequency rises, because the induced opposing
// net current grows with the rate of change of flux. A zero-net-current (open)
// shield has no such bulk mechanism available, so it stays near the unshielded
// value at every frequency. Asserting the trend, rather than a single operating
// point, is what separates the two constraints physically.
TEST_CASE("MQS shielding strengthens with frequency only when current may close",
          "[solvers][mqs][coupling][regions][axisymmetric]") {
    const std::string mesh_file = "test_mqs_axisym_shield_sweep.mesh";
    // The radial cell must stay below the aluminum skin depth at the top of the
    // sweep (delta is about 0.85 mm at 10 kHz). At 32 cells per band the
    // extracted port resistance is still negative; 64 resolves it.
    CreateShieldedTurnsMesh(mesh_file, /*r_min=*/0.05, /*r_max=*/0.20,
                            /*height=*/0.04, /*nz=*/8, /*cells_per_band=*/64);

    // Self inductance of turn A with the shield present, divided by the same
    // quantity with the shield replaced by air. Lower means more flux excluded.
    auto shielding_ratio = [&](double frequency, bool closed) {
        std::ostringstream tag;
        tag << "point_f1_" << frequency << "Hz";
        const std::string inductance_file = "inductance_matrix_" + tag.str() + ".csv";
        const std::string resistance_file = "resistance_matrix_" + tag.str() + ".csv";

        json config = MakeShieldedTurnsConfig(mesh_file, frequency);
        if (!closed) config["regions"][1]["current_constraint"] = "open";

        json air_config = MakeShieldedTurnsConfig(mesh_file, frequency);
        air_config["regions"][1]["material"] = "Air";

        mfem::Mesh air_mesh(mesh_file.c_str(), 1, 1);
        MagnetoquasistaticSolver air_solver(air_mesh, DecodeConfig(air_config));
        air_solver.Setup();
        air_solver.Run();
        air_solver.SaveAnalysis();
        const CsvMatrix air = ReadCsvMatrix(inductance_file);

        mfem::Mesh shield_mesh(mesh_file.c_str(), 1, 1);
        MagnetoquasistaticSolver shield_solver(shield_mesh, DecodeConfig(config));
        shield_solver.Setup();
        shield_solver.Run();
        shield_solver.SaveAnalysis();
        const CsvMatrix shielded = ReadCsvMatrix(inductance_file);
        const CsvMatrix resistance = ReadCsvMatrix(resistance_file);

        // Guard the discretization: negative extracted loss means the shield is
        // under-resolved against the skin depth and the ratio is meaningless.
        for (int turn = 0; turn < 2; ++turn) {
            REQUIRE(resistance.values[turn][turn] > 0.0);
        }

        fs::remove(inductance_file);
        fs::remove(resistance_file);
        return shielded.values[0][0] / air.values[0][0];
    };

    // A closed (short-circuited) shield excludes more flux as frequency rises.
    const double closed_low = shielding_ratio(100.0, /*closed=*/true);
    const double closed_high = shielding_ratio(10000.0, /*closed=*/true);
    REQUIRE(closed_low < 1.0);
    REQUIRE(closed_high < closed_low);

    // An open shield cannot carry net current, so no such trend develops: it
    // stays near the unshielded value at both ends of the sweep.
    const double open_low = shielding_ratio(100.0, /*closed=*/false);
    const double open_high = shielding_ratio(10000.0, /*closed=*/false);
    REQUIRE(std::abs(open_low - 1.0) < 0.05);
    REQUIRE(std::abs(open_high - 1.0) < 0.05);
    REQUIRE(closed_high < 0.5 * open_high);

    fs::remove(mesh_file);
}


TEST_CASE("AMR refines an axisymmetric coax and stays conforming", "[solvers][amr]") {
    const std::string mesh_file = "test_amr_coax.mesh";
    CreateCoaxMesh(mesh_file, /*r_inner=*/1.0, /*r_outer=*/4.0,
                   /*height=*/1.0, /*nr=*/6, /*nz=*/4);

    json config = MakeCoaxAmrConfig(mesh_file, /*max_iterations=*/4);

    mfem::Mesh mesh(mesh_file.c_str(), 1, 1);

    ElectrostaticAmrProbe solver(mesh, DecodeConfig(config));
    solver.Setup();
    REQUIRE_NOTHROW(solver.Run());

    const auto& history = solver.AmrHistory();
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
    const double last_peak = history.back().peak_field_magnitude;
    REQUIRE(std::isfinite(last_peak));
    REQUIRE(last_peak > 0.5 * analytic_peak);
    REQUIRE(last_peak < 1.2 * analytic_peak);

    fs::remove(mesh_file);
}

TEST_CASE("Axisymmetric electrostatic AMR applies permittivity once",
          "[solvers][amr][materials][m1]") {
    const std::string mesh_file = "test_amr_axisymmetric_es_materials.mesh";
    CreateCoaxMesh(mesh_file, 1.0, 4.0, 1.0, 6, 2, 3);

    json config = MakeCoaxAmrConfig(mesh_file, 1);
    config["simulation"]["amr"]["enabled"] = false;
    config["entity_groups"].push_back(
        {{"name", "OuterDielectric"}, {"dim", 2}, {"attribute_ids", {2}}});
    config["materials"].push_back(
        {{"name", "HighPermittivity"}, {"properties", {{"epsilon_r", 100.0}}}});
    config["regions"].push_back(
        {{"name", "OuterDielectric"}, {"entity_group", "OuterDielectric"},
         {"material", "HighPermittivity"}});

    mfem::Mesh mesh(mesh_file.c_str(), 1, 1);
    ElectrostaticAmrProbe solver(mesh, DecodeConfig(config));
    solver.Setup();
    solver.Run();

    const mfem::Vector actual = solver.LocalErrors();
    FieldExportSet fields = solver.CollectExportFields();
    mfem::GridFunction& potential = *FindField(fields, "V").primary;
    mfem::Vector epsilon_values(2);
    epsilon_values(0) = Constants::EPSILON_0;
    epsilon_values(1) = 100.0 * Constants::EPSILON_0;
    mfem::PWConstCoefficient epsilon(epsilon_values);
    const mfem::Vector intended =
        AxisymmetricDiffusionErrors(mesh, potential, epsilon, false);
    const mfem::Vector coefficient_cubed =
        AxisymmetricDiffusionErrors(mesh, potential, epsilon, true);

    RequireErrorsEqual(actual, intended);
    REQUIRE(ErrorsDiffer(actual, coefficient_cubed));

    mfem::Array<int> actual_marked;
    mfem::Array<int> intended_marked;
    amr::MarkElementsDorfler(actual, 0.7, actual_marked);
    amr::MarkElementsDorfler(intended, 0.7, intended_marked);
    REQUIRE(actual_marked == intended_marked);

    fs::remove(mesh_file);
}

TEST_CASE("Axisymmetric magnetostatic AMR applies reluctivity once",
          "[solvers][amr][materials][m1]") {
    const std::string mesh_file = "test_amr_axisymmetric_ms_materials.mesh";
    CreateCoaxMesh(mesh_file, 1.0, 4.0, 1.0, 6, 2, 3);

    json config = MakeCoaxAmrConfig(mesh_file, 1);
    config["simulation"]["physics_type"] = "magnetostatics";
    config["simulation"]["amr"]["enabled"] = false;
    config["entity_groups"].push_back(
        {{"name", "OuterMagnetic"}, {"dim", 2}, {"attribute_ids", {2}}});
    config["materials"] = json::array({
        {{"name", "LowPermeability"}, {"properties", {{"mu_r", 1.0}}}},
        {{"name", "HighPermeability"}, {"properties", {{"mu_r", 100.0}}}}
    });
    config["regions"] = json::array({
        {{"name", "InnerMagnetic"}, {"entity_group", "Dielectric"},
         {"material", "LowPermeability"}},
        {{"name", "OuterMagnetic"}, {"entity_group", "OuterMagnetic"},
         {"material", "HighPermeability"}}
    });
    config.erase("terminals");
    config["boundaries"] = json::array({
        {{"name", "Inner"}, {"type", "Dirichlet"},
         {"entity_group", "Inner"}, {"value", 1.0}},
        {{"name", "Outer"}, {"type", "Dirichlet"},
         {"entity_group", "Outer"}, {"value", 0.0}}
    });
    config["scenarios"] = json::array({
        {{"name", "energized"}, {"excitations", json::array()}}
    });

    mfem::Mesh mesh(mesh_file.c_str(), 1, 1);
    MagnetostaticAmrProbe solver(mesh, DecodeConfig(config));
    solver.Setup();
    solver.Run();

    const mfem::Vector actual = solver.LocalErrors();
    FieldExportSet fields = solver.CollectExportFields();
    mfem::GridFunction& potential = *FindField(fields, "A").primary;
    mfem::Vector reluctivity_values(2);
    reluctivity_values(0) = 1.0 / Constants::MU_0;
    reluctivity_values(1) = 1.0 / (100.0 * Constants::MU_0);
    mfem::PWConstCoefficient reluctivity(reluctivity_values);
    const mfem::Vector intended =
        AxisymmetricCurlCurlErrors(mesh, potential, reluctivity, false);
    const mfem::Vector coefficient_cubed =
        AxisymmetricCurlCurlErrors(mesh, potential, reluctivity, true);

    RequireErrorsEqual(actual, intended);
    REQUIRE(ErrorsDiffer(actual, coefficient_cubed));

    mfem::Array<int> actual_marked;
    mfem::Array<int> intended_marked;
    amr::MarkElementsDorfler(actual, 0.7, actual_marked);
    amr::MarkElementsDorfler(intended, 0.7, intended_marked);
    REQUIRE(actual_marked == intended_marked);

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

    mfem::Mesh mesh(mesh_in_tmp.string().c_str(), 1, 1);

    ElectrostaticSolver solver(mesh, DecodeConfig(config));
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
