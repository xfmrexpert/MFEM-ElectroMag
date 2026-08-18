// Copyright (c) 2026 T. C. Raymond
// SPDX-License-Identifier: MIT

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include "../src/input_parser.hpp"
#include <fstream>
#include <filesystem>

namespace fs = std::filesystem;

namespace {

json CanonicalConfig() {
    return json{
        {"simulation", {
            {"physics_type", "magnetoquasistatics"},
            {"geometry_type", "axisymmetric"},
            {"analysis_type", "coupling_matrix"},
            {"mesh", "test.msh"},
            {"order", 2},
            {"solver_tolerance", 1e-9},
            {"solver_max_iter", 321},
            {"solver_print_level", 2},
            {"output_paraview", true},
            {"output_gmsh", true},
            {"export_refine", 3},
            {"amr", {
                {"enabled", true},
                {"max_iterations", 4},
                {"max_dofs", 12345},
                {"error_fraction", 0.6},
                {"error_tolerance", 1e-5},
                {"conforming", true}
            }}
        }},
        {"entity_groups", json::array({
            {{"name", "Conductor"}, {"dim", 2}, {"attribute_ids", {1, 2}}},
            {{"name", "FarField"}, {"dim", 1}, {"attribute_ids", {3}}}
        })},
        {"regions", json::array({
            {{"entity_group", "Conductor"}, {"material", "Copper"},
             {"current_constraint", "open"}}
        })},
        {"materials", json::array({
            {{"name", "Copper"},
             {"properties", {{"sigma", 5.8e7}, {"epsilon_r", 2.5}, {"mu_r", 1.2}}}}
        })},
        {"terminals", json::array({
            {{"name", "Coil"}, {"excitation_type", "current"},
             {"conductor_type", "stranded"}, {"entity_group", "Conductor"}}
        })},
        {"boundaries", json::array({
            {{"type", "Robin"}, {"entity_group", "FarField"},
             {"value", 4.0}, {"robin_coefficient", 2.0}}
        })},
        {"scenarios", json::array({
            {{"name", "Second"}, {"frequency", 50.0}, {"excitations", json::array({
                {{"terminal", "Coil"}, {"value", 20.0}}
            })}},
            {{"name", "First"}, {"frequency", 50.0}, {"excitations", json::array({
                {{"terminal", "Coil"}, {"value", 10.0}, {"floating", false}}
            })}}
        })}
    };
}

} // namespace

TEST_CASE("InputParser decodes the canonical schema", "[input_parser]") {
    json source = CanonicalConfig();
    const ProblemConfig config = InputParser(source).GetProblemConfig();

    REQUIRE(config.PhysicsType == PhysicsType::Magnetoquasistatics);
    REQUIRE(config.GeometryType == GeometryType::Axisymmetric);
    REQUIRE(config.AnalysisType == AnalysisType::CouplingMatrix);
    REQUIRE(config.Order == 2);
    REQUIRE(config.SolverTolerance == Catch::Approx(1e-9));
    REQUIRE(config.SolverMaxIter == 321);
    REQUIRE(config.SolverPrintLevel == 2);
    REQUIRE(config.OutputParaview);
    REQUIRE(config.OutputGmsh);
    REQUIRE(config.ExportRefine == 3);
    REQUIRE(config.Amr.Enabled);
    REQUIRE(config.Amr.MaxIterations == 4);
    REQUIRE(config.Amr.MaxDofs == 12345);
    REQUIRE(config.Amr.ErrorFraction == Catch::Approx(0.6));
    REQUIRE(config.Amr.ErrorTolerance == Catch::Approx(1e-5));

    REQUIRE(config.EntityGroups.at("Conductor").Dim == EntityDim::Domain);
    REQUIRE((config.EntityGroups.at("Conductor").AttributeIds == std::vector<int>{1, 2}));
    REQUIRE(config.EntityGroups.at("FarField").Dim == EntityDim::Boundary);
    REQUIRE(config.Regions[0].CurrentConstraint == RegionCurrentConstraint::Open);
    REQUIRE(config.Regions.size() == 1);
    REQUIRE(config.Regions[0].EntityGroupName == "Conductor");
    REQUIRE(config.Regions[0].MaterialName == "Copper");
    REQUIRE(config.Materials.at("Copper").Conductivity == Catch::Approx(5.8e7));
    REQUIRE(config.Materials.at("Copper").RelPermittivity == Catch::Approx(2.5));
    REQUIRE(config.Materials.at("Copper").RelPermeability == Catch::Approx(1.2));

    REQUIRE(config.Terminals.at("Coil").ExcitationType == Quantity::Current);
    REQUIRE(config.Terminals.at("Coil").Conductor == ConductorType::Stranded);
    REQUIRE(config.Terminals.at("Coil").EntityGroupName == "Conductor");
    REQUIRE(config.BoundaryConditions.size() == 1);
    REQUIRE(config.BoundaryConditions[0].Type == BoundaryConditionType::Robin);
    REQUIRE(config.BoundaryConditions[0].Value == Catch::Approx(4.0));
    REQUIRE(config.BoundaryConditions[0].RobinCoeff == Catch::Approx(2.0));

    REQUIRE(config.Scenarios.size() == 2);
    REQUIRE(config.Scenarios[0].first == "Second");
    REQUIRE(config.Scenarios[0].second.Frequency == Catch::Approx(50.0));
    REQUIRE(config.Scenarios[0].second.Excitations[0].Value == Catch::Approx(20.0));
    REQUIRE(config.Scenarios[1].first == "First");
}

TEST_CASE("InputParser decodes typed Dirichlet and Neumann boundaries",
          "[input_parser][boundaries]") {
    json source = CanonicalConfig();
    source["boundaries"] = json::array({
        {{"type", "Dirichlet"}, {"entity_group", "FarField"}, {"value", 3.0}},
        {{"type", "Neumann"}, {"entity_group", "FarField"}, {"value", -2.0}}
    });

    const auto boundaries = InputParser(source).GetProblemConfig().BoundaryConditions;
    REQUIRE(boundaries.size() == 2);
    REQUIRE(boundaries[0].Type == BoundaryConditionType::Dirichlet);
    REQUIRE(boundaries[1].Type == BoundaryConditionType::Neumann);
    REQUIRE(boundaries[1].Value == Catch::Approx(-2.0));
}

TEST_CASE("InputParser expands MQS frequency sweeps", "[input_parser][mqs]") {
    json source = CanonicalConfig();

    SECTION("linear spacing is inclusive and preserves excitations") {
        source["scenarios"] = json::array({{
            {"name", "Linear"},
            {"frequency", {{"scale", "linear"}, {"start", 10.0},
                           {"stop", 30.0}, {"points", 3}}},
            {"excitations", json::array({{{"terminal", "Coil"}, {"value", 7.0}}})}
        }});
        const auto scenarios = InputParser(source).GetProblemConfig().Scenarios;
        REQUIRE(scenarios.size() == 3);
        REQUIRE(scenarios[0].first == "Linear_f1_10Hz");
        REQUIRE(scenarios[1].second.Frequency == Catch::Approx(20.0));
        REQUIRE(scenarios[2].second.Frequency == Catch::Approx(30.0));
        REQUIRE(scenarios[2].second.Excitations[0].Value == Catch::Approx(7.0));
    }

    SECTION("log spacing is inclusive") {
        source["scenarios"] = json::array({{
            {"name", "Log"},
            {"frequency", {{"scale", "log"}, {"start", 10.0},
                           {"stop", 1000.0}, {"points", 3}}}
        }});
        const auto scenarios = InputParser(source).GetProblemConfig().Scenarios;
        REQUIRE(scenarios.size() == 3);
        REQUIRE(scenarios[0].second.Frequency == Catch::Approx(10.0));
        REQUIRE(scenarios[1].second.Frequency == Catch::Approx(100.0));
        REQUIRE(scenarios[2].second.Frequency == Catch::Approx(1000.0));
    }

    SECTION("one point resolves to start") {
        source["scenarios"] = json::array({{
            {"name", "Single"},
            {"frequency", {{"scale", "linear"}, {"start", 25.0},
                           {"stop", 50.0}, {"points", 1}}}
        }});
        const auto scenarios = InputParser(source).GetProblemConfig().Scenarios;
        REQUIRE(scenarios.size() == 1);
        REQUIRE(scenarios[0].second.Frequency == Catch::Approx(25.0));
        REQUIRE(scenarios[0].first == "Single_f1_25Hz");
    }
}

TEST_CASE("InputParser throws on missing file", "[input_parser]") {
    REQUIRE_THROWS_AS(InputParser(std::string("nonexistent_file.json")), std::runtime_error);
}

TEST_CASE("InputParser resolves paths relative to the config file", "[input_parser]") {
    const fs::path directory = fs::temp_directory_path() / "mfem-electromag-parser-test";
    const fs::path config_path = directory / "config.json";
    fs::create_directories(directory);

    json source = CanonicalConfig();
    source["simulation"]["results_path"] = "results";
    {
        std::ofstream output(config_path);
        output << source;
    }

    const ProblemConfig config = InputParser(config_path.string()).GetProblemConfig();
    REQUIRE(fs::path(config.MeshPath) == directory / "test.msh");
    REQUIRE(fs::path(config.ResultsDirectory) == directory / "results");

    fs::remove_all(directory);
}

TEST_CASE("Results path configures an output directory", "[input_parser]") {
    json test_config = {
        {"simulation", {
            {"physics_type", "electrostatics"},
            {"mesh", "test.msh"}
        }}
    };

    SECTION("missing path preserves mesh-directory output") {
        InputParser parser(test_config);
        REQUIRE(parser.GetProblemConfig().ResultsDirectory.empty());
    }

    SECTION("empty path preserves mesh-directory output") {
        test_config["simulation"]["results_path"] = "";
        InputParser parser(test_config);
        REQUIRE(parser.GetProblemConfig().ResultsDirectory.empty());
    }

    SECTION("relative path is resolved from the config directory") {
        test_config["simulation"]["results_path"] = "results";
        InputParser parser(test_config);
        REQUIRE(fs::path(parser.GetProblemConfig().ResultsDirectory) == fs::path(".") / "results");
    }

    SECTION("absolute path is preserved") {
        const fs::path results_path = fs::temp_directory_path() / "mfem-electromag-results";
        test_config["simulation"]["results_path"] = results_path.string();
        InputParser parser(test_config);
        REQUIRE(fs::path(parser.GetProblemConfig().ResultsDirectory) == results_path);
    }
}

TEST_CASE("InputParser wraps decoding type failures", "[input_parser]") {
    json source = CanonicalConfig();
    source["simulation"]["order"] = "second";
    REQUIRE_THROWS_AS(InputParser(source).GetProblemConfig(), std::runtime_error);
}

TEST_CASE("InputParser requires terminal excitation_type", "[input_parser]") {
    SECTION("missing key is rejected") {
        json source = CanonicalConfig();
        source["terminals"][0].erase("excitation_type");
        REQUIRE_THROWS_AS(InputParser(source).GetProblemConfig(), std::runtime_error);
    }

    SECTION("null value is rejected") {
        json source = CanonicalConfig();
        source["terminals"][0]["excitation_type"] = nullptr;
        REQUIRE_THROWS_AS(InputParser(source).GetProblemConfig(), std::runtime_error);
    }

    SECTION("the legacy 'excitation' spelling does not satisfy the requirement") {
        json source = CanonicalConfig();
        source["terminals"][0].erase("excitation_type");
        source["terminals"][0]["excitation"] = "current";
        REQUIRE_THROWS_AS(InputParser(source).GetProblemConfig(), std::runtime_error);
    }

    SECTION("the error names the offending terminal") {
        json source = CanonicalConfig();
        source["terminals"][0].erase("excitation_type");
        REQUIRE_THROWS_WITH(InputParser(source).GetProblemConfig(),
                            Catch::Matchers::ContainsSubstring("Coil") &&
                            Catch::Matchers::ContainsSubstring("excitation_type"));
    }
}
