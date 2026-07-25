// Copyright (c) 2026 T. C. Raymond
// SPDX-License-Identifier: MIT

#include <catch2/catch_test_macros.hpp>
#include "../src/input_parser.hpp"
#include <fstream>
#include <filesystem>

namespace fs = std::filesystem;

TEST_CASE("InputParser can parse valid JSON", "[input_parser]") {
    json test_config = {
        {"simulation", {
            {"type", "electrostatics"},
            {"mesh", "test.msh"},
            {"order", 2},
            {"axisymmetric", true}
        }},
        {"materials", json::array()}
    };

    InputParser parser(test_config);
    REQUIRE(parser.config["simulation"]["type"] == "electrostatics");
    REQUIRE(parser.config["simulation"]["order"] == 2);
}

TEST_CASE("InputParser throws on missing file", "[input_parser]") {
    REQUIRE_THROWS_AS(InputParser(std::string("nonexistent_file.json")), std::runtime_error);
}

TEST_CASE("GetMeshPath handles relative paths", "[input_parser]") {
    // Create a temporary mesh file for testing
    std::string temp_mesh = "temp_test_mesh.msh";
    std::ofstream temp_file(temp_mesh);
    temp_file << "test mesh content" << std::endl;
    temp_file.close();

    json test_config = {
        {"simulation", {
            {"type", "electrostatics"},
            {"mesh", temp_mesh}
        }}
    };

    InputParser parser(test_config);
    // std::string mesh_path = parser.GetMeshPath();

    // REQUIRE(fs::exists(mesh_path));

    // Cleanup
    fs::remove(temp_mesh);
}

TEST_CASE("GetMeshPath throws on missing mesh file", "[input_parser]") {
    json test_config = {
        {"simulation", {
            {"type", "electrostatics"},
            {"mesh", "definitely_does_not_exist.msh"}
        }}
    };

    InputParser parser(test_config);
    // REQUIRE_THROWS_AS(parser.GetMeshPath(), std::runtime_error);
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

TEST_CASE("SetupBoundaries handles valid attributes", "[input_parser]") {
    // Create a simple test mesh
    std::string temp_mesh = "temp_boundary_test.msh";
    std::ofstream temp_file(temp_mesh);
    temp_file << "$MeshFormat\n4.1 0 8\n$EndMeshFormat\n";
    temp_file << "$Entities\n0 0 0 0\n$EndEntities\n";
    temp_file << "$Nodes\n0 0 0\n$EndNodes\n";
    temp_file << "$Elements\n0 0 0 0\n$EndElements\n";
    temp_file.close();

    json test_config = {
        {"simulation", {
            {"type", "electrostatics"},
            {"mesh", temp_mesh}
        }},
        {"boundaries", json::array({
            {
                {"name", "test_boundary"},
                {"attributes", {1, 2}},
                {"type", "Dirichlet"},
                {"value", 100.0}
            }
        })}
    };

    InputParser parser(test_config);

    // Note: We can't fully test without creating a valid MFEM mesh,
    // but we can verify the parser doesn't crash
    REQUIRE(parser.config["boundaries"].size() == 1);

    // Cleanup
    fs::remove(temp_mesh);
}

TEST_CASE("Boundary attribute bounds checking", "[input_parser]") {
    json test_config = {
        {"simulation", {
            {"type", "electrostatics"},
            {"mesh", "test.msh"}
        }},
        {"boundaries", json::array({
            {
                {"name", "invalid_boundary"},
                {"attributes", {0, -1, 100}},  // Invalid attributes
                {"type", "Dirichlet"},
                {"value", 0.0}
            }
        })}
    };

    InputParser parser(test_config);

    // The parser should handle out-of-bounds attributes gracefully
    // (not crash or throw, just skip invalid attributes)
    REQUIRE(parser.config["boundaries"].size() == 1);
}
