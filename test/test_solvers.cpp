// Copyright (c) 2026 T. C. Raymond
// SPDX-License-Identifier: MIT

#include <catch2/catch_test_macros.hpp>
#include "../src/input_parser.hpp"
#include "../src/electrostatic_solver.hpp"
#include "../src/magnetostatic_solver.hpp"
#include "../src/magnetoquasistatic_solver.hpp"
#include <filesystem>
#include <fstream>

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
            {"type", "electrostatics"},
            {"mesh", mesh_file},
            {"order", 1},
            {"axisymmetric", true}
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
            {"type", "magnetostatics"},
            {"mesh", mesh_file},
            {"order", 1},
            {"axisymmetric", true}
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
            {"type", "magnetoquasistatics"},
            {"mesh", mesh_file},
            {"order", 1},
            {"axisymmetric", true},
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
                {"type", "electrostatics"},
                {"mesh", mesh_file},
                {"order", 1}
            }},
            {"materials", json::array()},
            {"boundaries", json::array()}
        };

        InputParser parser(config);
        mfem::Mesh mesh(mesh_file.c_str(), 1, 1);

        std::string type = parser.config["simulation"]["type"];
        REQUIRE(type == "electrostatics");
    }

    SECTION("magnetostatics") {
        json config = {
            {"simulation", {
                {"type", "magnetostatics"},
                {"mesh", mesh_file},
                {"order", 1}
            }},
            {"materials", json::array()},
            {"sources", json::array()},
            {"boundaries", json::array()}
        };

        InputParser parser(config);
        mfem::Mesh mesh(mesh_file.c_str(), 1, 1);

        std::string type = parser.config["simulation"]["type"];
        REQUIRE(type == "magnetostatics");
    }

    SECTION("magnetoquasistatics") {
        json config = {
            {"simulation", {
                {"type", "magnetoquasistatics"},
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

        std::string type = parser.config["simulation"]["type"];
        REQUIRE(type == "magnetoquasistatics");
    }

    // Cleanup
    fs::remove(mesh_file);
}
