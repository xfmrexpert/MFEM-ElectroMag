// Copyright (c) 2026 T. C. Raymond
// SPDX-License-Identifier: MIT

#pragma once
#include "mfem.hpp"
#include "json.hpp" // nlohmann/json
#include "constants.hpp"
#include <fstream>
#include <iostream>
#include <unordered_map>
#include <optional>
#include <filesystem> // C++17

using json = nlohmann::json;
namespace fs = std::filesystem;

class InputParser {
public:
    json config;
    std::string config_dir;

    InputParser(const std::string &filename) {
        std::ifstream f(filename);
        if (!f.is_open()) {
            throw std::runtime_error("Could not open config file: " + filename);
        }
        f >> config;
        
        // Store directory of config file
        fs::path p(filename);
        config_dir = p.parent_path().string();
        if (config_dir.empty()) config_dir = ".";
    }

    // Allow construction from existing json
    InputParser(const json& c) : config(c), config_dir(".") {}

    // Get solver parameters with defaults
    [[nodiscard]] double GetSolverTolerance() const {
        if (config["simulation"].contains("solver_tolerance")) {
            return config["simulation"]["solver_tolerance"];
        }
        return Constants::DEFAULT_SOLVER_TOLERANCE;
    }

    [[nodiscard]] int GetSolverMaxIter() const {
        if (config["simulation"].contains("solver_max_iter")) {
            return config["simulation"]["solver_max_iter"];
        }
        return Constants::DEFAULT_SOLVER_MAX_ITER;
    }

    [[nodiscard]] int GetSolverPrintLevel() const {
        if (config["simulation"].contains("solver_print_level")) {
            return config["simulation"]["solver_print_level"];
        }
        return Constants::DEFAULT_SOLVER_PRINT_LEVEL;
    }
    
    [[nodiscard]] std::string GetMeshPath() {
        std::string mesh_path;

        if (config["simulation"].contains("mesh")) {
            mesh_path = config["simulation"]["mesh"];
            fs::path p(mesh_path);

            // If path is absolute, use it directly
            if (p.is_absolute()) {
                mesh_path = p.string();
            } else {
                // Otherwise, make it relative to the config file location
                mesh_path = (fs::path(config_dir) / p).string();
            }
        } else {
            mesh_path = "default.msh";
        }

        // Validate that the mesh file exists
        if (!fs::exists(mesh_path)) {
            throw std::runtime_error("Mesh file not found: " + mesh_path);
        }

        return mesh_path;
    }

    // --------------------------------------------------------
    // Material Setup (Reluctivity nu = 1/mu) for Magnetostatics
    // --------------------------------------------------------
    void SetupReluctivity(mfem::Mesh &mesh, mfem::Vector &nu_values) {
        int max_attr = mesh.attributes.Max();
        // PWConstCoefficient expects index = attr - 1
        // So we need exactly max_attr entries.
        nu_values.SetSize(max_attr);
        nu_values = 1.0 / Constants::MU_0; // Default nu0 (air)

        if (config.contains("materials")) {
            for (auto &mat : config["materials"]) {
                if (mat.contains("properties") && mat["properties"].contains("mu_r")) {
                    double mu_r = mat["properties"]["mu_r"];
                    double mu = mu_r * Constants::MU_0;
                    double nu = 1.0 / mu;

                    for (int attr : mat["attributes"]) {
                        // Protect against out-of-bounds keys
                        if (attr > 0 && attr <= max_attr) {
                            nu_values[attr - 1] = nu; // <--- FIX: attr - 1
                        }
                    }
                }
            }
        }
    }

    // --------------------------------------------------------
    // Material Setup (Conductivity) for Magnetoquasistatics
    // --------------------------------------------------------
    void SetupConductivity(mfem::Mesh &mesh, mfem::Vector &sigma_values) {
        int max_attr = mesh.attributes.Max();
        sigma_values.SetSize(max_attr);
        sigma_values = 0.0; 

        if (config.contains("materials")) {
            for (auto &mat : config["materials"]) {
                if (mat.contains("properties") && mat["properties"].contains("sigma")) {
                    double sigma = mat["properties"]["sigma"];

                    for (int attr : mat["attributes"]) {
                        if (attr > 0 && attr <= max_attr) {
                            sigma_values[attr - 1] = sigma; // <--- FIX: attr - 1
                        }
                    }
                }
            }
        }
    }

    // --------------------------------------------------------
    // Material Setup (Permittivity epsilon) for Electrostatics
    // --------------------------------------------------------
    void SetupPermittivity(mfem::Mesh &mesh, mfem::Vector &eps_values) {
        int max_attr = mesh.attributes.Max();
        eps_values.SetSize(max_attr);
        eps_values = 1.0; 

        if (config.contains("materials")) {
            for (auto &mat : config["materials"]) {
                if (mat.contains("properties") && mat["properties"].contains("epsilon_r")) {
                    double eps_r = mat["properties"]["epsilon_r"];

                    for (int attr : mat["attributes"]) {
                        if (attr > 0 && attr <= max_attr) {
                            eps_values[attr - 1] = eps_r; // <--- FIX: attr - 1
                        }
                    }
                }
            }
        }
    }

    // --------------------------------------------------------
    // Boundary Conditions
    // --------------------------------------------------------
    struct BoundaryCondition {
        std::string type;      // "Dirichlet", "Neumann", "Robin"
        mfem::Array<int> marker;
        double value;
        double robin_coeff;    // For Robin BCs: alpha * u + beta * du/dn = value

        BoundaryCondition(const std::string& t, const mfem::Array<int>& m, double v, double rc = 0.0)
            : type(t), marker(m), value(v), robin_coeff(rc) {}
    };

    void SetupBoundaries(mfem::Mesh &mesh,
                         std::vector<std::pair<mfem::Array<int>, double>> &dirichlet_bcs) {

        int max_bdr = mesh.bdr_attributes.Max();

        if (config.contains("boundaries")) {
            for (auto &bc : config["boundaries"]) {
                if (bc["type"] == "Dirichlet") {
                    mfem::Array<int> marker(max_bdr);
                    marker = 0;
                    double val = bc["value"];

                    for (int attr : bc["attributes"]) {
                        if (attr > 0 && attr <= max_bdr) marker[attr - 1] = 1;
                    }
                    dirichlet_bcs.push_back({marker, val});
                }
                // Neumann and Robin BCs are stored but not yet fully implemented in solvers
                // This provides the infrastructure for future implementation
            }
        }
    }

    // Extended version that handles all boundary types
    void SetupAllBoundaries(mfem::Mesh &mesh,
                           std::vector<BoundaryCondition> &all_bcs) {

        int max_bdr = mesh.bdr_attributes.Max();

        if (config.contains("boundaries")) {
            for (auto &bc : config["boundaries"]) {
                std::string bc_type = bc["type"];
                mfem::Array<int> marker(max_bdr);
                marker = 0;
                double val = bc["value"];
                double robin_coeff = bc.value("robin_coefficient", 1.0);

                for (int attr : bc["attributes"]) {
                    if (attr > 0 && attr <= max_bdr) marker[attr - 1] = 1;
                }

                all_bcs.emplace_back(bc_type, marker, val, robin_coeff);
            }
        }
    }

    // --------------------------------------------------------
    // Sources (Current Density J)
    // --------------------------------------------------------
    void SetupSources(mfem::Mesh &mesh, mfem::Vector &j_values) {
        int max_attr = mesh.attributes.Max();
        j_values.SetSize(max_attr);
        j_values = 0.0; 

        if (config.contains("sources")) {
            for (auto &src : config["sources"]) {
                if (src["type"] == "CurrentDensity") {
                    double val = src["value"];
                    for (int attr : src["attributes"]) {
                        if (attr > 0 && attr <= max_attr) {
                            j_values[attr - 1] = val; // <--- FIX: attr - 1
                        }
                    }
                }
            }
        }
    }
};