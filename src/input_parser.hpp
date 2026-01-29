// Copyright (c) 2026 T. C. Raymond
// SPDX-License-Identifier: MIT

#pragma once
#include "mfem.hpp"
#include "json.hpp" // nlohmann/json
#include <fstream>
#include <iostream>
#include <unordered_map>
#include <filesystem> // C++17

using json = nlohmann::json;
using namespace mfem;
using namespace std;
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
    explicit InputParser(json c) : config(c), config_dir(".") {}
    
    std::string GetMeshPath() {
        if (config["simulation"].contains("mesh")) {
            std::string mesh_path = config["simulation"]["mesh"];
            fs::path p(mesh_path);
            
            // If path is absolute, use it directly
            if (p.is_absolute()) {
                return mesh_path;
            }

            // Otherwise, make it relative to the config file location
            return (fs::path(config_dir) / p).string();
        }
        return "default.msh";
    }

    // --------------------------------------------------------
    // Material Setup (Reluctivity nu = 1/mu) for Magnetostatics
    // --------------------------------------------------------
    void SetupReluctivity(Mesh &mesh, Vector &nu_values) {
        int max_attr = mesh.attributes.Max();
        // PWConstCoefficient expects index = attr - 1
        // So we need exactly max_attr entries.
        nu_values.SetSize(max_attr); 
        nu_values = 1.0 / (4.0 * M_PI * 1e-7); // Default nu0 (air)

        if (config.contains("materials")) {
            for (auto &mat : config["materials"]) {
                if (mat.contains("properties") && mat["properties"].contains("mu_r")) {
                    double mu_r = mat["properties"]["mu_r"];
                    double mu = mu_r * 4.0 * M_PI * 1e-7;
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
    void SetupConductivity(Mesh &mesh, Vector &sigma_values) {
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
    void SetupPermittivity(Mesh &mesh, Vector &eps_values) {
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
    void SetupBoundaries(Mesh &mesh, 
                         std::vector<std::pair<Array<int>, double>> &dirichlet_bcs) {
        
        int max_bdr = mesh.bdr_attributes.Max();
        
        if (config.contains("boundaries")) {
            for (auto &bc : config["boundaries"]) {
                if (bc["type"] == "Dirichlet") {
                    Array<int> marker(max_bdr);
                    marker = 0;
                    double val = bc["value"];

                    for (int attr : bc["attributes"]) {
                        if (attr <= max_bdr) marker[attr - 1] = 1; 
                    }
                    dirichlet_bcs.push_back({marker, val});
                }
            }
        }
    }

    // --------------------------------------------------------
    // Sources (Current Density J)
    // --------------------------------------------------------
    void SetupSources(Mesh &mesh, Vector &j_values) {
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