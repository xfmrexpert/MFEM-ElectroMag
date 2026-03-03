// Copyright (c) 2026 T. C. Raymond
// SPDX-License-Identifier: MIT

#pragma once
#include "json.hpp" // nlohmann/json
#include "constants.hpp"
#include "problem_config.hpp"
#include <fstream>
#include <iostream>
#include <unordered_map>
#include <optional>
#include <memory>
#include <filesystem> // C++17

using json = nlohmann::json;
namespace fs = std::filesystem;

class InputParser {
    std::unique_ptr<json> owned_config;

public:
    const json& config;
    std::string config_dir;

    InputParser(const std::string &filename)
        : owned_config(std::make_unique<json>()),
          config(*owned_config) {
        std::ifstream f(filename);
        if (!f.is_open()) {
            throw std::runtime_error("Could not open config file: " + filename);
        }
        f >> *owned_config;

        // Store directory of config file
        fs::path p(filename);
        config_dir = p.parent_path().string();
        if (config_dir.empty()) config_dir = ".";
    }

    // Allow construction from existing json
    InputParser(const json& c) : config(c), config_dir(".") {}
    InputParser(json&&) = delete;

    const ProblemConfig GetProblemConfig() const {
        ProblemConfig prob_config;
        prob_config.SolverTolerance = GetSolverTolerance();
        prob_config.SolverMaxIter = GetSolverMaxIter();
        prob_config.SolverPrintLevel = GetSolverPrintLevel();
        prob_config.ModelType = GetModelType();
        prob_config.MeshPath = GetMeshPath();
        prob_config.Ports = GetPorts();
        prob_config.Regions = GetRegions();
        prob_config.Materials = GetMaterials();
        prob_config.BoundaryConditions = GetBoundaries();
        prob_config.Sources = GetSources();
        return prob_config;
    }

private:
    [[nodiscard]] ::ModelType GetModelType() const {
        if (config.contains("simulation") && config["simulation"].is_object() &&
            config["simulation"].contains("model_type")) {
            std::string type_str = config["simulation"]["model_type"];
            if (type_str == "axisymmetric") {
                return ::ModelType::Axisymmetric;
            } else if (type_str == "planar") {
                return ::ModelType::Planar; 
            }
        }
        return ::ModelType::Planar; // Default
    }

// Get solver parameters with defaults
    [[nodiscard]] double GetSolverTolerance() const {
        if (config.contains("simulation") && config["simulation"].is_object() &&
            config["simulation"].contains("solver_tolerance")) {
            return config["simulation"]["solver_tolerance"];
        }
        return Constants::DEFAULT_SOLVER_TOLERANCE;
    }

    [[nodiscard]] int GetSolverMaxIter() const {
        if (config.contains("simulation") && config["simulation"].is_object() &&
            config["simulation"].contains("solver_max_iter")) {
            return config["simulation"]["solver_max_iter"];
        }
        return Constants::DEFAULT_SOLVER_MAX_ITER;
    }

    [[nodiscard]] int GetSolverPrintLevel() const {
        if (config.contains("simulation") && config["simulation"].is_object() &&
            config["simulation"].contains("solver_print_level")) {
            return config["simulation"]["solver_print_level"];
        }
        return Constants::DEFAULT_SOLVER_PRINT_LEVEL;
    }
    
    [[nodiscard]] std::string GetMeshPath() const {
        std::string mesh_path;

        if (config.contains("simulation") && config["simulation"].is_object() &&
            config["simulation"].contains("mesh")) {
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

        return mesh_path;
    }

    std::vector<Port> GetPorts() const {
        std::vector<Port> ports;
        if (config.contains("ports")) {
            for (auto &port : config["ports"]) {
                if (port.contains("region")) {
                    int region = port["region"];
                    ports.push_back({region});
                }
            }
        }
        return ports;
    }

    std::vector<Region> GetRegions() const {
        std::vector<Region> regions;
        if (config.contains("regions")) {
            for (auto &region : config["regions"]) {
                Region _region;
                if (region.contains("attribute_ids")) {
                    for (int attr : region["attribute_ids"]) {
                        _region.AttributeIds.push_back(attr);
                    }
                }
                if (region.contains("material")) {
                    _region.Material = (int)region["material"] - 1; // Convert 1-based to 0-based
                }
                regions.push_back(_region);
            }
        }
        return regions;
    }

    std::vector<Material> GetMaterials() const {
        std::vector<Material> materials;
        if (config.contains("materials")) {
            for (auto &material : config["materials"]) {
                Material _material;
                if (material.contains("properties")) {
                    auto& props = material["properties"];
                    if (props.contains("sigma")) {
                        _material.Conductivity = props["sigma"];
                    }
                    if (props.contains("epsilon_r")) {
                        _material.RelPermittivity = props["epsilon_r"];
                    }
                    if (props.contains("mu_r")) {
                        _material.RelPermeability = props["mu_r"];
                    }
                }
                materials.push_back(_material);
            }
        }
        return materials;
    }

    // --------------------------------------------------------
    // Boundary Conditions
    // --------------------------------------------------------

    std::vector<BoundaryCondition> GetBoundaries() const {
        
        std::vector<BoundaryCondition> bcs;

        if (config.contains("boundaries")) {
            for (auto &bc : config["boundaries"]) {
                if (!bc.contains("type") || !bc.contains("value") || !bc.contains("attributes")) {
                     continue; // Skip invalid entries, let validator handle reporting
                }
                std::string bc_type = bc["type"];
                std::vector<int> markers;
                double val = bc["value"];
                double robin_coeff = bc.value("robin_coefficient", 1.0);

                for (int attr : bc["attributes"]) {
                    markers.push_back(attr);
                }

                bcs.emplace_back(bc_type, markers, val, robin_coeff);
            }
        }
        return bcs;
    }

    // --------------------------------------------------------
    // Sources (Current Density J)
    // --------------------------------------------------------
    std::vector<Source> GetSources() const {
        std::vector<Source> sources;

        if (config.contains("sources")) {
            for (auto &src : config["sources"]) {
                Source source;
                if (src.contains("type") && src["type"] == "CurrentDensity") {
                    if (src.contains("value")) {
                        source.CurrentDensity = src["value"];
                    }
                    if (src.contains("attributes")) {
                        for (int attr : src["attributes"]) {
                            source.Markers.push_back(attr);
                        }
                    }
                    sources.push_back(source);
                }
            }
        }
        return sources;
    }
};