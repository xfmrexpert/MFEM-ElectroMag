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
        prob_config.Order = GetOrder();
        prob_config.SolverTolerance = GetSolverTolerance();
        prob_config.SolverMaxIter = GetSolverMaxIter();
        prob_config.SolverPrintLevel = GetSolverPrintLevel();
        prob_config.ModelType = GetModelType();
        prob_config.StudyType = GetStudyType();
        prob_config.Frequency = GetFrequency();
        prob_config.MeshPath = GetMeshPath();
        prob_config.OutputParaview = GetOutputParaview();
        prob_config.OutputGmsh = GetOutputGmsh();
        prob_config.ResultsFile = GetResultsFile();
        prob_config.ExportRefine = GetExportRefine();
        prob_config.Regions = GetRegions();
        prob_config.Materials = GetMaterials();
        prob_config.Terminals = GetTerminals();
        prob_config.BoundaryConditions = GetBoundaries();
        prob_config.Scenarios = GetScenarios();
        return prob_config;
    }

private:
    [[nodiscard]] ::ModelType GetModelType() const {
        if (config.contains("simulation") && config["simulation"].is_object()) {
            const auto& sim = config["simulation"];

            // Preferred: explicit string "model_type": "axisymmetric" | "planar"
            if (sim.contains("model_type") && sim["model_type"].is_string()) {
                std::string type_str = sim["model_type"];
                if (type_str == "axisymmetric") {
                    return ::ModelType::Axisymmetric;
                } else if (type_str == "planar") {
                    return ::ModelType::Planar;
                }
            }

            // Convenience: boolean "axisymmetric": true | false
            if (sim.contains("axisymmetric") && sim["axisymmetric"].is_boolean()) {
                return sim["axisymmetric"].get<bool>()
                    ? ::ModelType::Axisymmetric
                    : ::ModelType::Planar;
            }
        }
        return ::ModelType::Planar; // Default
    }

    [[nodiscard]] ::StudyType GetStudyType() const {
        if (config.contains("simulation") && config["simulation"].is_object()) {
            const auto& sim = config["simulation"];
            if (sim.contains("study") && sim["study"].is_string()) {
                std::string s = sim["study"];
                if (s == "field")           return ::StudyType::Field;
                if (s == "coupling_matrix") return ::StudyType::CouplingMatrix;
            }
        }
        return ::StudyType::Field; // Default
    }

    [[nodiscard]] double GetFrequency() const {
        if (config.contains("simulation") && config["simulation"].is_object() &&
            config["simulation"].contains("frequency")) {
            return config["simulation"]["frequency"];
        }
        return 60.0; // Default (MQS only; ignored by ES/MS)
    }

    [[nodiscard]] int GetOrder() const {
        if (config.contains("simulation") && config["simulation"].is_object() &&
            config["simulation"].contains("order")) {
            int order = config["simulation"]["order"];
            if (order < 1) {
                std::cerr << "Warning: 'order' must be >= 1, got " << order
                          << ". Using order = 1." << std::endl;
                return 1;
            }
            return order;
        }
        return 1; // Default
    }

    [[nodiscard]] bool GetOutputParaview() const {
        if (config.contains("simulation") && config["simulation"].is_object() &&
            config["simulation"].contains("output_paraview") &&
            config["simulation"]["output_paraview"].is_boolean()) {
            return config["simulation"]["output_paraview"].get<bool>();
        }
        return false;
    }

    [[nodiscard]] bool GetOutputGmsh() const {
        if (config.contains("simulation") && config["simulation"].is_object() &&
            config["simulation"].contains("output_gmsh") &&
            config["simulation"]["output_gmsh"].is_boolean()) {
            return config["simulation"]["output_gmsh"].get<bool>();
        }
        return false;
    }

    [[nodiscard]] std::string GetResultsFile() const {
        if (config.contains("simulation") && config["simulation"].is_object() &&
            config["simulation"].contains("results_file")) {
            std::string p = config["simulation"]["results_file"];
            fs::path pp(p);
            if (pp.is_absolute()) return pp.string();
            return (fs::path(config_dir) / pp).string();
        }
        return {};
    }

    [[nodiscard]] int GetExportRefine() const {
        if (config.contains("simulation") && config["simulation"].is_object() &&
            config["simulation"].contains("export_refine")) {
            int n = config["simulation"]["export_refine"];
            if (n < 1) return 1;
            return n;
        }
        return -1; // Sentinel: caller should default to Order.
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

    std::vector<Terminal> GetTerminals() const {
        std::vector<Terminal> terminals;
        if (config.contains("terminals")) {
            for (auto &t : config["terminals"]) {
                Terminal terminal;
                if (t.contains("name")) {
                    terminal.Name = t["name"];
                }
                // "excitation": "voltage" (default) | "current"
                if (t.contains("excitation") && t["excitation"].is_string()) {
                    std::string d = t["excitation"];
                    terminal.Excitation = (d == "current") ? Quantity::Current
                                                           : Quantity::Voltage;
                }
                if (t.contains("attribute_ids")) {
                    for (int attr : t["attribute_ids"]) {
                        terminal.AttributeIds.push_back(attr);
                    }
                }
                terminals.push_back(terminal);
            }
        }
        return terminals;
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
    // Boundary Conditions (closures: far-field, symmetry, axis)
    // --------------------------------------------------------

    std::vector<BoundaryCondition> GetBoundaries() const {

        std::vector<BoundaryCondition> bcs;

        if (config.contains("boundaries")) {
            for (auto &bc : config["boundaries"]) {
                if (!bc.contains("type") || !bc.contains("value") || !bc.contains("attribute_ids")) {
                     continue; // Skip invalid entries, let validator handle reporting
                }
                std::string bc_type = bc["type"];
                std::vector<int> markers;
                double val = bc["value"];
                double robin_coeff = bc.value("robin_coefficient", 1.0);

                for (int attr : bc["attribute_ids"]) {
                    markers.push_back(attr);
                }

                bcs.emplace_back(bc_type, markers, val, robin_coeff);
            }
        }
        return bcs;
    }

    // --------------------------------------------------------
    // Scenarios (one solve each: parameters + per-terminal excitations)
    // --------------------------------------------------------
    std::vector<Scenario> GetScenarios() const {
        std::vector<Scenario> scenarios;

        if (config.contains("scenarios")) {
            for (auto &sc : config["scenarios"]) {
                Scenario scenario;
                if (sc.contains("name")) {
                    scenario.Name = sc["name"];
                }

                if (sc.contains("excitations")) {
                    for (auto &d : sc["excitations"]) {
                        Excitation excitation;
                        if (d.contains("terminal")) {
                            excitation.TerminalName = d["terminal"];
                        }
                        excitation.Value = d.value("value", 0.0);
                        excitation.Floating = d.value("floating", false);
                        scenario.Excitations.push_back(excitation);
                    }
                }
                scenarios.push_back(scenario);
            }
        }
        return scenarios;
    }
};