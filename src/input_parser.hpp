// Copyright (c) 2026 T. C. Raymond
// SPDX-License-Identifier: MIT

#pragma once
#include "json.hpp" // nlohmann/json
#include "constants.hpp"
#include "problem_config.hpp"
#include "status_reporter.hpp"
#include <fstream>
#include <iostream>
#include <unordered_map>
#include <optional>
#include <memory>
#include <filesystem> // C++17
#include <iterator>
#include <cctype>
#include <cmath>
#include <iomanip>
#include <sstream>
#include <stdexcept>

using json = nlohmann::json;
namespace fs = std::filesystem;

class InputParser {
    std::unique_ptr<json> owned_config;

public:
    const json& config;
    std::string config_dir;
    std::string config_path;   // Source file (empty when built from in-memory json)

    InputParser(const std::string &filename)
        : owned_config(std::make_unique<json>()),
          config(*owned_config) {
        std::ifstream f(filename);
        if (!f.is_open()) {
            throw std::runtime_error("Could not open config file: " + filename);
        }
        f >> *owned_config;

        // Store directory and path of config file (path used for diagnostics)
        config_path = filename;
        fs::path p(filename);
        config_dir = p.parent_path().string();
        if (config_dir.empty()) config_dir = ".";
    }

    // Allow construction from existing json
    InputParser(const json& c) : config(c), config_dir(".") {}
    InputParser(json&&) = delete;

	const ProblemConfig GetProblemConfig() const {
		try {
			return BuildProblemConfig();
		} catch (const json::exception& e) {
			// Re-throw with the config file path and, when available, a
			// resolved line/column so users can find the offending field.
			throw std::runtime_error(DescribeJsonError(e));
		}
	}

private:
	const ProblemConfig BuildProblemConfig() const {
		ProblemConfig prob_config;
		prob_config.Order = GetOrder();
		prob_config.PhysicsType = GetPhysicsType();
		prob_config.GeometryType = GetGeometryType();
		prob_config.AnalysisType = GetAnalysisType();
		prob_config.MeshPath = GetMeshPath();

		prob_config.EntityGroups = GetEntityGroups();
		prob_config.Regions = GetRegions();
		prob_config.Materials = GetMaterials();
		prob_config.Terminals = GetTerminals();
		prob_config.BoundaryConditions = GetBoundaries();
		prob_config.Scenarios = GetScenarios();

		prob_config.SolverTolerance = GetSolverTolerance();
		prob_config.SolverMaxIter = GetSolverMaxIter();
		prob_config.SolverPrintLevel = GetSolverPrintLevel();

		prob_config.OutputParaview = GetOutputParaview();
		prob_config.OutputGmsh = GetOutputGmsh();
        prob_config.ResultsDirectory = GetResultsDirectory();
		prob_config.ExportRefine = GetExportRefine();
		prob_config.Amr = GetAmrSettings();

		return prob_config;
	}

	// Build a human-friendly message from a nlohmann json exception.
	// With JSON_DIAGNOSTICS the message carries a JSON pointer like
	// "(/terminals/0/entity_group)" and with JSON_DIAGNOSTIC_POSITIONS it also
	// carries "(bytes X-Y)"; the byte range is translated into a line/column
	// against the source config file when one is available.
	[[nodiscard]] std::string DescribeJsonError(const json::exception& e) const {
		std::string what = e.what();
		std::string location;

		// Extract a "(bytes X-Y)" range emitted by JSON_DIAGNOSTIC_POSITIONS.
		const std::string marker = "(bytes ";
		const auto start = what.find(marker);
		if (start != std::string::npos) {
			const auto open = start + marker.size();
			const auto dash = what.find('-', open);
			if (dash != std::string::npos) {
				try {
					const std::size_t byte_off =
						static_cast<std::size_t>(std::stoull(what.substr(open, dash - open)));
					auto lc = ByteOffsetToLineCol(byte_off);
					if (lc.first > 0) {
						location = " (line " + std::to_string(lc.first) +
								   ", column " + std::to_string(lc.second) + ")";
					}
				} catch (const std::exception&) {
					// Leave location empty if the byte offset can't be parsed.
				}
			}
		}

		std::string prefix = "Failed to parse config";
		if (!config_path.empty()) {
			prefix += " '" + config_path + "'";
		}
		if (!location.empty()) {
			prefix += location;
		}
		return prefix + ": " + what;
	}

	// Translate a 0-based byte offset in the config file into a 1-based
	// (line, column). Returns {0, 0} when the file can't be read or the offset
	// is out of range (e.g. for in-memory configs with no backing file).
	// The file is read in text mode (no std::ios::binary) so CRLF is collapsed
	// to LF, matching the coordinate system nlohmann used while parsing (the
	// constructor reads via a text-mode ifstream); otherwise the column would
	// be skewed by one byte per preceding CRLF line ending on Windows.
	[[nodiscard]] std::pair<int, int> ByteOffsetToLineCol(std::size_t byte_off) const {
		if (config_path.empty()) {
			return {0, 0};
		}
		std::ifstream f(config_path);
		if (!f.is_open()) {
			return {0, 0};
		}
		const std::string contents((std::istreambuf_iterator<char>(f)),
								   std::istreambuf_iterator<char>());
		if (byte_off >= contents.size()) {
			byte_off = contents.empty() ? 0 : contents.size() - 1;
		}
		int line = 1;
		int col = 1;
		for (std::size_t i = 0; i < byte_off; ++i) {
			if (contents[i] == '\n') {
				++line;
				col = 1;
			} else {
				++col;
			}
		}
		return {line, col};
	}

    // "simulation.physics": electrostatics | magnetostatics | magnetoquasistatics.
    // Tolerant default; ConfigValidator enforces presence and validity.
    [[nodiscard]] ::PhysicsType GetPhysicsType() const {
        if (config.contains("simulation") && config["simulation"].is_object()) {
            const auto& sim = config["simulation"];
            if (sim.contains("physics_type") && sim["physics_type"].is_string()) {
                const std::string s = sim["physics_type"];
                if (s == "electrostatics")      return ::PhysicsType::Electrostatics;
                if (s == "magnetostatics")      return ::PhysicsType::Magnetostatics;
                if (s == "magnetoquasistatics") return ::PhysicsType::Magnetoquasistatics;
            }
        }
        return ::PhysicsType::Electrostatics; // Default
    }

    [[nodiscard]] ::GeometryType GetGeometryType() const {
        if (config.contains("simulation") && config["simulation"].is_object()) {
            const auto& sim = config["simulation"];

            // "geometry_type": "axisymmetric" | "planar"
            if (sim.contains("geometry_type") && sim["geometry_type"].is_string()) {
                std::string type_str = sim["geometry_type"];
                if (type_str == "axisymmetric") {
                    return ::GeometryType::Axisymmetric;
                } else if (type_str == "planar") {
                    return ::GeometryType::Planar;
                }
            }
        }
        return ::GeometryType::Planar; // Default
    }

    [[nodiscard]] ::AnalysisType GetAnalysisType() const {
        if (config.contains("simulation") && config["simulation"].is_object()) {
            const auto& sim = config["simulation"];
            if (sim.contains("analysis_type") && sim["analysis_type"].is_string()) {
                std::string s = sim["analysis_type"];
                if (s == "field")           return ::AnalysisType::Field;
                if (s == "coupling_matrix") return ::AnalysisType::CouplingMatrix;
            }
        }
        return ::AnalysisType::Field; // Default
    }

    [[nodiscard]] int GetOrder() const {
        if (config.contains("simulation") && config["simulation"].is_object() &&
            config["simulation"].contains("order")) {
            int order = config["simulation"]["order"];
            if (order < 1) {
                StatusReporter::Global().Warning(
                    "'order' must be >= 1, got " + std::to_string(order)
                    + ". Using order = 1.");
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

    [[nodiscard]] std::string GetResultsDirectory() const {
        if (config.contains("simulation") && config["simulation"].is_object() &&
            config["simulation"].contains("results_path")) {
            std::string p = config["simulation"]["results_path"];
            if (p.empty()) return {};
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

    // Parse the optional "simulation.amr" block. Missing keys fall back to the
    // AmrSettings defaults and unknown future keys are ignored, so the two sides
    // (C# requester / solver) can evolve independently. nlohmann::json::value()
    // parses numbers with the invariant '.' separator regardless of locale.
    [[nodiscard]] AmrSettings GetAmrSettings() const {
        AmrSettings amr; // defaults (disabled)
        if (config.contains("simulation") && config["simulation"].is_object() &&
            config["simulation"].contains("amr") &&
            config["simulation"]["amr"].is_object()) {
            const auto& a = config["simulation"]["amr"];
            amr.Enabled        = a.value("enabled",         amr.Enabled);
            amr.MaxIterations  = a.value("max_iterations",  amr.MaxIterations);
            amr.MaxDofs        = a.value("max_dofs",         amr.MaxDofs);
            amr.ErrorFraction  = a.value("error_fraction",  amr.ErrorFraction);
            amr.ErrorTolerance = a.value("error_tolerance", amr.ErrorTolerance);
            amr.Conforming     = a.value("conforming",       amr.Conforming);
        }
        return amr;
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

	std::unordered_map<std::string, EntityGroup> GetEntityGroups() const {
		std::unordered_map<std::string, EntityGroup> groups;
		if (config.contains("entity_groups")) {
			for (const auto& g : config["entity_groups"]) {
                EntityGroup group;
				std::string name = g.value("name", std::string{});
				if (g.contains("dim")) {
					int dim = g["dim"];
					group.Dim = (dim == 1) ? EntityDim::Boundary : EntityDim::Domain;
				}
				if (g.contains("attribute_ids")) {
					for (int attr : g["attribute_ids"]) {
						group.AttributeIds.push_back(attr);
					}
				}
				groups.emplace(std::move(name), std::move(group));
			}
		}
		return groups;
	}

    std::map<std::string, Terminal> GetTerminals() const {
        std::map<std::string, Terminal> terminals;
        if (config.contains("terminals")) {
            for (auto& t : config["terminals"]) {
                Terminal terminal;
				std::string name = t.value("name", std::string{});
                // "excitation": "voltage" (default) | "current"
                if (t.contains("excitation") && t["excitation"].is_string()) {
                    std::string d = t["excitation"];
                    terminal.Excitation = (d == "current") ? Quantity::Current
                                                           : Quantity::Voltage;
                }
				if (t.contains("conductor_type") && t["conductor_type"].is_string()) {
					std::string c = t["conductor_type"];
					terminal.Conductor = (c == "stranded") ? ConductorType::Stranded
						: ConductorType::Massive;
				}
                if (t.contains("entity_group")) {
                    terminal.EntityGroupName = t["entity_group"];
                }
                terminals.emplace(std::move(name), std::move(terminal));
            }
        }
        return terminals;
    }

    std::vector<Region> GetRegions() const {
        std::vector<Region> regions;
        if (config.contains("regions")) {
            for (auto &region : config["regions"]) {
                Region _region;
                if (region.contains("entity_group")) {
                    _region.EntityGroupName = region["entity_group"];
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
                if (!bc.contains("type") || !bc.contains("value") || !bc.contains("entity_group")) {
                     continue; // Skip invalid entries, let validator handle reporting
                }
                std::string bc_type = bc["type"];
                std::string group_name = bc["entity_group"];
                double val = bc["value"];
                double robin_coeff = bc.value("robin_coefficient", 1.0);

                bcs.emplace_back(bc_type, group_name, val, robin_coeff);
            }
        }
        return bcs;
    }

    // --------------------------------------------------------
    // Scenarios (one solve each: parameters + per-terminal excitations)
    // --------------------------------------------------------
    static std::string SweepScenarioName(const std::string& base_name,
                                         int point,
                                         double frequency) {
        std::ostringstream value;
        value << std::setprecision(12) << frequency;
        std::string token = value.str();
        for (char& c : token) {
            if (c == '.') c = 'p';
            else if (c == '+') c = '_';
            else if (c == '-') c = 'm';
        }
        return base_name + "_f" + std::to_string(point + 1) + "_" + token + "Hz";
    }

    static Scenario ParseScenarioExcitations(const json& source) {
        Scenario scenario;
        if (source.contains("excitations")) {
            for (const auto& d : source["excitations"]) {
                Excitation excitation;
                if (d.contains("terminal")) {
                    excitation.TerminalName = d["terminal"];
                }
                excitation.Value = d.value("value", 0.0);
                excitation.Floating = d.value("floating", false);
                scenario.Excitations.push_back(excitation);
            }
        }
        return scenario;
    }

    std::vector<std::pair<std::string, Scenario>> GetScenarios() const {
        std::vector<std::pair<std::string, Scenario>> scenarios;

        if (config.contains("scenarios")) {
            const bool is_mqs = GetPhysicsType() == PhysicsType::Magnetoquasistatics;
            for (const auto& sc : config["scenarios"]) {
                const std::string name = sc.value("name", "");
                Scenario scenario = ParseScenarioExcitations(sc);

                if (!is_mqs || !sc.contains("frequency")) {
                    scenarios.emplace_back(name, std::move(scenario));
                    continue;
                }

                const auto& frequency = sc["frequency"];
                if (frequency.is_number()) {
                    scenario.Frequency = frequency.get<double>();
                    scenarios.emplace_back(name, std::move(scenario));
                    continue;
                }

                const std::string scale = frequency.at("scale").get<std::string>();
                const double start = frequency.at("start").get<double>();
                const double stop = frequency.at("stop").get<double>();
                const int points = frequency.at("points").get<int>();
                const double log_start = std::log(start);
                const double log_stop = std::log(stop);

                for (int i = 0; i < points; ++i) {
                    Scenario point_scenario = scenario;
                    if (points == 1) {
                        point_scenario.Frequency = start;
                    }
                    else if (i == points - 1) {
                        point_scenario.Frequency = stop;
                    }
                    else {
                        const double fraction = static_cast<double>(i) / (points - 1);
                        point_scenario.Frequency = scale == "log"
                            ? std::exp(log_start + fraction * (log_stop - log_start))
                            : start + fraction * (stop - start);
                    }
                    scenarios.emplace_back(
                        SweepScenarioName(name, i, point_scenario.Frequency),
                        std::move(point_scenario));
                }
            }
        }
        return scenarios;
    }
};