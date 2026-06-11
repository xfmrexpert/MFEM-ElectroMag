// Copyright (c) 2026 T. C. Raymond
// SPDX-License-Identifier: MIT

#pragma once

#include <string>
#include <vector>
#include <set>
#include <stdexcept>
#include "json.hpp"
#include "mfem.hpp"

using json = nlohmann::json;

/**
 * @brief Comprehensive configuration validator
 */
class ConfigValidator {
public:
    struct ValidationError {
        std::string field;
        std::string message;

        ValidationError(const std::string& f, const std::string& m)
            : field(f), message(m) {}
    };

private:
    std::vector<ValidationError> errors;

    void AddError(const std::string& field, const std::string& message) {
        errors.emplace_back(field, message);
    }

    void ValidateSimulation(const json& config) {
        if (!config.contains("simulation")) {
            AddError("simulation", "Missing required 'simulation' section");
            return;
        }

        const auto& sim = config["simulation"];

        // Required fields
        if (!sim.contains("type")) {
            AddError("simulation.type", "Missing required field 'type'");
        } else {
            std::string type = sim["type"];
            if (type != "electrostatics" && type != "magnetostatics" && type != "magnetoquasistatics") {
                AddError("simulation.type", "Invalid type '" + type + "'. Must be 'electrostatics', 'magnetostatics', or 'magnetoquasistatics'");
            }

            // Type-specific requirements
            if (type == "magnetoquasistatics" && !sim.contains("frequency")) {
                AddError("simulation.frequency", "Magnetoquasistatic simulations require 'frequency' field");
            }
        }

        if (!sim.contains("mesh")) {
            AddError("simulation.mesh", "Missing required field 'mesh'");
        }

        // Optional fields with validation
        if (sim.contains("order")) {
            int order = sim["order"];
            if (order < 1 || order > 10) {
                AddError("simulation.order", "Order must be between 1 and 10");
            }
        }

        if (sim.contains("solver_tolerance")) {
            double tol = sim["solver_tolerance"];
            if (tol <= 0 || tol >= 1) {
                AddError("simulation.solver_tolerance", "Solver tolerance must be between 0 and 1");
            }
        }

        if (sim.contains("solver_max_iter")) {
            int max_iter = sim["solver_max_iter"];
            if (max_iter < 1) {
                AddError("simulation.solver_max_iter", "Max iterations must be at least 1");
            }
        }
    }

    void ValidateRegions(const json& config, const mfem::Mesh* mesh = nullptr) {
        if (!config.contains("regions")) {
            return;
        }

        const auto& regions = config["regions"];
        if (!regions.is_array()) {
            AddError("regions", "Regions must be an array");
            return;
        }

        int max_attr = mesh ? mesh->attributes.Max() : 0;
        size_t num_materials = 0;
        if (config.contains("materials") && config["materials"].is_array()) {
            num_materials = config["materials"].size();
        }

        for (size_t i = 0; i < regions.size(); ++i) {
            const auto& reg = regions[i];
            std::string prefix = "regions[" + std::to_string(i) + "]";

            if (!reg.contains("attribute_ids")) {
                AddError(prefix + ".attribute_ids", "Missing required field 'attribute_ids'");
            } else if (mesh) {
                for (auto attr : reg["attribute_ids"]) {
                    int attr_val = attr;
                    if (attr_val <= 0 || attr_val > max_attr) {
                        AddError(prefix + ".attribute_id", "Attribute " + std::to_string(attr_val) +
                                " is out of range [1, " + std::to_string(max_attr) + "]");
                    }
                }
            }

            if (!reg.contains("material")) {
                AddError(prefix + ".material", "Missing required field 'material'");
            } else {
                int mat_idx = reg["material"];
                if (mat_idx < 1 || mat_idx > static_cast<int>(num_materials)) {
                    AddError(prefix + ".material", "Material index " + std::to_string(mat_idx) +
                             " is out of range [1, " + std::to_string(num_materials) + "]");
                }
            }
        }
    }

    void ValidateMaterials(const json& config, const mfem::Mesh* mesh = nullptr) {
        if (!config.contains("materials")) {
            return; // Materials are optional for some physics types
        }

        const auto& materials = config["materials"];
        if (!materials.is_array()) {
            AddError("materials", "Materials must be an array");
            return;
        }

        std::string type = config["simulation"].value("type", "");

        for (size_t i = 0; i < materials.size(); ++i) {
            const auto& mat = materials[i];
            std::string prefix = "materials[" + std::to_string(i) + "]";

            if (!mat.contains("properties")) {
                AddError(prefix + ".properties", "Missing required field 'properties'");
                continue;
            }

            const auto& props = mat["properties"];

            // Type-specific property validation
            if (type == "electrostatics") {
                if (!props.contains("epsilon_r")) {
                    AddError(prefix + ".properties.epsilon_r", "Electrostatic materials require 'epsilon_r'");
                } else {
                    double eps_r = props["epsilon_r"];
                    if (eps_r <= 0) {
                        AddError(prefix + ".properties.epsilon_r", "Permittivity must be positive");
                    }
                }
            } else if (type == "magnetostatics" || type == "magnetoquasistatics") {
                if (!props.contains("mu_r")) {
                    AddError(prefix + ".properties.mu_r", "Magnetic materials require 'mu_r'");
                } else {
                    double mu_r = props["mu_r"];
                    if (mu_r <= 0) {
                        AddError(prefix + ".properties.mu_r", "Permeability must be positive");
                    }
                }

                if (type == "magnetoquasistatics" && !props.contains("sigma")) {
                    // Sigma is optional (default 0), but warn if missing
                }
            }
        }
    }

    void ValidateBoundaries(const json& config, const mfem::Mesh* mesh = nullptr) {
        if (!config.contains("boundaries")) {
            return; // Boundaries might be optional
        }

        const auto& boundaries = config["boundaries"];
        if (!boundaries.is_array()) {
            AddError("boundaries", "Boundaries must be an array");
            return;
        }

        int max_bdr = mesh ? mesh->bdr_attributes.Max() : 0;

        for (size_t i = 0; i < boundaries.size(); ++i) {
            const auto& bc = boundaries[i];
            std::string prefix = "boundaries[" + std::to_string(i) + "]";

            if (!bc.contains("type")) {
                AddError(prefix + ".type", "Missing required field 'type'");
            } else {
                std::string bc_type = bc["type"];
                if (bc_type != "Dirichlet" && bc_type != "Neumann" && bc_type != "Robin") {
                    AddError(prefix + ".type", "Invalid boundary type '" + bc_type + "'. Must be 'Dirichlet', 'Neumann', or 'Robin'");
                }
            }

            if (!bc.contains("attribute_ids")) {
                AddError(prefix + ".attribute_ids", "Missing required field 'attribute_ids'");
            } else if (!bc["attribute_ids"].is_array()) {
                AddError(prefix + ".attribute_ids", "Attribute ids must be an array");
            } else if (mesh) {
                for (int attr : bc["attribute_ids"]) {
                    if (attr < 0 || attr > max_bdr) {
                        AddError(prefix + ".attribute_ids", "Boundary attribute " + std::to_string(attr) +
                                " is out of range [0, " + std::to_string(max_bdr) + "]");
                    }
                }
            }

            if (!bc.contains("value")) {
                AddError(prefix + ".value", "Missing required field 'value'");
            }
        }
    }

    void ValidateTerminals(const json& config, const mfem::Mesh* mesh = nullptr) {
        std::string type = config["simulation"].value("type", "");

        // Magnetic simulations need at least one current terminal to excite the field.
        if ((type == "magnetostatics" || type == "magnetoquasistatics") && !config.contains("terminals")) {
            AddError("terminals", "Magnetic simulations require at least one terminal");
            return;
        }

        if (!config.contains("terminals")) {
            return;
        }

        const auto& terminals = config["terminals"];
        if (!terminals.is_array()) {
            AddError("terminals", "Terminals must be an array");
            return;
        }

        const int max_dom = mesh ? mesh->attributes.Max() : 0;
        const int max_bdr = mesh ? mesh->bdr_attributes.Max() : 0;

        for (size_t i = 0; i < terminals.size(); ++i) {
            const auto& t = terminals[i];
            std::string prefix = "terminals[" + std::to_string(i) + "]";

            if (!t.contains("name")) {
                AddError(prefix + ".name", "Missing required field 'name'");
            }

            std::string drive = t.value("drive", "voltage");
            if (drive != "voltage" && drive != "current") {
                AddError(prefix + ".drive", "Invalid drive '" + drive + "'. Must be 'voltage' or 'current'");
            }

            // Voltage terminals bind to boundary attributes; current terminals to domain attributes.
            const bool is_current = (drive == "current");
            const int max_attr = is_current ? max_dom : max_bdr;

            if (!t.contains("attribute_ids")) {
                AddError(prefix + ".attribute_ids", "Missing required field 'attribute_ids'");
            } else if (!t["attribute_ids"].is_array()) {
                AddError(prefix + ".attribute_ids", "Attribute ids must be an array");
            } else if (mesh) {
                for (int attr : t["attribute_ids"]) {
                    if (attr <= 0 || attr > max_attr) {
                        AddError(prefix + ".attribute_ids", "Terminal attribute " + std::to_string(attr) +
                                " is out of range [1, " + std::to_string(max_attr) + "]");
                    }
                }
            }
        }
    }

    void ValidateScenarios(const json& config, const mfem::Mesh* mesh = nullptr) {
        if (!config.contains("scenarios")) {
            return;
        }

        const auto& scenarios = config["scenarios"];
        if (!scenarios.is_array()) {
            AddError("scenarios", "Scenarios must be an array");
            return;
        }

        // Build the set of declared terminal names and which are current-driven,
        // so drives can be cross-checked (catches typo'd terminal references).
        std::set<std::string> terminal_names;
        std::set<std::string> current_terminals;
        if (config.contains("terminals") && config["terminals"].is_array()) {
            for (const auto& t : config["terminals"]) {
                std::string name = t.value("name", "");
                if (!name.empty()) {
                    terminal_names.insert(name);
                    if (t.value("drive", "voltage") == "current") {
                        current_terminals.insert(name);
                    }
                }
            }
        }

        for (size_t i = 0; i < scenarios.size(); ++i) {
            const auto& sc = scenarios[i];
            std::string prefix = "scenarios[" + std::to_string(i) + "]";

            if (!sc.contains("drives")) {
                continue;
            }
            if (!sc["drives"].is_array()) {
                AddError(prefix + ".drives", "Drives must be an array");
                continue;
            }

            for (size_t j = 0; j < sc["drives"].size(); ++j) {
                const auto& d = sc["drives"][j];
                std::string dprefix = prefix + ".drives[" + std::to_string(j) + "]";

                if (!d.contains("terminal")) {
                    AddError(dprefix + ".terminal", "Missing required field 'terminal'");
                    continue;
                }

                std::string tname = d["terminal"];
                if (terminal_names.find(tname) == terminal_names.end()) {
                    AddError(dprefix + ".terminal", "Unknown terminal '" + tname +
                            "'. No terminal with that name is declared");
                }

                if (d.value("floating", false) && current_terminals.count(tname)) {
                    AddError(dprefix + ".floating", "Terminal '" + tname +
                            "' is current-driven; 'floating' applies to voltage terminals only");
                }
            }
        }
    }

public:
    /**
     * @brief Validate a configuration against a mesh
     * @param config The JSON configuration
     * @param mesh Optional mesh pointer for attribute range checking
     * @return true if configuration is valid
     */
    bool Validate(const json& config, const mfem::Mesh* mesh = nullptr) {
        errors.clear();

        ValidateSimulation(config);
        ValidateMaterials(config, mesh);
        ValidateRegions(config, mesh);
        ValidateBoundaries(config, mesh);
        ValidateTerminals(config, mesh);
        ValidateScenarios(config, mesh);

        return errors.empty();
    }

    /**
     * @brief Get all validation errors
     */
    [[nodiscard]] const std::vector<ValidationError>& GetErrors() const {
        return errors;
    }

    /**
     * @brief Get formatted error message
     */
    [[nodiscard]] std::string GetErrorMessage() const {
        if (errors.empty()) {
            return "No errors";
        }

        std::string msg = "Configuration validation failed:\n";
        for (const auto& err : errors) {
            msg += "  - " + err.field + ": " + err.message + "\n";
        }
        return msg;
    }

    /**
     * @brief Validate and throw if invalid
     * @throws std::runtime_error with detailed error message
     */
    void ValidateOrThrow(const json& config, const mfem::Mesh* mesh = nullptr) {
        if (!Validate(config, mesh)) {
            throw std::runtime_error(GetErrorMessage());
        }
    }
};
