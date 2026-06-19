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

    // Entity group names declared in the "entity_groups" section, split by the
    // dimensionality the InputParser assigns (dim == 1 => boundary, otherwise
    // domain). Populated by ValidateEntityGroups and consulted when checking
    // that regions/terminals/boundaries reference a group of the right kind.
    std::set<std::string> entity_group_names_;
    std::set<std::string> boundary_group_names_;
    std::set<std::string> domain_group_names_;

    void AddError(const std::string& field, const std::string& message) {
        errors.emplace_back(field, message);
    }

    // "simulation.physics_type" as a string, or "" when absent.
    [[nodiscard]] std::string PhysicsType(const json& config) const {
        if (config.contains("simulation") && config["simulation"].is_object()) {
            return config["simulation"].value("physics_type", std::string{});
        }
        return {};
    }

    // Verify that an "entity_group" reference points at a declared group and,
    // when expected_kind is set, that the group has the matching dimensionality
    // (1 => boundary, 2 => domain; 0 => either). Mirrors how the solvers resolve
    // EntityGroupName -> EntityGroup via config.EntityGroups.at(name).
    void CheckEntityGroupRef(const json& node, const std::string& field, int expected_kind) {
        if (!node.is_string()) {
            AddError(field, "Entity group reference must be a string");
            return;
        }
        const std::string ref = node.get<std::string>();
        if (entity_group_names_.find(ref) == entity_group_names_.end()) {
            AddError(field, "Unknown entity group '" + ref +
                    "'. No entity group with that name is declared in 'entity_groups'");
            return;
        }
        if (expected_kind == 1 && boundary_group_names_.count(ref) == 0) {
            AddError(field, "Entity group '" + ref + "' must be a boundary group (dim = 1)");
        } else if (expected_kind == 2 && domain_group_names_.count(ref) == 0) {
            AddError(field, "Entity group '" + ref + "' must be a domain group (dim != 1)");
        }
    }

    void ValidateSimulation(const json& config) {
        if (!config.contains("simulation")) {
            AddError("simulation", "Missing required 'simulation' section");
            return;
        }

        const auto& sim = config["simulation"];

        // Required fields
        if (!sim.contains("physics_type")) {
            AddError("simulation.physics_type", "Missing required field 'physics_type'");
        } else {
            std::string type = sim["physics_type"];
            if (type != "electrostatics" && type != "magnetostatics" && type != "magnetoquasistatics") {
                AddError("simulation.physics_type", "Invalid physics_type '" + type + "'. Must be 'electrostatics', 'magnetostatics', or 'magnetoquasistatics'");
            }

            // Physics-specific requirements
            if (type == "magnetoquasistatics" && !sim.contains("frequency")) {
                AddError("simulation.frequency", "Magnetoquasistatic simulations require 'frequency' field");
            }
        }

        if (!sim.contains("mesh")) {
            AddError("simulation.mesh", "Missing required field 'mesh'");
        }

        // Optional enumerated fields
        if (sim.contains("geometry_type")) {
            std::string g = sim["geometry_type"];
            if (g != "axisymmetric" && g != "planar") {
                AddError("simulation.geometry_type", "Invalid geometry_type '" + g + "'. Must be 'axisymmetric' or 'planar'");
            }
        }

        if (sim.contains("analysis_type")) {
            std::string a = sim["analysis_type"];
            if (a != "field" && a != "coupling_matrix") {
                AddError("simulation.analysis_type", "Invalid analysis_type '" + a + "'. Must be 'field' or 'coupling_matrix'");
            }
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

    void ValidateEntityGroups(const json& config, const mfem::Mesh* mesh = nullptr) {
        entity_group_names_.clear();
        boundary_group_names_.clear();
        domain_group_names_.clear();

        if (!config.contains("entity_groups")) {
            return;
        }

        const auto& groups = config["entity_groups"];
        if (!groups.is_array()) {
            AddError("entity_groups", "Entity groups must be an array");
            return;
        }

        const int max_dom = mesh ? mesh->attributes.Max() : 0;
        const int max_bdr = mesh ? mesh->bdr_attributes.Max() : 0;

        for (size_t i = 0; i < groups.size(); ++i) {
            const auto& g = groups[i];
            std::string prefix = "entity_groups[" + std::to_string(i) + "]";

            std::string name = g.value("name", std::string{});
            if (name.empty()) {
                AddError(prefix + ".name", "Missing required field 'name'");
            } else if (!entity_group_names_.insert(name).second) {
                AddError(prefix + ".name", "Duplicate entity group name '" + name + "'");
            }

            // dim == 1 => boundary group, otherwise domain group (matches InputParser).
            bool is_boundary = false;
            if (!g.contains("dim")) {
                AddError(prefix + ".dim", "Missing required field 'dim'");
            } else {
                int dim = g["dim"];
                is_boundary = (dim == 1);
            }

            if (!name.empty()) {
                if (is_boundary) boundary_group_names_.insert(name);
                else             domain_group_names_.insert(name);
            }

            if (!g.contains("attribute_ids")) {
                AddError(prefix + ".attribute_ids", "Missing required field 'attribute_ids'");
            } else if (!g["attribute_ids"].is_array()) {
                AddError(prefix + ".attribute_ids", "Attribute ids must be an array");
            } else if (mesh) {
                const int max_attr = is_boundary ? max_bdr : max_dom;
                for (int attr : g["attribute_ids"]) {
                    if (attr <= 0 || attr > max_attr) {
                        AddError(prefix + ".attribute_ids", "Attribute " + std::to_string(attr) +
                                " is out of range [1, " + std::to_string(max_attr) + "]");
                    }
                }
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

        size_t num_materials = 0;
        if (config.contains("materials") && config["materials"].is_array()) {
            num_materials = config["materials"].size();
        }

        for (size_t i = 0; i < regions.size(); ++i) {
            const auto& reg = regions[i];
            std::string prefix = "regions[" + std::to_string(i) + "]";

            if (!reg.contains("entity_group")) {
                AddError(prefix + ".entity_group", "Missing required field 'entity_group'");
            } else {
                CheckEntityGroupRef(reg["entity_group"], prefix + ".entity_group", /*domain*/2);
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

        std::string type = PhysicsType(config);

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

            if (!bc.contains("entity_group")) {
                AddError(prefix + ".entity_group", "Missing required field 'entity_group'");
            } else {
                CheckEntityGroupRef(bc["entity_group"], prefix + ".entity_group", /*boundary*/1);
            }

            if (!bc.contains("value")) {
                AddError(prefix + ".value", "Missing required field 'value'");
            }
        }
    }

    void ValidateTerminals(const json& config, const mfem::Mesh* mesh = nullptr) {
        std::string type = PhysicsType(config);

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

        for (size_t i = 0; i < terminals.size(); ++i) {
            const auto& t = terminals[i];
            std::string prefix = "terminals[" + std::to_string(i) + "]";

            if (!t.contains("name")) {
                AddError(prefix + ".name", "Missing required field 'name'");
            }

            std::string excitation = t.value("excitation", "voltage");
            if (excitation != "voltage" && excitation != "current") {
                AddError(prefix + ".excitation", "Invalid excitation '" + excitation + "'. Must be 'voltage' or 'current'");
            }

            // Voltage terminals bind to boundary groups; current terminals to domain groups.
            const bool is_current = (excitation == "current");
            if (!t.contains("entity_group")) {
                AddError(prefix + ".entity_group", "Missing required field 'entity_group'");
            } else {
                CheckEntityGroupRef(t["entity_group"], prefix + ".entity_group",
                                    is_current ? /*domain*/2 : /*boundary*/1);
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
        // so excitations can be cross-checked (catches typo'd terminal references).
        std::set<std::string> terminal_names;
        std::set<std::string> current_terminals;
        if (config.contains("terminals") && config["terminals"].is_array()) {
            for (const auto& t : config["terminals"]) {
                std::string name = t.value("name", "");
                if (!name.empty()) {
                    terminal_names.insert(name);
                    if (t.value("excitation", "voltage") == "current") {
                        current_terminals.insert(name);
                    }
                }
            }
        }

        for (size_t i = 0; i < scenarios.size(); ++i) {
            const auto& sc = scenarios[i];
            std::string prefix = "scenarios[" + std::to_string(i) + "]";

            if (!sc.contains("excitations")) {
                continue;
            }
            if (!sc["excitations"].is_array()) {
                AddError(prefix + ".excitations", "Excitations must be an array");
                continue;
            }

            for (size_t j = 0; j < sc["excitations"].size(); ++j) {
                const auto& d = sc["excitations"][j];
                std::string dprefix = prefix + ".excitations[" + std::to_string(j) + "]";

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
        ValidateEntityGroups(config, mesh);
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
