// Copyright (c) 2026 T. C. Raymond
// SPDX-License-Identifier: MIT

#pragma once

#include <string>
#include <vector>
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
        int max_attr = mesh ? mesh->attributes.Max() : 0;

        for (size_t i = 0; i < materials.size(); ++i) {
            const auto& mat = materials[i];
            std::string prefix = "materials[" + std::to_string(i) + "]";

            if (!mat.contains("attributes")) {
                AddError(prefix + ".attributes", "Missing required field 'attributes'");
            } else if (!mat["attributes"].is_array()) {
                AddError(prefix + ".attributes", "Attributes must be an array");
            } else if (mesh) {
                // Validate attribute ranges
                for (int attr : mat["attributes"]) {
                    if (attr <= 0 || attr > max_attr) {
                        AddError(prefix + ".attributes", "Attribute " + std::to_string(attr) +
                                " is out of range [1, " + std::to_string(max_attr) + "]");
                    }
                }
            }

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

            if (!bc.contains("attributes")) {
                AddError(prefix + ".attributes", "Missing required field 'attributes'");
            } else if (!bc["attributes"].is_array()) {
                AddError(prefix + ".attributes", "Attributes must be an array");
            } else if (mesh) {
                for (int attr : bc["attributes"]) {
                    if (attr <= 0 || attr > max_bdr) {
                        AddError(prefix + ".attributes", "Boundary attribute " + std::to_string(attr) +
                                " is out of range [1, " + std::to_string(max_bdr) + "]");
                    }
                }
            }

            if (!bc.contains("value")) {
                AddError(prefix + ".value", "Missing required field 'value'");
            }
        }
    }

    void ValidateSources(const json& config, const mfem::Mesh* mesh = nullptr) {
        std::string type = config["simulation"].value("type", "");

        // Sources are required for magnetostatic/magnetoquasistatic
        if ((type == "magnetostatics" || type == "magnetoquasistatics") && !config.contains("sources")) {
            AddError("sources", "Magnetic simulations require at least one source");
            return;
        }

        if (!config.contains("sources")) {
            return;
        }

        const auto& sources = config["sources"];
        if (!sources.is_array()) {
            AddError("sources", "Sources must be an array");
            return;
        }

        int max_attr = mesh ? mesh->attributes.Max() : 0;

        for (size_t i = 0; i < sources.size(); ++i) {
            const auto& src = sources[i];
            std::string prefix = "sources[" + std::to_string(i) + "]";

            if (!src.contains("type")) {
                AddError(prefix + ".type", "Missing required field 'type'");
            }

            if (!src.contains("attributes")) {
                AddError(prefix + ".attributes", "Missing required field 'attributes'");
            } else if (!src["attributes"].is_array()) {
                AddError(prefix + ".attributes", "Attributes must be an array");
            } else if (mesh) {
                for (int attr : src["attributes"]) {
                    if (attr <= 0 || attr > max_attr) {
                        AddError(prefix + ".attributes", "Source attribute " + std::to_string(attr) +
                                " is out of range [1, " + std::to_string(max_attr) + "]");
                    }
                }
            }

            if (!src.contains("value")) {
                AddError(prefix + ".value", "Missing required field 'value'");
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
        ValidateBoundaries(config, mesh);
        ValidateSources(config, mesh);

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
