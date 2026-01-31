// Copyright (c) 2026 T. C. Raymond
// SPDX-License-Identifier: MIT

#pragma once

#include <memory>
#include <unordered_map>
#include <functional>
#include <string>
#include <vector>
#include <stdexcept>
#include "mfem.hpp"
#include "json.hpp"
#include "physics_solver.hpp"
#include "electrostatic_solver.hpp"
#include "magnetostatic_solver.hpp"
#include "magnetoquasistatic_solver.hpp"

using json = nlohmann::json;

/**
 * @brief Factory for creating physics solvers based on simulation type
 */
class SolverFactory {
public:
    using SolverCreator = std::function<std::unique_ptr<PhysicsSolver>(mfem::Mesh&, json&)>;

private:
    std::unordered_map<std::string, SolverCreator> registry;

    SolverFactory() {
        // Register all available solvers
        Register("electrostatics",
            [](mfem::Mesh& mesh, json& config) -> std::unique_ptr<PhysicsSolver> {
                return std::make_unique<ElectrostaticSolver>(mesh, config);
            });

        Register("magnetostatics",
            [](mfem::Mesh& mesh, json& config) -> std::unique_ptr<PhysicsSolver> {
                return std::make_unique<MagnetostaticSolver>(mesh, config);
            });

        Register("magnetoquasistatics",
            [](mfem::Mesh& mesh, json& config) -> std::unique_ptr<PhysicsSolver> {
                return std::make_unique<MagnetoquasistaticSolver>(mesh, config);
            });
    }

public:
    /**
     * @brief Get the singleton instance of the factory
     */
    static SolverFactory& Instance() {
        static SolverFactory instance;
        return instance;
    }

    /**
     * @brief Register a new solver type
     * @param type The string identifier for the physics type
     * @param creator Factory function that creates the solver
     */
    void Register(const std::string& type, SolverCreator creator) {
        registry[type] = creator;
    }

    /**
     * @brief Create a solver based on the simulation type
     * @param type The physics type string (from config)
     * @param mesh The mesh object
     * @param config The JSON configuration
     * @return Unique pointer to the created solver
     * @throws std::runtime_error if the solver type is not registered
     */
    [[nodiscard]] std::unique_ptr<PhysicsSolver> Create(
        const std::string& type,
        mfem::Mesh& mesh,
        json& config) const {

        auto it = registry.find(type);
        if (it == registry.end()) {
            throw std::runtime_error("Unknown physics type: " + type);
        }
        return it->second(mesh, config);
    }

    /**
     * @brief Check if a solver type is registered
     * @param type The physics type string
     * @return true if the type is registered
     */
    [[nodiscard]] bool IsRegistered(const std::string& type) const {
        return registry.find(type) != registry.end();
    }

    /**
     * @brief Get list of all registered solver types
     * @return Vector of registered type names
     */
    [[nodiscard]] std::vector<std::string> GetRegisteredTypes() const {
        std::vector<std::string> types;
        types.reserve(registry.size());
        for (const auto& [type, _] : registry) {
            types.push_back(type);
        }
        return types;
    }

    // Delete copy constructor and assignment operator
    SolverFactory(const SolverFactory&) = delete;
    SolverFactory& operator=(const SolverFactory&) = delete;
};
