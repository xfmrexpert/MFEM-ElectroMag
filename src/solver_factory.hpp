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
    using SolverCreator = std::function<std::unique_ptr<PhysicsSolver>(mfem::Mesh&, const json&)>;

private:
    std::unordered_map<PhysicsType, SolverCreator> registry;

    SolverFactory() {
        // Register all available solvers
        Register(PhysicsType::Electrostatics,
            [](mfem::Mesh& mesh, const json& config) -> std::unique_ptr<PhysicsSolver> {
                return std::make_unique<ElectrostaticSolver>(mesh, config);
            });
        
        Register(PhysicsType::Magnetostatics,
            [](mfem::Mesh& mesh, const json& config) -> std::unique_ptr<PhysicsSolver> {
                return std::make_unique<MagnetostaticSolver>(mesh, config);
            });

        Register(PhysicsType::Magnetoquasistatics,
            [](mfem::Mesh& mesh, const json& config) -> std::unique_ptr<PhysicsSolver> {
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
     * @param physics The physics formulation this creator handles
     * @param creator Factory function that creates the solver
     */
    void Register(PhysicsType physics, SolverCreator creator) {
        registry[physics] = creator;
    }

    /**
     * @brief Create a solver based on the physics formulation
     * @param physics The physics formulation (from ProblemConfig)
     * @param mesh The mesh object
     * @param config The JSON configuration
     * @return Unique pointer to the created solver
     * @throws std::runtime_error if the physics type is not registered
     */
    [[nodiscard]] std::unique_ptr<PhysicsSolver> Create(
        PhysicsType physics,
        mfem::Mesh& mesh,
        const json& config) const {

        auto it = registry.find(physics);
        if (it == registry.end()) {
            throw std::runtime_error(std::string("Unregistered physics: ") + ToString(physics));
        }
        return it->second(mesh, config);
    }

    /**
     * @brief Check if a solver type is registered
     * @param physics The physics formulation
     * @return true if the type is registered
     */
    [[nodiscard]] bool IsRegistered(PhysicsType physics) const {
        return registry.find(physics) != registry.end();
    }

    // Delete copy constructor and assignment operator
    SolverFactory(const SolverFactory&) = delete;
    SolverFactory& operator=(const SolverFactory&) = delete;
};
