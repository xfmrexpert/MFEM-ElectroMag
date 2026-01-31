// Copyright (c) 2026 T. C. Raymond
// SPDX-License-Identifier: MIT

#include <iostream>
#include <memory>
#include <string>
#include "input_parser.hpp"
#include "physics_solver.hpp"
#include "solver_factory.hpp"
#include "config_validator.hpp"

int main(int argc, char *argv[]) {
    try {
        // 1. Shared Infrastructure
        std::string config_file = (argc > 1) ? argv[1] : "config.json";
        InputParser parser(config_file);

        // 2. Validate Configuration (basic validation before mesh loading)
        ConfigValidator validator;
        validator.ValidateOrThrow(parser.config);

        // 3. Load Mesh
        mfem::Mesh mesh(parser.GetMeshPath(), 1, 1);

        // 4. Validate Configuration (with mesh for attribute checking)
        validator.ValidateOrThrow(parser.config, &mesh);

        // 5. Factory Logic - Create Solver
        std::string type = parser.config["simulation"]["type"];
        auto solver = SolverFactory::Instance().Create(type, mesh, parser.config);

        // 6. Execution
        solver->Setup();
        solver->Run();
        solver->Save();

        return 0;
    }
    catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    catch (...) {
        std::cerr << "Error: Unknown exception occurred" << std::endl;
        return 2;
    }
}