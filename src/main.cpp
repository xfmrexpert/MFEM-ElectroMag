#include <unordered_set>
#include <memory>
#include "input_parser.hpp"
#include "physics_solver.hpp"
#include "electrostatic_solver.hpp"
#include "magnetostatic_solver.hpp"
#include "magnetoquasistatic_solver.hpp"

using namespace std;
using namespace mfem;

int main(int argc, char *argv[]) {
    // 1. Shared Infrastructure
    std::string config_file = (argc > 1) ? argv[1] : "config.json";
    InputParser parser(config_file);
    Mesh mesh(parser.GetMeshPath(), 1, 1);

    // 2. Factory Logic
    std::unique_ptr<PhysicsSolver> solver;
    std::string type = parser.config["simulation"]["type"];

    if (type == "electrostatics") {
        solver = std::make_unique<ElectrostaticSolver>(mesh, parser.config);
    } 
    else if (type == "magnetostatics") {
        solver = std::make_unique<MagnetostaticSolver>(mesh, parser.config);
    }
    else if (type == "magnetoquasistatics") {
        solver = std::make_unique<MagnetoquasistaticSolver>(mesh, parser.config);
    }
    else {
        std::cerr << "Unknown physics type: " << type << std::endl;
        return 1;
    }

    // 3. Execution
    solver->Setup();
    solver->Run();
    solver->Save();

    return 0;
}