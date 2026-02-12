#include "mfem.hpp"
#include "../src/axisymmetric_mass_integrator.hpp"
#include <chrono>
#include <iostream>

using namespace mfem;
using namespace std;

int main() {
    // 1. Setup Mesh: Simple 10x10 quad mesh
    // Note: older MFEM constructors might be deprecated, using a simpler one or ignoring warning
    Mesh mesh(10, 10, Element::QUADRILATERAL, 1, 1.0, 1.0);

    // 2. Setup Finite Element Space
    H1_FECollection fec(1, mesh.Dimension());
    FiniteElementSpace fes(&mesh, &fec);

    // 3. Setup Integrator
    ConstantCoefficient one(1.0);
    AxisymmetricMassIntegrator integrator(one);

    // 4. Prepare for loop
    DenseMatrix elmat;
    ElementTransformation *trans;
    const FiniteElement *el;

    // Use a large number of repetitions to get measurable time
    int num_repeats = 10000;
    int num_elements = fes.GetNE();

    auto start = chrono::high_resolution_clock::now();

    for (int r = 0; r < num_repeats; r++) {
        for (int i = 0; i < num_elements; i++) {
            trans = fes.GetElementTransformation(i);
            el = fes.GetFE(i);
            integrator.AssembleElementMatrix(*el, *trans, elmat);
        }
    }

    auto end = chrono::high_resolution_clock::now();
    chrono::duration<double> diff = end - start;

    cout << "Total time: " << diff.count() << " s" << endl;
    double avg_ns = (diff.count() / (double)(num_repeats * num_elements)) * 1e9;
    cout << "Average time per AssembleElementMatrix: " << avg_ns << " ns" << endl;

    return 0;
}
