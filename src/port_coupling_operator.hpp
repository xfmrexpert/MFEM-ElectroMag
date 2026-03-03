#pragma once

#include "mfem.hpp"

class PortCouplingOperator : public mfem::Operator
{
private:
    int num_ports;
    mfem::Array<mfem::Vector*> port_forms;

public:
    PortCouplingOperator(int n_dofs, int n_ports, const mfem::Array<mfem::Vector*> forms)
        : Operator(n_dofs, n_ports), num_ports(n_ports), port_forms(forms)
    { }

    // Computes y = S_Av * x
    // Maps the vector of port voltages 'x' to the FE space vector 'y'
    virtual void Mult(const mfem::Vector &x, mfem::Vector &y) const override
    {
        MFEM_ASSERT(x.Size() == num_ports, "Input vector size must match number of ports.");
        MFEM_ASSERT(y.Size() == height, "Output vector size mismatch.");

        y = 0.0;
        for (int i = 0; i < num_ports; ++i) {
            // y += x_i * port_form_i
            y.Add(x(i), *port_forms[i]); 
        }
    }

    // Computes y = S_vA * x 
    // Maps the FE space vector 'x' to the vector of induced currents 'y'
    virtual void MultTranspose(const mfem::Vector &x, mfem::Vector &y) const override
    {
        MFEM_ASSERT(x.Size() == height, "Input vector size mismatch.");
        MFEM_ASSERT(y.Size() == num_ports, "Output vector size must match number of ports.");

        for (int i = 0; i < num_ports; ++i) {
            // y_i = inner product of port_form_i and x
            y(i) = (*port_forms[i]) * x; 
        }
    }

};