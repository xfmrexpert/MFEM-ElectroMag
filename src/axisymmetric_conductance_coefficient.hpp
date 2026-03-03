#pragma once

#include "mfem.hpp"

// Custom coefficient for axisymmetric G_dc = sigma / (2 * pi * r)
class AxisymmetricConductanceCoeff : public mfem::Coefficient
{
private:
    double sigma;
public:
    AxisymmetricConductanceCoeff(double s) : sigma(s) { }

    virtual double Eval(mfem::ElementTransformation &T, const mfem::IntegrationPoint &ip) override
    {
        mfem::Vector transip;
        T.Transform(ip, transip);
        double r = transip(0); // Assuming 'r' is mapped to the first coordinate (x)
        
        // Prevent division by zero if the domain touches the axis of symmetry
        if (r < 1e-12) return 0.0; 
        
        return sigma / (2.0 * M_PI * r);
    }
};