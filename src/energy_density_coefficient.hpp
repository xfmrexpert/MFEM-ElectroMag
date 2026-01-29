// Copyright (c) 2026 T. C. Raymond
// SPDX-License-Identifier: MIT

#pragma once

#include "mfem.hpp"

// -----------------------------------------------------------------------------
// 3. Energy Density Coefficient
// -----------------------------------------------------------------------------
// Computes local energy density: w = 0.5 * eps * |E|^2
// Useful for visualization. Note: This does NOT include 2*pi*r. 
// It represents the physical density at that point in space.
class EnergyDensityCoefficient : public mfem::Coefficient
{
private:
   mfem::GridFunction *Phi;
   mfem::Coefficient *Eps;

public:
   EnergyDensityCoefficient(mfem::GridFunction &phi, mfem::Coefficient &eps) 
      : Phi(&phi), Eps(&eps) { }

   double Eval(mfem::ElementTransformation &T, 
               const mfem::IntegrationPoint &ip) override
   {
      mfem::Vector grad(2);
      Phi->GetGradient(T, grad); // Returns (Er, Ez)
      
      double mag_sq = grad * grad; // |E|^2
      double eps_val = Eps->Eval(T, ip);

      return 0.5 * eps_val * mag_sq;
   }
};
