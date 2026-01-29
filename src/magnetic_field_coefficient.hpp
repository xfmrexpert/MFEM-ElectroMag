// Copyright (c) 2026 T. C. Raymond
// SPDX-License-Identifier: MIT

#pragma once

#include "mfem.hpp"

/**
 * @brief Computes B = Curl(A) in 2D Axisymmetry
 * B_r = -dA/dz
 * B_z = A/r + dA/dr
 */
class MagneticFieldCoefficient : public mfem::VectorCoefficient
{
private:
   mfem::GridFunction *A; // The solution (Magnetic Vector Potential)
   static constexpr double r_tol = 1e-12;

public:
   MagneticFieldCoefficient(mfem::GridFunction &a_gf) 
      : mfem::VectorCoefficient(2), A(&a_gf) { }

   void Eval(mfem::Vector &B, mfem::ElementTransformation &T, 
             const mfem::IntegrationPoint &ip) override
   {
      mfem::Vector pos(2); // 2D axisymmetric
      T.Transform(ip, pos);
      double r = pos(0);

      // Get Value (A) and Gradient (dA/dr, dA/dz)
      double A_val = A->GetValue(T, ip);
      
      mfem::Vector grad_A(2);
      A->GetGradient(T, grad_A);
      double dA_dr = grad_A(0);
      double dA_dz = grad_A(1);

      B.SetSize(2);

      // B_r = -dA/dz
      B(0) = -dA_dz;

      // B_z = A/r + dA/dr
      if (r > r_tol) {
         B(1) = (A_val / r) + dA_dr;
      } else {
         // Limit r->0
         B(1) = 2.0 * dA_dr;
      }
   }
};