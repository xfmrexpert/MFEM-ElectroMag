// Copyright (c) 2026 T. C. Raymond
// SPDX-License-Identifier: MIT

#pragma once

#include "mfem.hpp"

/**
 * @brief Computes B = Curl(A) in 2D Axisymmetry
 * B_r = -dA/dz
 * B_z = A/r + dA/dr, with the exact axis limit B_z -> 2 dA/dr at r = 0.
 */
class MagneticFieldCoefficient : public mfem::VectorCoefficient
{
private:
   mfem::GridFunction *A; // The solution (Magnetic Vector Potential)

public:
   explicit MagneticFieldCoefficient(mfem::GridFunction* a_gf)
      : mfem::VectorCoefficient(2), A(a_gf) { }

   void Eval(mfem::Vector &B, mfem::ElementTransformation &T, 
             const mfem::IntegrationPoint &ip) override
   {
      T.SetIntPoint(&ip);

      mfem::Vector pos(2); // 2D axisymmetric
      T.Transform(ip, pos);
      const double r = pos(0);

      // Get Value (A) and Gradient (dA/dr, dA/dz)
      const double A_val = A->GetValue(T, ip);

      mfem::Vector grad_A(2);
      A->GetGradient(T, grad_A);

      B.SetSize(2);
      // Regularity forces A_phi -> 0 linearly as r -> 0, so l'Hopital gives
      // the finite limit B_z -> 2 dA/dr on the axis.
      B(0) = -grad_A(1);
      B(1) = (r == 0.0) ? (2.0 * grad_A(0))
                        : (grad_A(0) + A_val / r);
   }
};

class PlanarMagneticFieldCoefficient : public mfem::VectorCoefficient
{
private:
   mfem::GridFunction* A;

public:
   explicit PlanarMagneticFieldCoefficient(mfem::GridFunction* a_gf)
      : mfem::VectorCoefficient(2), A(a_gf) { }

   void Eval(mfem::Vector& B, mfem::ElementTransformation& T,
             const mfem::IntegrationPoint& ip) override
   {
      T.SetIntPoint(&ip);
      mfem::Vector grad_A(2);
      A->GetGradient(T, grad_A);

      B.SetSize(2);
      B(0) = grad_A(1);
      B(1) = -grad_A(0);
   }
};
