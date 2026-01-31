// Copyright (c) 2026 T. C. Raymond
// SPDX-License-Identifier: MIT

#pragma once

#include "mfem.hpp"
#include "constants.hpp"

/**
 * @brief Axisymmetric Mass Integrator for Eddy Currents
 * Solves: Integral( sigma * A * v * 2*pi*r * dr * dz )
 * Used for the mass term in time-harmonic problems (j * omega * sigma * A)
 */
class AxisymmetricMassIntegrator : public mfem::BilinearFormIntegrator
{
private:
   mfem::Coefficient *Q; // Conductivity (sigma)
   static constexpr double factor = Constants::TWO_PI;

public:
   AxisymmetricMassIntegrator(mfem::Coefficient &q) : Q(&q) 
   {
      MFEM_ASSERT(Q != nullptr, "Coefficient cannot be null");
   }

   void AssembleElementMatrix(const mfem::FiniteElement &el,
                               mfem::ElementTransformation &Trans,
                               mfem::DenseMatrix &elmat) override
   {
      int nd = el.GetDof();
      elmat.SetSize(nd);
      elmat = 0.0;

      mfem::Vector shape(nd);
      mfem::Vector pos(2);

      // Standard order is fine here (no 1/r singularity)
      // r dependency is linear, so Order + 1 is sufficient
      int order = 2 * el.GetOrder() + 1;
      const mfem::IntegrationRule *ir = &mfem::IntRules.Get(el.GetGeomType(), order);

      for (int i = 0; i < ir->GetNPoints(); i++)
      {
         const mfem::IntegrationPoint &ip = ir->IntPoint(i);
         Trans.SetIntPoint(&ip);
         Trans.Transform(ip, pos);
         
         double r = pos(0);

         double val = Q->Eval(Trans, ip); // Conductivity sigma
         
         // Weight = w * det(J) * (2*pi*r) * sigma
         double w = ip.weight * Trans.Weight() * (factor * r) * val;

         el.CalcShape(ip, shape);

         for (int j = 0; j < nd; j++)
         {
             for (int k = 0; k <= j; k++) // Symmetry
             {
                 double entry = w * shape(j) * shape(k);
                 elmat(j, k) += entry;
                 if (j != k) elmat(k, j) += entry;
             }
         }
      }
   }
};