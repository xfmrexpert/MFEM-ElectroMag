// Copyright (c) 2026 T. C. Raymond
// SPDX-License-Identifier: MIT

#pragma once

#include "mfem.hpp"
#include "constants.hpp"

// -----------------------------------------------------------------------------
// 2. Linear Form Integrator: Current Density (J_phi)
// -----------------------------------------------------------------------------
// Solves: Integral( J_phi * v * r * dr * dz )
class AxisymmetricLFIntegrator : public mfem::LinearFormIntegrator
{
private:
   mfem::Coefficient *src; // e.g., J_phi (magnetostatics) or rho (electrostatics)

public:
   AxisymmetricLFIntegrator(mfem::Coefficient &q) : src(&q)
   {
      MFEM_ASSERT(src != nullptr, "Coefficient cannot be null");
   }

   void AssembleRHSElementVect(const mfem::FiniteElement &el,
                               mfem::ElementTransformation &Trans,
                               mfem::Vector &elvect) override
   {
      const int nd = el.GetDof();
      elvect.SetSize(nd);
      elvect = 0.0;

      mfem::Vector shape(nd);
      mfem::Vector pos(Trans.GetSpaceDim());

      const int order = 2 * el.GetOrder() + 1;
      const mfem::IntegrationRule &ir =
         mfem::IntRules.Get(el.GetGeomType(), order);

      for (int i = 0; i < ir.GetNPoints(); i++)
      {
         const mfem::IntegrationPoint &ip = ir.IntPoint(i);
         Trans.SetIntPoint(&ip);
         Trans.Transform(ip, pos);

         const double r = pos(0);
         const double val = src->Eval(Trans, ip);

         const double w = ip.weight * Trans.Weight() * r * val;

         el.CalcShape(ip, shape);
         elvect.Add(w, shape);
      }
   }
};
