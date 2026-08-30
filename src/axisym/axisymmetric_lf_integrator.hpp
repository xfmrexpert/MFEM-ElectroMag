// Copyright (c) 2026 T. C. Raymond
// SPDX-License-Identifier: MIT

#pragma once

#include "mfem.hpp"
#include "axisymmetric_measure.hpp"

// -----------------------------------------------------------------------------
// 2. Linear Form Integrator: Current Density (J_phi)
// -----------------------------------------------------------------------------
// Solves: Integral( J_phi * v * 2*pi*r * dr * dz )
// Carries the full axisymmetric measure (see axisymmetric_measure.hpp).
class AxisymmetricLFIntegrator : public mfem::LinearFormIntegrator
{
private:
   mfem::Coefficient *src; // e.g., J_phi (magnetostatics) or rho (electrostatics)

public:
   explicit AxisymmetricLFIntegrator(
      mfem::Coefficient &q, const mfem::IntegrationRule *ir = nullptr)
      : mfem::LinearFormIntegrator(ir), src(&q)
   {
      MFEM_ASSERT(src != nullptr, "Coefficient cannot be null");
   }

   static const mfem::IntegrationRule &GetRule(
      const mfem::FiniteElement &el,
      const mfem::ElementTransformation &Trans)
   {
      // The source is piecewise constant; the remaining factors are the test
      // function, radial coordinate, and Jacobian determinant.
      const int order = el.GetOrder() + Trans.Order() + Trans.OrderW();

      if (el.Space() == mfem::FunctionSpace::rQk)
      {
         return mfem::RefinedIntRules.Get(el.GetGeomType(), order);
      }
      return mfem::IntRules.Get(el.GetGeomType(), order);
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

      const mfem::IntegrationRule *ir = GetIntegrationRule(el, Trans);

      for (int i = 0; i < ir->GetNPoints(); i++)
      {
         const mfem::IntegrationPoint &ip = ir->IntPoint(i);
         Trans.SetIntPoint(&ip);
         Trans.Transform(ip, pos);

         const double r = pos(0);
         const double val = src->Eval(Trans, ip);

         const double w = ip.weight * Trans.Weight() * Axisymmetric::Measure(r) * val;

         el.CalcShape(ip, shape);
         elvect.Add(w, shape);
      }
   }

protected:
   const mfem::IntegrationRule *GetDefaultIntegrationRule(
      const mfem::FiniteElement &trial_fe,
      const mfem::FiniteElement &,
      const mfem::ElementTransformation &Trans) const override
   {
      return &GetRule(trial_fe, Trans);
   }
};
