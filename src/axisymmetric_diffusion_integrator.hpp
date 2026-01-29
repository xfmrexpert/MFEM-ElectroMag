// Copyright (c) 2026 T. C. Raymond
// SPDX-License-Identifier: MIT

#pragma once

#include "mfem.hpp"

// -----------------------------------------------------------------------------
// 1. Stiffness Integrator: -div( eps * grad(u) )
// -----------------------------------------------------------------------------
// Solves: Integral( eps * grad(u) . grad(v) * 2*pi*r * dr * dz )
class AxisymmetricDiffusionIntegrator : public mfem::BilinearFormIntegrator
{
private:
   mfem::Coefficient *Q;
   static constexpr double factor = 2.0 * M_PI;

public:
   AxisymmetricDiffusionIntegrator(mfem::Coefficient &q) : Q(&q) 
   {
      MFEM_ASSERT(Q != nullptr, "Coefficient cannot be null");
   }

   void AssembleElementMatrix(const mfem::FiniteElement &el,
                               mfem::ElementTransformation &Trans,
                               mfem::DenseMatrix &elmat) override
   {
      int nd = el.GetDof();
      int dim = el.GetDim();
      double w;

      elmat.SetSize(nd);
      elmat = 0.0;

      mfem::DenseMatrix dshape(nd, dim);
      mfem::DenseMatrix dshapedxt(nd, dim);
      mfem::Vector pos(2); // 2D axisymmetric

      // Increase order by 1 for the linear 'r' term in the measure
      int order = 2 * el.GetOrder() + Trans.OrderGrad(&el) + 1;
      const mfem::IntegrationRule *ir = &mfem::IntRules.Get(el.GetGeomType(), order);

      for (int i = 0; i < ir->GetNPoints(); i++)
      {
         const mfem::IntegrationPoint &ip = ir->IntPoint(i);
         Trans.SetIntPoint(&ip);
         Trans.Transform(ip, pos);
         
         double r = pos(0); // Radius is X
         
         // Weight = quad_weight * det(J) * (2 * pi * r) * epsilon
         w = ip.weight * Trans.Weight() * (factor * r) * Q->Eval(Trans, ip);

         el.CalcDShape(ip, dshape);
         Mult(dshape, Trans.InverseJacobian(), dshapedxt);
         AddMult_a_AAt(w, dshapedxt, elmat);
      }
   }
};