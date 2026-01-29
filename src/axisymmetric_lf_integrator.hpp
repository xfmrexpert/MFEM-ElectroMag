#pragma once

#include "mfem.hpp"

// -----------------------------------------------------------------------------
// 2. Linear Form Integrator: Space Charge (rho)
// -----------------------------------------------------------------------------
// Solves: Integral( rho * v * 2*pi*r * dr * dz )
class AxisymmetricLFIntegrator : public mfem::LinearFormIntegrator
{
private:
   mfem::Coefficient *Q; // Represents Space Charge Density (rho)
   static constexpr double factor = 2.0 * M_PI;

public:
   AxisymmetricLFIntegrator(mfem::Coefficient &q) : Q(&q) 
   {
      MFEM_ASSERT(Q != nullptr, "Coefficient cannot be null");
   }

   void AssembleRHSElementVect(const mfem::FiniteElement &el,
                                mfem::ElementTransformation &Trans,
                                mfem::Vector &elvect) override
   {
      int nd = el.GetDof();
      elvect.SetSize(nd);
      elvect = 0.0;

      mfem::Vector shape(nd);
      mfem::Vector pos(2); // 2D axisymmetric

      // Increase order by 1 for 'r' term
      int order = 2 * el.GetOrder() + 1; 
      const mfem::IntegrationRule *ir = &mfem::IntRules.Get(el.GetGeomType(), order);

      for (int i = 0; i < ir->GetNPoints(); i++)
      {
         const mfem::IntegrationPoint &ip = ir->IntPoint(i);
         Trans.SetIntPoint(&ip);
         Trans.Transform(ip, pos);

         double r = pos(0);

         // Weight = quad_weight * det(J) * (2 * pi * r) * rho
         double val = Q->Eval(Trans, ip);
         double w = ip.weight * Trans.Weight() * (factor * r) * val;

         el.CalcShape(ip, shape);
         add(elvect, w, shape, elvect);
      }
   }
};