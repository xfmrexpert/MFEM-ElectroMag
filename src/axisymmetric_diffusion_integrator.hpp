// Copyright (c) 2026 T. C. Raymond
// SPDX-License-Identifier: MIT

#pragma once

#include "mfem.hpp"
#include "constants.hpp"

// -----------------------------------------------------------------------------
// 1. Stiffness Integrator: -div( eps * grad(u) )
// -----------------------------------------------------------------------------
// Solves: Integral( eps * grad(u) . grad(v) * 2*pi*r * dr * dz )
//
// Implements the Zienkiewicz-Zhu flux/energy hooks (ComputeElementFlux /
// ComputeFluxEnergy) so MFEM's recovery-based error estimators
// (ZienkiewiczZhuEstimator) drive AMR with an axisymmetrically-consistent error.
// The pointwise flux is the physical density eps*grad(u); the energy integral
// carries the same 2*pi*r measure as AssembleElementMatrix so refinement
// concentrates by physical (r-z) field energy rather than a planar approximation.
class AxisymmetricDiffusionIntegrator : public mfem::BilinearFormIntegrator
{
private:
   mfem::Coefficient *Q;
   static constexpr double factor = Constants::TWO_PI;

public:
   explicit AxisymmetricDiffusionIntegrator(
      mfem::Coefficient &q, const mfem::IntegrationRule *ir = nullptr)
      : mfem::BilinearFormIntegrator(ir), Q(&q)
   {
      MFEM_ASSERT(Q != nullptr, "Coefficient cannot be null");
   }

   static const mfem::IntegrationRule &GetRule(
      const mfem::FiniteElement &trial_fe,
      const mfem::FiniteElement &test_fe,
      const mfem::ElementTransformation &Trans)
   {
      // Curved inverse mappings make the transformed-gradient integrand
      // rational, so this is a conservative order estimate.
      const int order = Trans.OrderGrad(&trial_fe)
         + Trans.OrderGrad(&test_fe) + Trans.Order();

      if (trial_fe.Space() == mfem::FunctionSpace::rQk)
      {
         return mfem::RefinedIntRules.Get(trial_fe.GetGeomType(), order);
      }
      return mfem::IntRules.Get(trial_fe.GetGeomType(), order);
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

      const mfem::IntegrationRule *ir = GetIntegrationRule(el, Trans);

      for (int i = 0; i < ir->GetNPoints(); i++)
      {
         const mfem::IntegrationPoint &ip = ir->IntPoint(i);
         Trans.SetIntPoint(&ip);
         Trans.Transform(ip, pos);

         double r = pos(0); // Radius is X

         // Weight = quad_weight * det(J) * (2 * pi * r) * epsilon
         w = ip.weight * Trans.Weight() * (factor * r) * Q->Eval(Trans, ip);

         el.CalcDShape(ip, dshape);
         // dN/ds * J^-1
         Mult(dshape, Trans.InverseJacobian(), dshapedxt);
         // Integral of grad u dot grad v r dr dz
         AddMult_a_AAt(w, dshapedxt, elmat);
      }
   }

   // Pointwise discrete flux flux = eps*grad(u), evaluated at the flux element's
   // nodes (or @p ir when supplied). Mirrors mfem::DiffusionIntegrator's
   // scalar-coefficient path, including the component-major layout
   // flux(fnd*j + i) the ZZ estimator expects. The flux is a physical density,
   // so it carries NO geometric measure here (the 2*pi*r weight lives in the
   // energy integral below, matching AssembleElementMatrix).
   void ComputeElementFlux(const mfem::FiniteElement &el,
                           mfem::ElementTransformation &Trans,
                           mfem::Vector &u,
                           const mfem::FiniteElement &fluxelem,
                           mfem::Vector &flux, bool with_coef = true,
                           const mfem::IntegrationRule *ir = nullptr) override
   {
      const int nd = el.GetDof();
      const int dim = el.GetDim();
      const int spaceDim = Trans.GetSpaceDim();

      mfem::DenseMatrix dshape(nd, dim);
      mfem::DenseMatrix invdfdx(dim, spaceDim);
      mfem::Vector vec(dim);
      mfem::Vector vecdxt(spaceDim);

      if (!ir)
      {
         ir = &fluxelem.GetNodes();
      }
      const int fnd = ir->GetNPoints();
      flux.SetSize(fnd * spaceDim);

      for (int i = 0; i < fnd; i++)
      {
         const mfem::IntegrationPoint &ip = ir->IntPoint(i);
         el.CalcDShape(ip, dshape);
         dshape.MultTranspose(u, vec);

         Trans.SetIntPoint(&ip);
         mfem::CalcInverse(Trans.Jacobian(), invdfdx);
         invdfdx.MultTranspose(vec, vecdxt);

         if (with_coef)
         {
            vecdxt *= Q->Eval(Trans, ip); // flux = eps * grad(u)
         }
         for (int j = 0; j < spaceDim; j++)
         {
            flux(fnd * j + i) = vecdxt(j);
         }
      }
   }

   // Energy norm of a flux expansion: Integral( eps * |flux|^2 * 2*pi*r ).
   // Mirrors mfem::DiffusionIntegrator (scalar Q) but applies the axisymmetric
   // 2*pi*r measure so the ZZ error indicator sqrt(energy) reflects the physical
   // r-z field. Anisotropic splitting (d_energy) is not supported.
   double ComputeFluxEnergy(const mfem::FiniteElement &fluxelem,
                            mfem::ElementTransformation &Trans,
                            mfem::Vector &flux,
                            mfem::Vector *d_energy = nullptr) override
   {
      const int nd = fluxelem.GetDof();
      const int spaceDim = Trans.GetSpaceDim();

      mfem::Vector shape(nd);
      mfem::Vector pointflux(spaceDim);
      mfem::Vector pos(2); // 2D axisymmetric

      const int order = 2 * fluxelem.GetOrder()
         + Trans.Order() + Trans.OrderW();
      const mfem::IntegrationRule *ir = GetIntRule();
      if (!ir)
      {
         ir = &mfem::IntRules.Get(fluxelem.GetGeomType(), order);
      }

      double energy = 0.0;
      if (d_energy) { *d_energy = 0.0; } // anisotropic estimation unsupported

      for (int i = 0; i < ir->GetNPoints(); i++)
      {
         const mfem::IntegrationPoint &ip = ir->IntPoint(i);
         Trans.SetIntPoint(&ip);
         fluxelem.CalcPhysShape(Trans, shape);

         pointflux = 0.0;
         for (int k = 0; k < spaceDim; k++)
         {
            for (int j = 0; j < nd; j++)
            {
               pointflux(k) += flux(k * nd + j) * shape(j);
            }
         }

         Trans.Transform(ip, pos);
         const double r = pos(0); // Radius is X

         // Axisymmetric measure: 2*pi*r, consistent with AssembleElementMatrix.
         const double w = Trans.Weight() * ip.weight * (factor * r);

         double e = (pointflux * pointflux);
         e *= Q->Eval(Trans, ip); // eps
         energy += w * e;
      }

      return energy;
   }

protected:
   const mfem::IntegrationRule *GetDefaultIntegrationRule(
      const mfem::FiniteElement &trial_fe,
      const mfem::FiniteElement &test_fe,
      const mfem::ElementTransformation &Trans) const override
   {
      return &GetRule(trial_fe, test_fe, Trans);
   }
};