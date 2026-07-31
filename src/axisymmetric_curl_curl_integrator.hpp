// Copyright (c) 2026 T. C. Raymond
// SPDX-License-Identifier: MIT

#pragma once

#include <algorithm>
#include <limits>

#include "mfem.hpp"

/**
 * @brief Thread-safe axisymmetric curl-curl bilinear form integrator for magnetostatics
 *        with A = A_phi(r,z) e_phi.
 *
 * Assembles (up to a global 2π factor, omitted consistently):
 *
 *   ∫ ν [ (∂A/∂r)(∂v/∂r) + (∂A/∂z)(∂v/∂z) + (A/r)(v/r) ] * r dr dz
 *
 * where ν = 1/μ.
 *
 * IMPORTANT:
 *   Enforce the essential BC A_phi = 0 on the symmetry axis r = 0 (regularity).
 *
 * Notes:
 *   - This class is THREAD-SAFE: all scratch storage is local to AssembleElementMatrix().
 *   - A tiny r clamp is used only as a safety net against pathological quadrature/mappings.
 *     Proper behavior near the axis should come from essential BC elimination.
 */
class AxisymmetricCurlCurlIntegrator : public mfem::BilinearFormIntegrator
{
public:
   explicit AxisymmetricCurlCurlIntegrator(
      mfem::Coefficient &reluctivity,
      const mfem::IntegrationRule *ir = nullptr)
      : mfem::BilinearFormIntegrator(ir), nu_(&reluctivity)
   {
      MFEM_ASSERT(nu_ != nullptr, "Reluctivity coefficient cannot be null");
   }

   static const mfem::IntegrationRule &GetRule(
      const mfem::FiniteElement &trial_fe,
      const mfem::FiniteElement &test_fe,
      const mfem::ElementTransformation &Trans)
   {
      // The transformed-gradient and 1/r contributions have different order
      // requirements. The latter is non-polynomial, so use the larger
      // geometry-aware estimate rather than claiming exact integration.
      const int gradient_order = Trans.OrderGrad(&trial_fe)
         + Trans.OrderGrad(&test_fe) + Trans.Order();
      const int radial_reaction_order = trial_fe.GetOrder()
         + test_fe.GetOrder() + Trans.Order() + Trans.OrderW();
      const int order = std::max(gradient_order, radial_reaction_order);

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
      const int nd  = el.GetDof();
      const int dim = el.GetDim();

      MFEM_ASSERT(dim == 2, "AxisymmetricCurlCurlIntegrator expects a 2D (r,z) finite element.");

      elmat.SetSize(nd);
      elmat = 0.0;

      // Thread-safe scratch (local)
      mfem::Vector      shape(nd);
      mfem::DenseMatrix dshape_ref(nd, dim);
      mfem::DenseMatrix dshape_phys(nd, dim);
      mfem::Vector      pos(dim);

      const mfem::IntegrationRule *ir = GetIntegrationRule(el, Trans);

      for (int i = 0; i < ir->GetNPoints(); i++)
      {
         const mfem::IntegrationPoint &ip = ir->IntPoint(i);
         Trans.SetIntPoint(&ip);

         // Physical coordinates: pos(0)=r, pos(1)=z
         Trans.Transform(ip, pos);
         const double r_raw = pos(0);

         // Safety net only: avoids division by zero in pathological cases.
         // Using a larger tolerance (e.g. 1e-12) prevents Infinity in the matrix
         // which can cause solver crashes. (1e-12)^(-2) = 1e24 which fits in double.
         const double r = std::max(r_raw, 1e-12);

         const double nu = nu_->Eval(Trans, ip);

         // Axisymmetric weight (global 2π omitted): ip.weight * detJ * r * nu
         const double w = ip.weight * Trans.Weight() * r * nu;

         el.CalcShape(ip, shape);
         el.CalcDShape(ip, dshape_ref);

         // Map reference derivatives -> physical derivatives.
         // This assumes MFEM convention: dshape_phys = dshape_ref * InvJ.
         // If your unit test shows transpose, switch to InvJ^T.
         Mult(dshape_ref, Trans.InverseJacobian(), dshape_phys);

         for (int j = 0; j < nd; j++)
         {
            const double Nj     = shape(j);
            const double dNj_dr = dshape_phys(j, 0);
            const double dNj_dz = dshape_phys(j, 1);

            for (int k = j; k < nd; k++)
            {
               const double Nk     = shape(k);
               const double dNk_dr = dshape_phys(k, 0);
               const double dNk_dz = dshape_phys(k, 1);

               const double Br_j = -dNj_dz;
               const double Bz_j = dNj_dr + Nj / r;

               const double Br_k = -dNk_dz;
               const double Bz_k = dNk_dr + Nk / r;

               const double val = Br_j * Br_k + Bz_j * Bz_k;

               elmat(j, k) += w * val;
            }
         }
      }

      // Symmetrize
      for (int j = 0; j < nd; j++)
      {
         for (int k = 0; k < j; k++)
         {
            elmat(j, k) = elmat(k, j);
         }
      }
   }

   void ComputeElementFlux(const mfem::FiniteElement& el,
       mfem::ElementTransformation& Trans,
       mfem::Vector& u,
       const mfem::FiniteElement& flux_elem,
       mfem::Vector& flux,
       bool with_coef = false,
       const mfem::IntegrationRule* ir = nullptr) override
   {
       const int nd = el.GetDof();
       const int dim = el.GetDim();
       const int space_dim = Trans.GetSpaceDim();

       MFEM_ASSERT(dim == 2 && space_dim == 2,
           "AxisymmetricCurlCurlIntegrator expects a 2D (r,z) mesh.");
       MFEM_ASSERT(u.Size() == nd,
           "Element solution has an unexpected size.");

       if (!ir)
       {
           ir = &flux_elem.GetNodes();
       }

       const int flux_nd = ir->GetNPoints();
       flux.SetSize(flux_nd * space_dim);

       mfem::Vector shape(nd);
       mfem::DenseMatrix dshape_ref(nd, dim);
       mfem::DenseMatrix dshape_phys(nd, space_dim);
       mfem::Vector grad_ref(dim);
       mfem::Vector grad_phys(space_dim);
       mfem::Vector pos(space_dim);

       for (int i = 0; i < flux_nd; ++i)
       {
           const mfem::IntegrationPoint& ip = ir->IntPoint(i);
           Trans.SetIntPoint(&ip);

           el.CalcShape(ip, shape);
           el.CalcDShape(ip, dshape_ref);

           dshape_ref.MultTranspose(u, grad_ref);
           mfem::CalcInverse(Trans.Jacobian(), dshape_phys);
           dshape_phys.MultTranspose(grad_ref, grad_phys);

           Trans.Transform(ip, pos);
           const double r = std::max(pos(0), 1e-12);
           const double A_phi = shape * u;

           // curl(A_phi e_phi) in the (r,z) component ordering.
           double B_r = -grad_phys(1);
           double B_z = grad_phys(0) + A_phi / r;

           if (with_coef)
           {
               const double nu = nu_->Eval(Trans, ip);
               B_r *= nu;
               B_z *= nu;
           }

           // MFEM vector GridFunction element data is component-major.
           flux(i) = B_r;
           flux(flux_nd + i) = B_z;
       }
   }

   double ComputeFluxEnergy(const mfem::FiniteElement& flux_elem,
       mfem::ElementTransformation& Trans,
       mfem::Vector& flux,
       mfem::Vector* d_energy = nullptr) override
   {
       const int nd = flux_elem.GetDof();
       const int space_dim = Trans.GetSpaceDim();

       MFEM_ASSERT(space_dim == 2,
           "AxisymmetricCurlCurlIntegrator expects a 2D (r,z) mesh.");
       MFEM_ASSERT(flux.Size() == nd * space_dim,
           "Flux vector has an unexpected size.");

       mfem::Vector shape(nd);
       mfem::Vector point_flux(space_dim);
       mfem::Vector pos(space_dim);

       const int order = 2 * flux_elem.GetOrder()
           + Trans.Order() + Trans.OrderW();

       const mfem::IntegrationRule* ir = GetIntRule();
       if (!ir)
       {
           ir = &mfem::IntRules.Get(flux_elem.GetGeomType(), order);
       }

       if (d_energy)
       {
           d_energy->SetSize(0);
       }

       double energy = 0.0;

       for (int i = 0; i < ir->GetNPoints(); ++i)
       {
           const mfem::IntegrationPoint& ip = ir->IntPoint(i);
           Trans.SetIntPoint(&ip);

           flux_elem.CalcPhysShape(Trans, shape);

           point_flux = 0.0;
           for (int component = 0; component < space_dim; ++component)
           {
               for (int j = 0; j < nd; ++j)
               {
                   point_flux(component) +=
                       flux(component * nd + j) * shape(j);
               }
           }

           Trans.Transform(ip, pos);
           const double r = std::max(pos(0), 1e-12);
           const double nu = nu_->Eval(Trans, ip);
           const double weight = ip.weight * Trans.Weight() * r;

           energy += weight * nu * (point_flux * point_flux);
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

private:
   mfem::Coefficient *nu_;
};

