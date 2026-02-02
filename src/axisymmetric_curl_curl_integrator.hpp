// Copyright (c) 2026 T. C. Raymond
// SPDX-License-Identifier: MIT

#pragma once

#include "mfem.hpp"
#include "constants.hpp"

/**
 * @brief Axisymmetric Curl-Curl Integrator for Magnetostatics
 * Solves: Integral( 1/mu * Curl(A) . Curl(v) * 2*pi*r * dr * dz )
 * Handles the r=0 singularity by substituting A/r -> dA/dr
 */
class AxisymmetricCurlCurlIntegrator : public mfem::BilinearFormIntegrator
{
private:
   mfem::Coefficient *Q; // Reluctivity (1/mu)
   static constexpr double factor = Constants::TWO_PI;
   static constexpr double r_tol = Constants::AXIS_TOLERANCE; // Tolerance for axis detection

public:
   AxisymmetricCurlCurlIntegrator(mfem::Coefficient &q) : Q(&q) 
   {
      MFEM_ASSERT(Q != nullptr, "Coefficient cannot be null");
   }

   void AssembleElementMatrix(const mfem::FiniteElement &el,
                               mfem::ElementTransformation &Trans,
                               mfem::DenseMatrix &elmat) override
   {
      int nd = el.GetDof();
      double w;

      elmat.SetSize(nd);
      elmat = 0.0;

      mfem::Vector shape(nd);
      mfem::DenseMatrix dshape(nd, el.GetDim());      // Gradient in Reference Space
      mfem::DenseMatrix dshape_phys(nd, el.GetDim()); // Gradient in Physical Space
      mfem::Vector pos(2); // 2D axisymmetric

      // Increase order slightly to capture 1/r curvature near axis
      int order = 2 * el.GetOrder() + Trans.OrderGrad(&el);
      const mfem::IntegrationRule *ir = &mfem::IntRules.Get(el.GetGeomType(), order);

      for (int i = 0; i < ir->GetNPoints(); i++)
      {
         const mfem::IntegrationPoint &ip = ir->IntPoint(i);
         Trans.SetIntPoint(&ip);
         Trans.Transform(ip, pos);
         
         double r = pos(0); 

         // 1. Material Property (1/mu)
         double nu = Q->Eval(Trans, ip);

         // 2. Integration Weight
         // w = alpha * weight_ip * det(J) * (2*pi*r)
         double geometry_factor = factor * r;
         
         // Safety: If r is effectively zero, the volume element vanishes,
         // but we still need the limit for the stiffness term stability.
         // However, standard FEM integration usually relies on r in the weight.
         // If r ~ 0, weight ~ 0, so contribution is minimal unless the term blows up.
         // We handle the blow-up below.
         w = ip.weight * Trans.Weight() * geometry_factor * nu;

         // 3. Shape Functions & Gradients
         el.CalcShape(ip, shape);
         el.CalcDShape(ip, dshape);
         Mult(dshape, Trans.InverseJacobian(), dshape_phys);

         // 4. Assemble Matrix
         if (r > r_tol)
         {
             // --- ROBUST SINGULARITY HANDLING (Far Field) ---
             // curl_z = (shape / r) + dN_dr
             for (int j = 0; j < nd; j++)
             {
                 // Physical derivatives of shape function j
                 double dNj_dr = dshape_phys(j, 0); // x-derivative
                 double dNj_dz = dshape_phys(j, 1); // y-derivative (mapped to z)

                 // Calculate the Azimuthal Curl components
                 // Curl(A)_r = -dA/dz
                 // Curl(A)_z = A/r + dA/dr

                 double curl_j_r = -dNj_dz;
                 double curl_j_z = (shape(j) / r) + dNj_dr;

                 // Exploit symmetry: only compute upper triangle
                 for (int k = j; k < nd; k++)
                 {
                     double dNk_dr = dshape_phys(k, 0);
                     double dNk_dz = dshape_phys(k, 1);

                     double curl_k_r = -dNk_dz;
                     double curl_k_z = (shape(k) / r) + dNk_dr;

                     // Dot Product of the two curl vectors
                     double val = curl_j_r * curl_k_r + curl_j_z * curl_k_z;

                     elmat(j, k) += w * val;
                     if (k != j) {
                         elmat(k, j) += w * val;
                     }
                 }
             }
         }
         else
         {
             // --- ROBUST SINGULARITY HANDLING (Near Field r->0) ---
             // Limit as r->0: A/r -> dA/dr
             // Total term becomes 2 * dA/dr
             for (int j = 0; j < nd; j++)
             {
                 double dNj_dr = dshape_phys(j, 0);
                 double dNj_dz = dshape_phys(j, 1);

                 double curl_j_r = -dNj_dz;
                 double curl_j_z = 2.0 * dNj_dr;

                 for (int k = j; k < nd; k++)
                 {
                     double dNk_dr = dshape_phys(k, 0);
                     double dNk_dz = dshape_phys(k, 1);

                     double curl_k_r = -dNk_dz;
                     double curl_k_z = 2.0 * dNk_dr;

                     double val = curl_j_r * curl_k_r + curl_j_z * curl_k_z;

                     elmat(j, k) += w * val;
                     if (k != j) {
                         elmat(k, j) += w * val;
                     }
                 }
             }
         }
      }
   }
};