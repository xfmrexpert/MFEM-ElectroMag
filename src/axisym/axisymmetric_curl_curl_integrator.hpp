// Copyright (c) 2026 T. C. Raymond
// SPDX-License-Identifier: MIT

#pragma once

#include <algorithm>
#include <cmath>
#include <limits>

#include "mfem.hpp"
#include "axisymmetric_measure.hpp"
#include "axisymmetric_field_relations.hpp"

/**
 * @brief Thread-safe axisymmetric curl-curl bilinear form integrator for magnetostatics
 *        with A = A_phi(r,z) e_phi.
 *
 * Assembles, with the full axisymmetric measure (see axisymmetric_measure.hpp):
 *
 *   integral nu * [dA/dz * dv/dz
 *                  + (dA/dr + A/r) * (dv/dr + v/r)] * 2*pi*r dr dz
 *
 * This follows from
 *
 *   curl(A_phi e_phi) = -dA_phi/dz e_r
 *                       + (dA_phi/dr + A_phi/r) e_z.
 *
 * IMPORTANT:
 *   Enforce the essential BC A_phi = 0 on the symmetry axis r = 0 (regularity).
 *   This is required whenever the domain closure reaches r = 0; annular domains
 *   (r_min > 0) need no axis condition at all.
 *
 * Notes:
 *   - This class is THREAD-SAFE: all scratch storage is local to AssembleElementMatrix().
 *   - Assembly does NOT clamp r. Standard interior quadrature keeps r > 0 even for
 *     elements touching the axis, so a zero radius there signals a bad custom rule
 *     or a bad mesh and is reported rather than papered over. The 1/r term must be
 *     kept exact per basis function: individual shape functions need not vanish on
 *     the axis, so no per-basis limit exists. Regularity is a property of the
 *     constrained solution and is delivered by essential BC elimination.
 *   - Flux/energy recovery, which does see the constrained solution, uses the exact
 *     limit B_z -> 2 dA/dr on the axis (see ComputeElementFlux). Both that path and
 *     MagneticFieldCoefficient call axisym::AxialFluxDensity with the mesh's
 *     scale-relative axis tolerance, so recovery, postprocessing, and geometry
 *     classification cannot disagree about what counts as the axis.
 */
class AxisymmetricCurlCurlIntegrator : public mfem::BilinearFormIntegrator
{
public:
   // @param axis_tolerance  Scale-relative radius below which flux recovery takes
   //                        the axis limit; use axisym::AxisGeometry::tolerance so the
   //                        integrator shares the mesh's axis policy.
   explicit AxisymmetricCurlCurlIntegrator(
      mfem::Coefficient &reluctivity,
      mfem::real_t axis_tolerance,
      const mfem::IntegrationRule *ir = nullptr)
      : mfem::BilinearFormIntegrator(ir), nu_(&reluctivity),
        axis_tolerance_(axis_tolerance)
   {
      MFEM_ASSERT(nu_ != nullptr, "Reluctivity coefficient cannot be null");
   }

   // Relative accuracy targeted for the non-polynomial 1/r term. This is a
   // quadrature-design target, not a field value, so it stays double: the
   // Bernstein-ellipse analysis in RadialExtraOrder is done in double
   // regardless of MFEM's storage precision, and 1e-10 is not representable
   // as a meaningful float target.
   static constexpr double kRadialQuadratureTolerance = 1.0e-10;

   // Ceiling on the extra points added for the 1/r term. Reached only for
   // elements whose inner radius is a tiny fraction of their radial width;
   // see GetRule for what that costs and why it is acceptable.
   //
   // Note this ceiling is not what actually limits the rule in practice:
   // kMaxPositiveWeightSimplexOrder below binds first on simplices, and by a
   // wide margin. The cap here only matters on tensor-product geometries.
   static constexpr int kMaxRadialExtraOrder = 120;

   // Highest simplex order for which MFEM tabulates a positive-weight rule.
   // Above this, IntRules.Get falls back to Grundmann-Moller, whose weights
   // alternate in sign and grow without bound. Those rules are unusable here:
   // the 1/r factor is largest exactly where the negative weights sit, so the
   // element matrix loses its positive definiteness to catastrophic
   // cancellation rather than merely losing digits. A lower-order rule with
   // positive weights is strictly the better failure mode, so GetRule clamps.
   //
   // Measured against MFEM 4.10 by walking IntRules.Get(TRIANGLE, o): the
   // first rule carrying a negative weight is order 26, and the first with a
   // point on the element boundary is order 16. Both are properties of the
   // rule tables, so this constant must be rechecked when MFEM is upgraded.
   static constexpr int kMaxPositiveWeightSimplexOrder = 25;

   // Smallest r_min/width at which the capped rule still meets
   // kRadialQuadratureTolerance. Measured: relative energy error stays near
   // 1e-11 down to this ratio, then degrades (about 1e-8 at 5e-3, 1e-5 at
   // 2e-3). Solvers use this to warn about under-resolved elements.
   //
   // This does NOT bracket the clamp above. Inverting RadialExtraOrder for
   // kRadialQuadratureTolerance = 1e-10 puts the order-26 crossing at
   // r_min/width ~= 0.21, about 20x looser than this ratio. An element can
   // therefore hit kMaxPositiveWeightSimplexOrder while still looking well
   // resolved by this measure. The two thresholds answer different questions:
   // this one is about accuracy, that one about whether a usable rule exists.
   static constexpr double kResolvedRadiusRatio = 1.0e-2;

   /**
    * @brief Additional integration order needed to resolve the 1/r term.
    *
    * The integrand splits into a polynomial part and the single non-polynomial
    * term N_j N_k / r. On an element spanning r in [a, a+h], mapping to the
    * reference interval x in [-1, 1] gives r = a + h(1+x)/2, so 1/r has a pole
    * at x0 = -(1 + 2s) with s = a/h. Gauss-Legendre applied to a function whose
    * nearest singularity lies at x0 converges geometrically like rho^-n, where
    *
    *     rho = |x0| + sqrt(x0^2 - 1)
    *
    * is the Bernstein-ellipse parameter. Achieving a relative accuracy eps
    * therefore needs roughly ln(1/eps)/ln(rho) points, which blows up as
    * s -> 0 because rho -> 1.
    *
    * This is why a fixed polynomial-order heuristic cannot work: it depends
    * only on the basis degree, while the actual difficulty is set by the
    * geometric ratio s. Measured convergence orders match this estimate closely
    * across s in [1e-3, 3].
    *
    * Returns 0 when the element touches the axis. There s = 0 and no finite
    * rule converges, because the exact integral of N_j N_k / r diverges
    * logarithmically for basis functions that do not vanish at r = 0. That
    * divergence is a property of individual basis functions, not of the
    * solution: axis regularity forces A_phi -> 0 linearly, and the essential
    * A_phi = 0 constraint removes exactly the offending directions. Spending
    * quadrature there would refine a quantity that elimination discards.
    */
   static int RadialExtraOrder(double min_radius, double radial_width)
   {
      if (!(radial_width > 0.0) || !(min_radius > 0.0)) { return 0; }

      const double s = min_radius / radial_width;
      const double x0 = 1.0 + 2.0 * s;
      const double rho = x0 + std::sqrt(x0 * x0 - 1.0);
      if (!(rho > 1.0)) { return kMaxRadialExtraOrder; }

      const double points =
         std::log(1.0 / kRadialQuadratureTolerance) / std::log(rho);
      if (!(points > 0.0)) { return 0; }
      if (points >= static_cast<double>(kMaxRadialExtraOrder))
      {
         return kMaxRadialExtraOrder;
      }
      return static_cast<int>(std::ceil(points));
   }

   // Radial extent of an element, sampled through its transformation so curved
   // geometry is respected.
   static void RadialExtent(const mfem::ElementTransformation &Trans,
                            mfem::real_t &min_radius, mfem::real_t &radial_width)
   {
      auto &T = const_cast<mfem::ElementTransformation &>(Trans);
      const mfem::IntegrationRule &vertices =
         *mfem::Geometries.GetVertices(T.GetGeometryType());

      mfem::real_t min_r = std::numeric_limits<mfem::real_t>::max();
      mfem::real_t max_r = std::numeric_limits<mfem::real_t>::lowest();
      mfem::Vector pos(T.GetSpaceDim());
      for (int i = 0; i < vertices.GetNPoints(); ++i)
      {
         const mfem::IntegrationPoint &ip = vertices.IntPoint(i);
         T.SetIntPoint(&ip);
         T.Transform(ip, pos);
         min_r = std::min(min_r, pos(0));
         max_r = std::max(max_r, pos(0));
      }

      min_radius = min_r;
      radial_width = max_r - min_r;
   }

   static const mfem::IntegrationRule &GetRule(
      const mfem::FiniteElement &trial_fe,
      const mfem::FiniteElement &test_fe,
      const mfem::ElementTransformation &Trans)
   {
      // Polynomial part: integrated exactly by the usual order estimate.
      const int gradient_order = Trans.OrderGrad(&trial_fe)
         + Trans.OrderGrad(&test_fe) + Trans.Order();
      const int radial_reaction_order = trial_fe.GetOrder()
         + test_fe.GetOrder() + Trans.Order() + Trans.OrderW();
      const int polynomial_order =
         std::max(gradient_order, radial_reaction_order);

      // Non-polynomial 1/r part: cost is set by the element's geometry, not by
      // the basis degree, so it must be added on top.
      mfem::real_t min_radius = 0.0;
      mfem::real_t radial_width = 0.0;
      RadialExtent(Trans, min_radius, radial_width);
      int order = polynomial_order
         + RadialExtraOrder(min_radius, radial_width);

      if (trial_fe.Space() == mfem::FunctionSpace::rQk)
      {
         return mfem::RefinedIntRules.Get(trial_fe.GetGeomType(), order);
      }

      // On simplices, refuse to cross into the negative-weight fallback rules;
      // see kMaxPositiveWeightSimplexOrder. Clamping silently under-integrates
      // the 1/r term, but it keeps the element matrix positive definite, which
      // is the property the solver actually depends on. Elements that reach
      // here are already flagged to the user by the kResolvedRadiusRatio check
      // in the solvers, so the accuracy loss is reported through that path.
      const mfem::Geometry::Type geom = trial_fe.GetGeomType();
      if ((geom == mfem::Geometry::TRIANGLE || geom == mfem::Geometry::TETRAHEDRON)
          && order > kMaxPositiveWeightSimplexOrder)
      {
         order = kMaxPositiveWeightSimplexOrder;
      }
      return mfem::IntRules.Get(geom, order);
   }

   void AssembleElementMatrix(const mfem::FiniteElement &el,
                              mfem::ElementTransformation &Trans,
                              mfem::DenseMatrix &elmat) override
   {
      const int nd  = el.GetDof();
      const int dim = el.GetDim();

      MFEM_VERIFY(dim == 2 && Trans.GetSpaceDim() == 2,
         "AxisymmetricCurlCurlIntegrator expects a 2D (r,z) finite element.");

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
         const mfem::real_t r = pos(0);

         // Interior quadrature keeps r > 0 even for elements touching the axis.
         // A non-positive radius here means an invalid mesh or a custom rule
         // with boundary points; clamping would silently deform the geometry.
         MFEM_VERIFY(r > 0.0,
            "AxisymmetricCurlCurlIntegrator evaluated at a non-positive radius (r = "
            << r << "). The 1/r term is singular there: use an interior "
            "integration rule and a mesh with non-negative radii.");

         const mfem::real_t nu = nu_->Eval(Trans, ip);

         // Axisymmetric weight: ip.weight * detJ * 2*pi*r * nu
         const mfem::real_t w = ip.weight * Trans.Weight()
            * Axisymmetric::Measure(r) * nu;

         el.CalcShape(ip, shape);
         el.CalcDShape(ip, dshape_ref);

         // Map row-oriented reference derivatives to physical derivatives
         // using MFEM's dshape_phys = dshape_ref * InvJ convention.
         Mult(dshape_ref, Trans.InverseJacobian(), dshape_phys);

         for (int j = 0; j < nd; j++)
         {
            const mfem::real_t Nj     = shape(j);
            const mfem::real_t dNj_dr = dshape_phys(j, 0);
            const mfem::real_t dNj_dz = dshape_phys(j, 1);

            for (int k = j; k < nd; k++)
            {
               const mfem::real_t Nk     = shape(k);
               const mfem::real_t dNk_dr = dshape_phys(k, 0);
               const mfem::real_t dNk_dz = dshape_phys(k, 1);

               const mfem::real_t Br_j = -dNj_dz;
               const mfem::real_t Bz_j = dNj_dr + Nj / r;

               const mfem::real_t Br_k = -dNk_dz;
               const mfem::real_t Bz_k = dNk_dr + Nk / r;

               const mfem::real_t val = Br_j * Br_k + Bz_j * Bz_k;

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

       MFEM_VERIFY(dim == 2 && space_dim == 2,
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
       // Inverse Jacobian: (dim x space_dim), NOT shaped like dshape.
       mfem::DenseMatrix inv_jacobian(dim, space_dim);
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
           mfem::CalcInverse(Trans.Jacobian(), inv_jacobian);
           inv_jacobian.MultTranspose(grad_ref, grad_phys);

           Trans.Transform(ip, pos);
           const mfem::real_t r = pos(0);
           const mfem::real_t A_phi = shape * u;

           // curl(A_phi e_phi) in the (r,z) component ordering, with the axis
           // limit applied under the shared scale-relative tolerance.
           mfem::real_t B_r = -grad_phys(1);
           mfem::real_t B_z = axisym::AxialFluxDensity(A_phi, grad_phys(0), r,
                                                 axis_tolerance_);

           if (with_coef)
           {
               const mfem::real_t nu = nu_->Eval(Trans, ip);
               B_r *= nu;
               B_z *= nu;
           }

           // MFEM vector GridFunction element data is component-major.
           flux(i) = B_r;
           flux(flux_nd + i) = B_z;
       }
   }

   // Returns mfem::real_t to match the base class declaration exactly; see the
   // note on AxisymmetricDiffusionIntegrator::ComputeFluxEnergy for why a
   // double here would silently stop overriding in a single-precision build.
   mfem::real_t ComputeFluxEnergy(const mfem::FiniteElement& flux_elem,
       mfem::ElementTransformation& Trans,
       mfem::Vector& flux,
       mfem::Vector* d_energy = nullptr) override
   {
       const int nd = flux_elem.GetDof();
       const int space_dim = Trans.GetSpaceDim();

       MFEM_VERIFY(space_dim == 2,
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

       mfem::real_t energy = 0.0;

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
           // No clamp needed: this term carries the measure r, never 1/r, so
           // the contribution simply vanishes on the axis.
           const mfem::real_t r = pos(0);
           const mfem::real_t nu = nu_->Eval(Trans, ip);
           const mfem::real_t weight = ip.weight * Trans.Weight()
               * Axisymmetric::Measure(r);

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
   mfem::real_t axis_tolerance_;
};

