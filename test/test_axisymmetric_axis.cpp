// Copyright (c) 2026 T. C. Raymond
// SPDX-License-Identifier: MIT

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include "axisymmetric_curl_curl_integrator.hpp"
#include "axisymmetric_mesh.hpp"
#include "magnetic_field_coefficient.hpp"

#include <cmath>
#include <memory>

namespace {

// Cartesian mesh over r in [r_min, r_min + r_extent], z in [0, 1].
std::unique_ptr<mfem::Mesh> MakeRadialMesh(int n, double r_min,
										   double r_extent = 1.0)
{
   auto mesh = std::make_unique<mfem::Mesh>(
	  mfem::Mesh::MakeCartesian2D(n, n, mfem::Element::QUADRILATERAL,
								  true, r_extent, 1.0));
   for (int v = 0; v < mesh->GetNV(); ++v)
   {
	  mesh->GetVertex(v)[0] += r_min;
   }
   return mesh;
}

// A_phi = r * (c0 + c1*z + c2*r^2) satisfies the axis regularity condition
// A_phi(0, z) = 0 exactly, and has the closed-form curl
//    B_r = -dA/dz          = -c1 * r
//    B_z = dA/dr + A/r     = 2*c0 + 2*c1*z + 4*c2*r^2
// whose axis value 2*c0 + 2*c1*z is precisely the 2*dA/dr limit.
constexpr double c0 = 0.7;
constexpr double c1 = -0.4;
constexpr double c2 = 0.25;

double ExactA(const mfem::Vector &x)
{
   const double r = x(0);
   const double z = x(1);
   return r * (c0 + c1 * z + c2 * r * r);
}

double ExactBr(double r, double /*z*/) { return -c1 * r; }

double ExactBz(double r, double z)
{
   return 2.0 * c0 + 2.0 * c1 * z + 4.0 * c2 * r * r;
}

} // namespace

TEST_CASE("Axisymmetric radial validation classifies the mesh",
		  "[axisymmetric][axis]")
{
   SECTION("domain touching the axis is detected")
   {
	  auto mesh = MakeRadialMesh(2, 0.0);
	  const axisym::MeshInfo info = axisym::ValidateMesh(*mesh);
	  REQUIRE(info.relation == axisym::AxisRelation::TouchesAxis);
	  REQUIRE(info.TouchesAxis());
	  REQUIRE(info.tolerance > 0.0);
	  REQUIRE(info.axis_boundary.Size() == mesh->bdr_attributes.Max());
	  REQUIRE(info.axis_boundary[3] == 1);
   }

   SECTION("annular domain needs no axis handling")
   {
	  auto mesh = MakeRadialMesh(2, 1.0);
	  const axisym::MeshInfo info = axisym::ValidateMesh(*mesh);
	  REQUIRE(info.relation == axisym::AxisRelation::Annular);
	  REQUIRE_FALSE(info.TouchesAxis());
	  REQUIRE(info.min_r == Catch::Approx(1.0));
	  REQUIRE(info.axis_boundary.Max() == 0);
   }

   SECTION("tolerance scales with the mesh, not with an absolute constant")
   {
	  auto small = MakeRadialMesh(2, 0.0, 1.0e-6);
	  auto large = MakeRadialMesh(2, 0.0, 1.0e6);
	  REQUIRE(axisym::InspectMesh(*small).tolerance <
			  axisym::InspectMesh(*large).tolerance);
   }

   SECTION("negative radius is rejected")
   {
	  auto mesh = MakeRadialMesh(2, -0.5);
	  REQUIRE(axisym::InspectMesh(*mesh).relation ==
			  axisym::AxisRelation::NegativeRadius);
   }

   SECTION("a curved edge bulging across the axis is caught")
   {
	  // Both endpoints of the left edge sit at r = 0, but the quadratic
	  // geometry dips to negative r in between. A vertex-only scan misses this.
	  auto mesh = MakeRadialMesh(1, 0.0);
	  mesh->SetCurvature(2);
	  mfem::VectorFunctionCoefficient deformation(
		 2, [](const mfem::Vector &x, mfem::Vector &mapped)
		 {
			mapped(0) = x(0) - 0.2 * (1.0 - x(0)) * x(1) * (1.0 - x(1));
			mapped(1) = x(1);
		 });
	  mesh->Transform(deformation);

	  REQUIRE(axisym::InspectMesh(*mesh).relation ==
			  axisym::AxisRelation::NegativeRadius);
   }
}

// The r -> 0 limit B_z -> 2 dA/dr is no longer a standalone function; it is
// inlined into AxisymmetricCurlCurlIntegrator::ComputeElementFlux and into
// MagneticFieldCoefficient::Eval. Both are exercised below against a
// manufactured solution on a mesh that reaches the axis, and each of those
// tests asserts that an on-axis point was actually sampled, so the limit
// cannot silently regress in either copy.

TEST_CASE("Axisymmetric curl-curl recovers B on a mesh touching the axis",
		  "[axisymmetric][axis]")
{
   auto mesh = MakeRadialMesh(4, 0.0);
	  const axisym::MeshInfo info = axisym::ValidateMesh(*mesh);
   REQUIRE(info.TouchesAxis());

   // Order 3 represents A_phi = r*(c0 + c1*z + c2*r^2) exactly.
   mfem::H1_FECollection collection(3, mesh->Dimension());
   mfem::FiniteElementSpace space(mesh.get(), &collection);
   mfem::GridFunction A(&space);
   mfem::FunctionCoefficient exact_A(ExactA);
   A.ProjectCoefficient(exact_A);

   mfem::ConstantCoefficient one(1.0);
   AxisymmetricCurlCurlIntegrator integrator(one, info.tolerance);

   // A Gauss-Lobatto basis is used deliberately: the default Gauss-Legendre
   // nodes are strictly interior, so they would never land on r = 0 and the
   // axis limit would go untested.
   mfem::L2_FECollection flux_collection(3, mesh->Dimension(),
                                         mfem::BasisType::GaussLobatto);
   mfem::FiniteElementSpace flux_space(mesh.get(), &flux_collection,
									   mesh->SpaceDimension());

   bool sampled_axis = false;

   for (int e = 0; e < space.GetNE(); ++e)
   {
	  mfem::Array<int> dofs;
	  space.GetElementDofs(e, dofs);
	  mfem::Vector u;
	  A.GetSubVector(dofs, u);

	  mfem::ElementTransformation &trans = *mesh->GetElementTransformation(e);
	  const mfem::FiniteElement &flux_fe = *flux_space.GetFE(e);
	  const mfem::IntegrationRule &nodes = flux_fe.GetNodes();

	  mfem::Vector flux;
	  integrator.ComputeElementFlux(*space.GetFE(e), trans, u, flux_fe, flux);

	  const int fnd = nodes.GetNPoints();
	  mfem::Vector pos(2);
	  for (int i = 0; i < fnd; ++i)
	  {
		 const mfem::IntegrationPoint &ip = nodes.IntPoint(i);
		 trans.SetIntPoint(&ip);
		 trans.Transform(ip, pos);
		 const double r = pos(0);
		 const double z = pos(1);

		 if (r == 0.0) { sampled_axis = true; }

		 REQUIRE(flux(i) == Catch::Approx(ExactBr(r, z)).margin(1.0e-9));
		 REQUIRE(flux(fnd + i) == Catch::Approx(ExactBz(r, z)).margin(1.0e-9));
	  }
   }

   // The regression this guards: on the axis the old clamp returned dA/dr
   // instead of 2*dA/dr, so the axis nodes must actually be visited.
   REQUIRE(sampled_axis);
}

TEST_CASE("MagneticFieldCoefficient agrees with the curl-curl flux on the axis",
		  "[axisymmetric][axis]")
{
   auto mesh = MakeRadialMesh(4, 0.0);
	  const axisym::MeshInfo info = axisym::ValidateMesh(*mesh);

   mfem::H1_FECollection collection(3, mesh->Dimension());
   mfem::FiniteElementSpace space(mesh.get(), &collection);
   mfem::GridFunction A(&space);
   mfem::FunctionCoefficient exact_A(ExactA);
   A.ProjectCoefficient(exact_A);

	  MagneticFieldCoefficient B(&A, info.tolerance);

   mfem::Vector value(2);
   mfem::Vector pos(2);
   bool sampled_axis = false;
   for (int e = 0; e < space.GetNE(); ++e)
   {
	  mfem::ElementTransformation &trans = *mesh->GetElementTransformation(e);
	  const mfem::IntegrationRule &nodes = space.GetFE(e)->GetNodes();
	  for (int i = 0; i < nodes.GetNPoints(); ++i)
	  {
		 const mfem::IntegrationPoint &ip = nodes.IntPoint(i);
		 trans.SetIntPoint(&ip);
		 trans.Transform(ip, pos);
		 B.Eval(value, trans, ip);

		 if (pos(0) == 0.0) { sampled_axis = true; }

		 REQUIRE(value(0) == Catch::Approx(ExactBr(pos(0), pos(1))).margin(1.0e-9));
		 REQUIRE(value(1) == Catch::Approx(ExactBz(pos(0), pos(1))).margin(1.0e-9));
	  }
   }

   // Guard the coverage itself: without this the axis branch could stop being
   // reached and every remaining assertion would still pass.
   REQUIRE(sampled_axis);
}

TEST_CASE("Near-zero axis coordinates recover the same limit as exact zero",
          "[axisymmetric][axis]")
{
   // Mesh generators emit axis coordinates with round-off, so a mesh that
   // geometry classification calls "touching the axis" can have r = 1e-18
   // rather than exactly 0. An exact r == 0.0 test in field recovery would then
   // evaluate A_phi / r at a near-zero radius, which is unbounded rather than
   // merely inaccurate. Recovery must use the same scale-relative axis policy
   // as the classification.
   auto mesh = MakeRadialMesh(4, 0.0);

   // Perturb only the axis vertices, by far less than the mesh tolerance.
   const double perturbation = 1.0e-18;
   for (int v = 0; v < mesh->GetNV(); ++v)
   {
      double *vertex = mesh->GetVertex(v);
      if (vertex[0] == 0.0) { vertex[0] = perturbation; }
   }

   const axisym::MeshInfo info = axisym::ValidateMesh(*mesh);
   REQUIRE(info.TouchesAxis());
   REQUIRE(info.min_r > 0.0);
   REQUIRE(info.min_r < info.tolerance);

   mfem::H1_FECollection collection(3, mesh->Dimension());
   mfem::FiniteElementSpace space(mesh.get(), &collection);
   mfem::GridFunction A(&space);
   mfem::FunctionCoefficient exact_A(ExactA);
   A.ProjectCoefficient(exact_A);

   mfem::ConstantCoefficient one(1.0);
   AxisymmetricCurlCurlIntegrator integrator(one, info.tolerance);
   MagneticFieldCoefficient B(&A, info.tolerance);

   mfem::L2_FECollection flux_collection(3, mesh->Dimension(),
                                         mfem::BasisType::GaussLobatto);
   mfem::FiniteElementSpace flux_space(mesh.get(), &flux_collection,
                                       mesh->SpaceDimension());

   bool sampled_axis = false;
   mfem::Vector value(2);
   mfem::Vector pos(2);

   for (int e = 0; e < space.GetNE(); ++e)
   {
      mfem::Array<int> dofs;
      space.GetElementDofs(e, dofs);
      mfem::Vector u;
      A.GetSubVector(dofs, u);

      mfem::ElementTransformation &trans = *mesh->GetElementTransformation(e);
      const mfem::FiniteElement &flux_fe = *flux_space.GetFE(e);
      const mfem::IntegrationRule &nodes = flux_fe.GetNodes();

      mfem::Vector flux;
      integrator.ComputeElementFlux(*space.GetFE(e), trans, u, flux_fe, flux);

      const int fnd = nodes.GetNPoints();
      for (int i = 0; i < fnd; ++i)
      {
         const mfem::IntegrationPoint &ip = nodes.IntPoint(i);
         trans.SetIntPoint(&ip);
         trans.Transform(ip, pos);
         const double r = pos(0);
         const double z = pos(1);

         if (info.IsOnAxisGeometry(r)) { sampled_axis = true; }

         // The perturbation is negligible physically, so the exact-zero values
         // remain the reference. Under an exact r == 0.0 test the axis nodes
         // would instead produce values of order A_phi / 1e-18.
         REQUIRE(flux(i) == Catch::Approx(ExactBr(r, z)).margin(1.0e-9));
         REQUIRE(flux(fnd + i) == Catch::Approx(ExactBz(r, z)).margin(1.0e-9));

         B.Eval(value, trans, ip);
         REQUIRE(value(0) == Catch::Approx(ExactBr(r, z)).margin(1.0e-9));
         REQUIRE(value(1) == Catch::Approx(ExactBz(r, z)).margin(1.0e-9));
      }
   }

   REQUIRE(sampled_axis);
}

