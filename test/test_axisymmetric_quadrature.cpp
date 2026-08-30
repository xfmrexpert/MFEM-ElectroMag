// Copyright (c) 2026 T. C. Raymond
// SPDX-License-Identifier: MIT

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include "../src/axisymmetric_curl_curl_integrator.hpp"
#include "../src/axisymmetric_diffusion_integrator.hpp"
#include "../src/axisymmetric_lf_integrator.hpp"
#include "../src/axisymmetric_mass_integrator.hpp"
#include "../src/axisymmetric_boundary_lf_integrator.hpp"
#include "../src/constants.hpp"

#include <algorithm>
#include <cmath>
#include <memory>

namespace {

std::unique_ptr<mfem::Mesh> MakeAxisymmetricMesh(bool curved)
{
   auto mesh = std::make_unique<mfem::Mesh>(
	  mfem::Mesh::MakeCartesian2D(1, 1, mfem::Element::QUADRILATERAL,
								  true, 1.0, 1.0));

   for (int v = 0; v < mesh->GetNV(); ++v)
   {
	  mesh->GetVertex(v)[0] += 1.0;
   }

   if (curved)
   {
	  mesh->SetCurvature(2);
	  mfem::VectorFunctionCoefficient deformation(
		 2, [](const mfem::Vector &x, mfem::Vector &mapped)
		 {
			mapped(0) = x(0);
			mapped(1) = x(1) + 0.1 * (x(0) - 1.0) * (2.0 - x(0));
		 });
	  mesh->Transform(deformation);
   }

   return mesh;
}

double RelativeDifference(const mfem::DenseMatrix &actual,
						  const mfem::DenseMatrix &reference)
{
   double difference_squared = 0.0;
   double reference_squared = 0.0;
   for (int row = 0; row < actual.Height(); ++row)
   {
	  for (int column = 0; column < actual.Width(); ++column)
	  {
		 const double difference = actual(row, column) - reference(row, column);
		 difference_squared += difference * difference;
		 reference_squared += reference(row, column) * reference(row, column);
	  }
   }

   return std::sqrt(difference_squared / reference_squared);
}

double RelativeDifference(const mfem::Vector &actual,
						  const mfem::Vector &reference)
{
   mfem::Vector difference(actual);
   difference -= reference;
   return difference.Norml2() / reference.Norml2();
}

void CheckAutomaticQuadrature(bool curved)
{
   auto mesh = MakeAxisymmetricMesh(curved);
   mfem::H1_FECollection collection(2, mesh->Dimension());
   mfem::FiniteElementSpace space(mesh.get(), &collection);

   const mfem::FiniteElement &element = *space.GetFE(0);
   mfem::ElementTransformation &transformation =
	  *mesh->GetElementTransformation(0);
   const int reference_order = 2 * element.GetOrder()
	  + transformation.Order() + transformation.OrderW() + 6;
   const mfem::IntegrationRule &reference_rule =
	  mfem::IntRules.Get(element.GetGeomType(), reference_order);
   mfem::ConstantCoefficient coefficient(1.0);

   mfem::DenseMatrix automatic_matrix;
   mfem::DenseMatrix reference_matrix;

   AxisymmetricMassIntegrator automatic_mass(coefficient);
	  AxisymmetricMassIntegrator reference_mass(coefficient, &reference_rule);
   automatic_mass.AssembleElementMatrix(element, transformation, automatic_matrix);
   reference_mass.AssembleElementMatrix(element, transformation, reference_matrix);
   REQUIRE(RelativeDifference(automatic_matrix, reference_matrix) < 1.0e-11);

   AxisymmetricDiffusionIntegrator automatic_diffusion(coefficient);
   AxisymmetricDiffusionIntegrator reference_diffusion(coefficient);
   reference_diffusion.SetIntRule(&reference_rule);
   automatic_diffusion.AssembleElementMatrix(
		element, transformation, automatic_matrix);
   reference_diffusion.AssembleElementMatrix(
		element, transformation, reference_matrix);
   REQUIRE(RelativeDifference(automatic_matrix, reference_matrix) < 1.0e-9);

   AxisymmetricCurlCurlIntegrator automatic_curl(coefficient, 0.0);
	 AxisymmetricCurlCurlIntegrator reference_curl(coefficient, 0.0);
   // The curl-curl rule now adapts to element geometry and can exceed the
   // shared reference order, so it needs a strictly finer reference.
   const mfem::IntegrationRule &curl_reference_rule =
	  mfem::IntRules.Get(element.GetGeomType(), 2 * reference_order + 40);
   reference_curl.SetIntRule(&curl_reference_rule);
   automatic_curl.AssembleElementMatrix(element, transformation, automatic_matrix);
   reference_curl.AssembleElementMatrix(element, transformation, reference_matrix);
	  // The 1/r reaction term is non-polynomial, but the automatic rule adds a
   // geometry-dependent order for it, so on this well-separated element
   // (r_min/width = 1) it is resolved to round-off rather than merely close.
   REQUIRE(RelativeDifference(automatic_matrix, reference_matrix) < 1.0e-11);

   mfem::Vector automatic_vector;
   mfem::Vector reference_vector;
   AxisymmetricLFIntegrator automatic_load(coefficient);
   AxisymmetricLFIntegrator reference_load(coefficient);
   reference_load.SetIntRule(&reference_rule);
   automatic_load.AssembleRHSElementVect(
	  element, transformation, automatic_vector);
   reference_load.AssembleRHSElementVect(
	  element, transformation, reference_vector);
   REQUIRE(RelativeDifference(automatic_vector, reference_vector) < 1.0e-11);
}

} // namespace

TEST_CASE("Axisymmetric quadrature accounts for element transformations",
		  "[axisymmetric][quadrature]")
{
   SECTION("affine geometry")
   {
	  CheckAutomaticQuadrature(false);
   }

   SECTION("curved quadratic geometry")
   {
	  CheckAutomaticQuadrature(true);
   }
}

namespace {

// One element spanning r in [s*h, s*h + h], z in [0, 1].
std::unique_ptr<mfem::Mesh> MakeRadialBand(double r_min, double width)
{
   auto mesh = std::make_unique<mfem::Mesh>(
	  mfem::Mesh::MakeCartesian2D(1, 1, mfem::Element::QUADRILATERAL,
								  true, width, 1.0));
   for (int v = 0; v < mesh->GetNV(); ++v)
   {
	  mesh->GetVertex(v)[0] += r_min;
   }
   return mesh;
}

double ConstantPotential(const mfem::Vector &) { return 1.0; }

// u^T K u for A_phi = 1, the magnetic energy the stiffness matrix represents.
double CurlCurlEnergy(mfem::Mesh &mesh)
{
   mfem::H1_FECollection collection(2, mesh.Dimension());
   mfem::FiniteElementSpace space(&mesh, &collection);
   mfem::ConstantCoefficient one(1.0);

   mfem::BilinearForm a(&space);
   a.AddDomainIntegrator(new AxisymmetricCurlCurlIntegrator(one, 0.0));
   a.Assemble();
   a.Finalize();

   mfem::GridFunction u(&space);
   mfem::FunctionCoefficient potential(ConstantPotential);
   u.ProjectCoefficient(potential);

   mfem::Vector Ku(u.Size());
   a.SpMat().Mult(u, Ku);
   return Ku * u;
}

} // namespace

// The 1/r term is non-polynomial and its quadrature difficulty is governed by
// the element's geometry, s = r_min/width, not by the basis degree: mapped to
// the reference interval, 1/r has a pole at -(1 + 2s), which approaches the
// integration interval as s -> 0. A basis-degree-only rule is therefore blind
// to the hard case, and a thin band close to the axis was integrated with tens
// of percent of error. A_phi = 1 does not vanish near the axis, so it exercises
// the term directly, and its energy has the closed form 2*pi*ln(b/a).
TEST_CASE("Near-axis annular curl-curl quadrature meets its accuracy target",
		  "[axisymmetric][quadrature]")
{
   const double width = 1.0;

   // kResolvedRadiusRatio is the documented limit of the capped rule.
   const double ratios[] = {
	  3.0, 1.0, 0.1, 0.03,
	  AxisymmetricCurlCurlIntegrator::kResolvedRadiusRatio};

   for (double s : ratios)
   {
	  const double a = s * width;
	  const double b = a + width;
	  auto mesh = MakeRadialBand(a, width);

	  const double exact = Constants::TWO_PI * std::log(b / a);
	  const double relative_error =
		 std::abs(CurlCurlEnergy(*mesh) - exact) / exact;

	  INFO("r_min/width = " << s);
	  REQUIRE(relative_error < 1.0e-10);
   }
}

// Below the documented ratio the rule is capped and only approximate. That is
// a deliberate cost bound, not an accident: the tensor-product point count
// grows quadratically in the order. Pinning the degradation here keeps the
// published limit honest and detects any silent change to the cap.
TEST_CASE("Curl-curl quadrature degrades only past its documented ratio",
		  "[axisymmetric][quadrature]")
{
   const double width = 1.0;
   const double s = 1.0e-3;
   const double a = s * width;
   auto mesh = MakeRadialBand(a, width);

   const double exact = Constants::TWO_PI * std::log((a + width) / a);
   const double relative_error =
	  std::abs(CurlCurlEnergy(*mesh) - exact) / exact;

   REQUIRE(s < AxisymmetricCurlCurlIntegrator::kResolvedRadiusRatio);
   // Far better than the ~40% the basis-degree-only rule produced here, but
   // short of the 1e-10 target, which is exactly what the warning reports.
   REQUIRE(relative_error > 1.0e-10);
   REQUIRE(relative_error < 1.0e-3);
}

TEST_CASE("Axisymmetric boundary load includes radial measure",
		  "[axisymmetric][quadrature][boundary]")
{
   auto mesh = MakeAxisymmetricMesh(false); // unit square shifted to 1 <= r <= 2
   mfem::H1_FECollection collection(1, mesh->Dimension());
   mfem::FiniteElementSpace space(mesh.get(), &collection);
   mfem::ConstantCoefficient unit_load(1.0);
   mfem::LinearForm load(&space);
   load.AddBoundaryIntegrator(new AxisymmetricBoundaryLFIntegrator(unit_load));
   load.Assemble();

   // Partition of unity gives the full measure integral_boundary 2*pi*r ds:
   // bottom + top + left + right = 3/2 + 3/2 + 1 + 2 = 6, times 2*pi.
   REQUIRE(load.Sum()
           == Catch::Approx(Constants::TWO_PI * 6.0).epsilon(1.0e-12));
}
