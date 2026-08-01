// Copyright (c) 2026 T. C. Raymond
// SPDX-License-Identifier: MIT

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include "../src/axisymmetric_curl_curl_integrator.hpp"
#include "../src/axisymmetric_diffusion_integrator.hpp"
#include "../src/axisymmetric_lf_integrator.hpp"
#include "../src/axisymmetric_mass_integrator.hpp"
#include "../src/axisymmetric_boundary_lf_integrator.hpp"

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

   AxisymmetricCurlCurlIntegrator automatic_curl(coefficient);
	 AxisymmetricCurlCurlIntegrator reference_curl(coefficient);
   reference_curl.SetIntRule(&reference_rule);
   automatic_curl.AssembleElementMatrix(element, transformation, automatic_matrix);
   reference_curl.AssembleElementMatrix(element, transformation, reference_matrix);
	  // The 1/r reaction term is non-polynomial, so automatic quadrature should
   // converge closely to over-integration without being expected to match it.
   REQUIRE(RelativeDifference(automatic_matrix, reference_matrix) < 1.0e-5);

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

   // Partition of unity gives integral_boundary r ds:
   // bottom + top + left + right = 3/2 + 3/2 + 1 + 2 = 6.
   REQUIRE(load.Sum() == Catch::Approx(6.0).epsilon(1.0e-12));
}
