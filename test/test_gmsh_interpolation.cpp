// Copyright (c) 2026 T. C. Raymond
// SPDX-License-Identifier: MIT
//
// Verifies the $InterpolationScheme emitted by gmsh_results_writer: the shape
// functions shipped in the file must actually reconstruct the field. These
// tests exercise the same code path a downstream consumer (the C# results
// pipeline) would implement from the file alone, so a regression here means a
// consumer would silently read wrong field values.

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>
#include <vector>

#include "mfem.hpp"
#include "io/gmsh_results_writer.hpp"

using gmsh_results::detail::BuildInterpScheme;
using gmsh_results::detail::GetHoLayout;
using gmsh_results::detail::InterpScheme;

namespace {

// Evaluate shape function i at (u, v) exactly the way a consumer would:
// sum over monomial terms of coeff * u^a * v^b.
double EvalShape(const InterpScheme& s, int i, double u, double v) {
	double acc = 0.0;
	for (size_t j = 0; j < s.exponents.size(); ++j) {
		acc += s.coeffs[i][j] * std::pow(u, s.exponents[j][0]) *
								std::pow(v, s.exponents[j][1]);
	}
	return acc;
}

// Reconstruct a field from nodal values, as a consumer would.
double Interpolate(const InterpScheme& s, const std::vector<double>& nodal,
				   double u, double v) {
	double acc = 0.0;
	for (size_t i = 0; i < nodal.size(); ++i) {
		acc += nodal[i] * EvalShape(s, i, u, v);
	}
	return acc;
}

} // namespace

TEST_CASE("Interpolation scheme satisfies the Lagrange property", "[gmsh][interp]") {
	// phi_i(node_k) == delta_ik. This is the defining property; if it fails,
	// nodal values are not the coefficients the consumer thinks they are.
	for (auto geom : { mfem::Geometry::TRIANGLE, mfem::Geometry::SQUARE }) {
		for (int order = 1; order <= 4; ++order) {
			const InterpScheme s = BuildInterpScheme(geom, order);
			const auto& layout = GetHoLayout(geom, order);
			const int n = static_cast<int>(layout.ref.size());

			for (int i = 0; i < n; ++i) {
				for (int k = 0; k < n; ++k) {
					const double got =
						EvalShape(s, i, layout.ref[k][0], layout.ref[k][1]);
					const double want = (i == k) ? 1.0 : 0.0;
					REQUIRE_THAT(got, Catch::Matchers::WithinAbs(want, 1e-9));
				}
			}
		}
	}
}

TEST_CASE("Interpolation scheme forms a partition of unity", "[gmsh][interp]") {
	// Shape functions must sum to 1 everywhere, otherwise a constant field
	// would not reconstruct as that constant.
	const double pts[][2] = {
		{ 0.1, 0.1 }, { 0.25, 0.5 }, { 0.5, 0.25 }, { 0.33, 0.33 }
	};
	for (auto geom : { mfem::Geometry::TRIANGLE, mfem::Geometry::SQUARE }) {
		for (int order = 1; order <= 4; ++order) {
			const InterpScheme s = BuildInterpScheme(geom, order);
			const int n = static_cast<int>(s.coeffs.size());
			for (const auto& p : pts) {
				double sum = 0.0;
				for (int i = 0; i < n; ++i) { sum += EvalShape(s, i, p[0], p[1]); }
				REQUIRE_THAT(sum, Catch::Matchers::WithinAbs(1.0, 1e-9));
			}
		}
	}
}

TEST_CASE("Interpolation scheme reproduces polynomials exactly", "[gmsh][interp]") {
	// An order-p basis must reproduce any polynomial of degree <= p exactly.
	// This is the property streamline integration depends on: sampling BETWEEN
	// nodes has to be right, not just at nodes.
	auto f = [](int deg, double u, double v) {
		// A polynomial of total degree exactly `deg`, with mixed terms.
		return 1.0 + 2.0 * u - 3.0 * v + std::pow(u + 0.5 * v, deg);
	};

	for (auto geom : { mfem::Geometry::TRIANGLE, mfem::Geometry::SQUARE }) {
		for (int order = 1; order <= 4; ++order) {
			const InterpScheme s = BuildInterpScheme(geom, order);
			const auto& layout = GetHoLayout(geom, order);

			std::vector<double> nodal;
			nodal.reserve(layout.ref.size());
			for (const auto& r : layout.ref) {
				nodal.push_back(f(order, r[0], r[1]));
			}

			// Sample at interior points that are NOT nodes.
			const double pts[][2] = {
				{ 0.137, 0.221 }, { 0.4, 0.35 }, { 0.05, 0.6 }, { 0.29, 0.11 }
			};
			for (const auto& p : pts) {
				const double got  = Interpolate(s, nodal, p[0], p[1]);
				const double want = f(order, p[0], p[1]);
				REQUIRE_THAT(got, Catch::Matchers::WithinAbs(want, 1e-8));
			}
		}
	}
}

TEST_CASE("Interpolation scheme matches MFEM shape functions", "[gmsh][interp]") {
	// The strongest check: the basis written to the file must agree with the
	// basis MFEM used to compute the solution, accounting for the Gmsh node
	// permutation. If these disagree, exported values and the exported basis
	// describe different fields.
	for (auto geom : { mfem::Geometry::TRIANGLE, mfem::Geometry::SQUARE }) {
		const mfem::Element::Type etype = (geom == mfem::Geometry::TRIANGLE)
			? mfem::Element::TRIANGLE : mfem::Element::QUADRILATERAL;
		mfem::Mesh mesh = mfem::Mesh::MakeCartesian2D(1, 1, etype);

		for (int order = 1; order <= 4; ++order) {
			mfem::H1_FECollection fec(order, 2, mfem::BasisType::ClosedUniform);
			mfem::FiniteElementSpace fes(&mesh, &fec);
			const mfem::FiniteElement& fe = *fes.GetFE(0);

			const InterpScheme s = BuildInterpScheme(geom, order);
			const auto& layout = GetHoLayout(geom, order);
			const auto perm =
				gmsh_results::detail::BuildDofPermutation(fe, layout);

			const double pts[][2] = { { 0.2, 0.3 }, { 0.45, 0.15 } };
			for (const auto& p : pts) {
				mfem::IntegrationPoint ip;
				ip.Set2(p[0], p[1]);
				mfem::Vector mfem_shape(fe.GetDof());
				fe.CalcShape(ip, mfem_shape);

				for (int k = 0; k < static_cast<int>(layout.ref.size()); ++k) {
					const double ours = EvalShape(s, k, p[0], p[1]);
					const double theirs = mfem_shape(perm[k]);
					REQUIRE_THAT(ours,
								 Catch::Matchers::WithinAbs(theirs, 1e-9));
				}
			}
		}
	}
}
