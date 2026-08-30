// Copyright (c) 2026 T. C. Raymond
// SPDX-License-Identifier: MIT

#pragma once

#include <memory>
#include <string>
#include <vector>

#include "mfem.hpp"
#include "physics_solver.hpp"
#include "../axisym/axisymmetric_curl_curl_integrator.hpp"
#include "../axisym/magnetic_axis_boundary.hpp"

/**
 * @brief Base class for solvers formulated in the magnetic vector potential.
 *
 * Holds what the magnetostatic and magnetoquasistatic solvers share by virtue
 * of solving for the same unknown -- A_phi (axisymmetric) or A_z (planar) --
 * rather than by coincidence: the reluctivity coefficient, the curl-curl
 * stiffness term built from it, terminal current density, and the axis
 * regularity condition. None of this applies to an electrostatic run, which is
 * why it does not belong in PhysicsSolver.
 *
 * The solution field itself stays in the derived classes: magnetostatics holds
 * a real GridFunction, the time-harmonic solver a ComplexGridFunction.
 */
class MagneticSolver : public PhysicsSolver {
protected:
	// nu = 1/mu (reluctivity), keyed by mesh DOMAIN attribute. Assigned by each
	// derived Setup() via MaterialCoefficient(); everything below reads it.
	std::unique_ptr<mfem::PWConstCoefficient> nu_coeff;

	// Boundary attributes lying entirely on r = 0. Discovered here rather than
	// during geometric classification because only an A_phi formulation needs a
	// dedicated axis attribute; an electrostatic run on the same mesh does not.
	mfem::Array<int> axis_boundary;

	MagneticSolver(mfem::Mesh& m, const ProblemConfig& c) : PhysicsSolver(m, c) {}

	// Axis regularity, imposition half: locate the dedicated r=0 boundary
	// attribute and merge it into the essential marker so A_phi = 0 there. Call
	// after ess_bdr is built from the closure BCs and before BuildOperators().
	void ImposeAxisRegularity() {
		if (geometry != GeometryType::Axisymmetric) { return; }

		axis_boundary = axisym::FindAxisBoundaryMarker(mesh, axisymmetric_mesh);

		MFEM_VERIFY(ess_bdr.Size() == axis_boundary.Size(),
			"Axis boundary marker does not match the mesh boundary attributes.");
		for (int i = 0; i < ess_bdr.Size(); ++i) {
			ess_bdr[i] = ess_bdr[i] || axis_boundary[i];
		}
	}

	// Axis regularity, verification half: a nonzero Dirichlet value on the axis
	// contradicts the A_phi = 0 constraint imposed above. The constraint would
	// silently win, so the configuration is rejected instead. Requires the FE
	// space, so call after BuildOperators().
	void ValidateMagneticAxisBoundaryValues() const {
		if (geometry != GeometryType::Axisymmetric ||
			!axisymmetric_mesh.TouchesAxis()) return;

		MFEM_VERIFY(fespace,
			"Magnetic axis boundary validation requires a finite element space.");

		mfem::Array<int> axis_tdofs;
		fespace->GetEssentialTrueDofs(axis_boundary, axis_tdofs);
		mfem::Array<int> is_axis_tdof(fespace->GetTrueVSize());
		is_axis_tdof = 0;
		for (int i = 0; i < axis_tdofs.Size(); ++i) {
			is_axis_tdof[axis_tdofs[i]] = 1;
		}

		for (const auto& bc : closure_bcs) {
			if (bc.Condition.Type != BoundaryConditionType::Dirichlet ||
				bc.Condition.Value == 0.0) continue;

			mfem::Array<int> boundary_tdofs;
			fespace->GetEssentialTrueDofs(bc.Marker, boundary_tdofs);
			for (int i = 0; i < boundary_tdofs.Size(); ++i) {
				const int tdof = boundary_tdofs[i];
				MFEM_VERIFY(!is_axis_tdof[tdof],
					"Boundary group '" + bc.Condition.EntityGroupName +
					"' assigns a nonzero Dirichlet value at true DOF " +
					std::to_string(tdof) + " on the magnetic symmetry axis. "
					"Axis regularity requires A_phi = 0 at r = 0.");
			}
		}
	}

	// Stiffness term: axisymmetric curl-curl (nu * curl A * curl A, carrying the
	// 1/r factor) or planar diffusion (nu * grad A * grad A). A fresh instance is
	// returned each call so the solve's bilinear form and the AMR error estimator
	// can own separate copies.
	mfem::BilinearFormIntegrator* MakeStiffnessIntegrator() const {
		if (geometry == GeometryType::Axisymmetric) {
			return new AxisymmetricCurlCurlIntegrator(
				*nu_coeff, axisymmetric_mesh.tolerance);
		}
		else {
			return new mfem::DiffusionIntegrator(*nu_coeff);
		}
	}

	// Uniform current density I/area over the terminal's domain attributes,
	// laid out per mesh attribute for a PWConstCoefficient.
	mfem::Vector BuildTerminalCurrentDensity(
		const std::string& terminal_name, double current) const {
		const Terminal& term = config.Terminals.at(terminal_name);
		const EntityGroup& group = config.EntityGroups.at(term.EntityGroupName);
		const double area = CalculateRegionArea(group.AttributeIds);
		MFEM_VERIFY(area > 0.0,
			"Current terminal '" + terminal_name + "' has zero cross-section.");

		mfem::Vector current_density(mesh.attributes.Max());
		current_density = 0.0;
		const double density = current / area;
		for (int attr : group.AttributeIds) {
			if (attr > 0 && attr <= current_density.Size()) {
				current_density[attr - 1] = density;
			}
		}
		return current_density;
	}
};
