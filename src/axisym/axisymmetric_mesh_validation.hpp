// Copyright (c) 2026 T. C. Raymond
// SPDX-License-Identifier: MIT

#pragma once

#include <string>

#include "mfem.hpp"
#include "axis_geometry.hpp"

/**
 * @brief Input validation for axisymmetric (r,z) meshes.
 *
 * Distinct from classification: axis_geometry.hpp reports what the mesh IS,
 * while this decides whether that is acceptable input. Keeping the two apart
 * lets callers -- notably tests -- inspect a malformed mesh without tripping an
 * abort.
 */
namespace axisym {

/**
 * @brief Validate an axisymmetric mesh and return its radial extent.
 *
 * Aborts with a descriptive message when the mesh contains materially negative
 * radii, which is always a modelling error rather than something to clamp away.
 *
 * The result is cached by the solver at setup. It remains valid across
 * conforming AMR, which subdivides elements inside the existing bounding box
 * and therefore cannot change the tolerance or the axis relation.
 */
inline AxisGeometry ValidateMesh(mfem::Mesh &mesh)
{
   AxisGeometry info = InspectAxisGeometry(mesh);

   MFEM_VERIFY(info.relation != AxisRelation::NegativeRadius,
			   "Axisymmetric mesh extends to a negative radius (min r = "
			   + std::to_string(info.min_r)
			   + "). The r coordinate is the mesh x coordinate and must be "
				 "non-negative; check the mesh orientation or units.");

   return info;
}

} // namespace axisym
