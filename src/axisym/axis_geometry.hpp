// Copyright (c) 2026 T. C. Raymond
// SPDX-License-Identifier: MIT

#pragma once

#include <algorithm>
#include <cmath>
#include <limits>

#include "mfem.hpp"

/**
 * @brief Geometric classification of a 2D (r,z) mesh relative to the symmetry
 *        axis r = 0.
 *
 * This header answers exactly one question: where does the mesh sit relative to
 * the axis, and at what scale should "on the axis" be judged? It is deliberately
 * free of physics and of any boundary-condition policy -- what a solver does
 * about an axis-touching domain depends on the formulation, so that decision
 * lives with the solver that makes it.
 *
 * Design rules encoded here:
 *   - The radial extent is inspected ONCE, at solver setup, rather than clamping
 *     r at every quadrature point. Clamping silently deforms the geometry and
 *     hides genuinely invalid input.
 *   - Tolerances are relative to the mesh bounding box, so the same model works
 *     whether it is expressed in nanometres or kilometres.
 */
namespace axisym {

// Relative tolerance used to decide whether a coordinate "is" on the axis.
// Loose enough to absorb mesh-generator round-off, tight enough that a real
// gap between the domain and the axis is never mistaken for contact.
inline constexpr double kRelativeGeometryTolerance = 1.0e-10;

/// True when @p r is on the symmetry axis to within @p tolerance.
inline bool IsOnAxisGeometry(double r, double tolerance)
{
   return std::abs(r) <= tolerance;
}

/// How a mesh's radial extent relates to the symmetry axis.
enum class AxisRelation
{
   NegativeRadius, ///< Part of the domain lies at r < 0: not a valid (r,z) mesh.
   TouchesAxis,    ///< The closure of the domain reaches r = 0.
   Annular         ///< The domain is bounded away from the axis (r_min > 0).
};

/**
 * @brief Result of the radial scan of an axisymmetric mesh.
 *
 * Cached by the solver at setup. This stays valid across conforming AMR:
 * refinement adds vertices inside existing elements, so the bounding box -- and
 * therefore @c tolerance and @c relation -- cannot change.
 */
struct AxisGeometry
{
   double min_r = 0.0;
   double max_r = 0.0;
   /// Mesh-scale-relative tolerance for axis proximity tests.
   double tolerance = 0.0;
   AxisRelation relation = AxisRelation::Annular;

   [[nodiscard]] bool TouchesAxis() const
   {
	  return relation == AxisRelation::TouchesAxis;
   }

   /// True when @p r represents axis geometry at this mesh's scale.
   [[nodiscard]] bool IsOnAxisGeometry(double r) const
   {
	  return axisym::IsOnAxisGeometry(r, tolerance);
   }
};

/**
 * @brief Scan the physical radial coordinate over a 2D (r,z) mesh.
 *
 * Sampling is done through each element's transformation at the element's own
 * nodal points, so curved / high-order meshes are handled correctly. Checking
 * only mesh vertices is not sufficient: a curved edge can bulge across r = 0
 * while both of its endpoints stay non-negative.
 *
 * Pure classification: this never aborts on an unusual result, so callers can
 * inspect a negative-radius mesh and report it themselves.
 */
inline AxisGeometry InspectAxisGeometry(mfem::Mesh &mesh)
{
   MFEM_VERIFY(mesh.Dimension() == 2,
			   "Axisymmetric geometry requires a 2D (r,z) mesh.");

   AxisGeometry info;
   if (mesh.GetNE() == 0)
   {
	  return info;
   }

   double min_r = std::numeric_limits<double>::max();
   double max_r = std::numeric_limits<double>::lowest();
   double min_z = std::numeric_limits<double>::max();
   double max_z = std::numeric_limits<double>::lowest();

   mfem::Vector pos(mesh.SpaceDimension());
   const mfem::FiniteElementSpace *nodal_fes = mesh.GetNodalFESpace();

   for (int e = 0; e < mesh.GetNE(); ++e)
   {
	  mfem::ElementTransformation *trans = mesh.GetElementTransformation(e);

	  // Curved meshes are sampled at their geometry nodes; linear meshes fall
	  // back to the reference element's vertices.
	  const mfem::IntegrationRule &nodes =
		 nodal_fes
			? nodal_fes->GetFE(e)->GetNodes()
			: *mfem::Geometries.GetVertices(mesh.GetElementBaseGeometry(e));

	  for (int i = 0; i < nodes.GetNPoints(); ++i)
	  {
		 const mfem::IntegrationPoint &ip = nodes.IntPoint(i);
		 trans->SetIntPoint(&ip);
		 trans->Transform(ip, pos);

		 min_r = std::min(min_r, pos(0));
		 max_r = std::max(max_r, pos(0));
		 min_z = std::min(min_z, pos(1));
		 max_z = std::max(max_z, pos(1));
	  }
   }

   info.min_r = min_r;
   info.max_r = max_r;

   const double scale = std::max({ std::abs(max_r), std::abs(min_r),
								   max_z - min_z, 1.0e-300 });
   info.tolerance = kRelativeGeometryTolerance * scale;

   if (min_r < -info.tolerance)
   {
	  info.relation = AxisRelation::NegativeRadius;
   }
   else if (min_r <= info.tolerance)
   {
	  info.relation = AxisRelation::TouchesAxis;
   }
   else
   {
	  info.relation = AxisRelation::Annular;
   }

   return info;
}

} // namespace axisym
