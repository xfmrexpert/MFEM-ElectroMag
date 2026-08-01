// Copyright (c) 2026 T. C. Raymond
// SPDX-License-Identifier: MIT

#pragma once

#include <algorithm>
#include <cmath>
#include <limits>
#include <string>

#include "mfem.hpp"

/**
  * @brief Setup-time validation and classification of an axisymmetric (r,z)
 *        mesh, including discovery of a dedicated symmetry-axis boundary.
 *
 * Design rules encoded here:
 *   - The radial extent of a mesh is inspected ONCE, at solver setup, rather
 *     than clamping r at every quadrature point. Clamping silently deforms the
 *     geometry and hides genuinely invalid input.
 *   - Tolerances are relative to the mesh bounding box, so the same model works
 *     whether it is expressed in nanometres or kilometres.
 *   - A domain is only subject to axis regularity handling if its closure
 *     actually reaches r = 0. Annular domains (r_min > 0) need nothing special.
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

/// Result of the setup-time radial scan of an axisymmetric mesh.
struct MeshInfo
{
   double min_r = 0.0;
   double max_r = 0.0;
   /// Mesh-scale-relative tolerance for axis proximity tests.
   double tolerance = 0.0;
   AxisRelation relation = AxisRelation::Annular;
   mfem::Array<int> axis_boundary;

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
 */
inline MeshInfo InspectMesh(mfem::Mesh &mesh)
{
   MFEM_VERIFY(mesh.Dimension() == 2,
			   "Axisymmetric geometry requires a 2D (r,z) mesh.");

	  MeshInfo info;
   const int num_bdr_attributes = mesh.bdr_attributes.Size()
	  ? mesh.bdr_attributes.Max() : 0;
   info.axis_boundary.SetSize(num_bdr_attributes);
   info.axis_boundary = 0;

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

	  if (!info.TouchesAxis())
   {
	  return info;
   }

   MFEM_VERIFY(num_bdr_attributes > 0,
			   "The axisymmetric domain reaches r = 0 but the mesh has no "
			   "boundary attributes. Tag the axis as its own boundary attribute "
			   "so A_phi = 0 can be enforced there.");

   mfem::Array<int> on_axis_count(num_bdr_attributes);
   mfem::Array<int> total_count(num_bdr_attributes);
   on_axis_count = 0;
   total_count = 0;

   const mfem::FiniteElementSpace *boundary_nodal_fes = mesh.GetNodalFESpace();
   for (int be = 0; be < mesh.GetNBE(); ++be)
   {
	  const int attr = mesh.GetBdrAttribute(be);
	  if (attr <= 0 || attr > num_bdr_attributes) { continue; }
	  total_count[attr - 1]++;

	  mfem::ElementTransformation *trans = mesh.GetBdrElementTransformation(be);
	  const mfem::IntegrationRule &nodes =
		 boundary_nodal_fes
			? boundary_nodal_fes->GetBE(be)->GetNodes()
			: *mfem::Geometries.GetVertices(mesh.GetBdrElementGeometry(be));

	  bool on_axis = true;
	  for (int i = 0; i < nodes.GetNPoints() && on_axis; ++i)
	  {
		 const mfem::IntegrationPoint &ip = nodes.IntPoint(i);
		 trans->SetIntPoint(&ip);
		 trans->Transform(ip, pos);
		 on_axis = info.IsOnAxisGeometry(pos(0));
	  }

	  if (on_axis) { on_axis_count[attr - 1]++; }
   }

   bool found_axis = false;
   for (int i = 0; i < num_bdr_attributes; ++i)
   {
	  if (total_count[i] > 0 && on_axis_count[i] == total_count[i])
	  {
		 info.axis_boundary[i] = 1;
		 found_axis = true;
	  }
	  else if (on_axis_count[i] > 0)
	  {
		 MFEM_ABORT("Boundary attribute " + std::to_string(i + 1) + " has "
				   + std::to_string(on_axis_count[i]) + " of "
				   + std::to_string(total_count[i])
				   + " elements on the symmetry axis. Give the axis its own "
					 "boundary attribute so A_phi = 0 is not applied to "
					 "unrelated boundaries.");
	  }
   }

   MFEM_VERIFY(found_axis,
			   "The axisymmetric domain reaches r = 0 but no boundary attribute "
			   "lies entirely on the symmetry axis. Tag the axis as its own "
			   "boundary attribute so A_phi = 0 can be enforced there.");

   return info;
}

/**
 * @brief Validate an axisymmetric mesh and return its radial extent.
 *
 * Aborts with a descriptive message when the mesh contains materially negative
 * radii, which is always a modelling error rather than something to clamp away.
 */
inline MeshInfo ValidateMesh(mfem::Mesh &mesh)
{
	  MeshInfo info = InspectMesh(mesh);

	  MFEM_VERIFY(info.relation != AxisRelation::NegativeRadius,
			   "Axisymmetric mesh extends to a negative radius (min r = "
			   + std::to_string(info.min_r)
			   + "). The r coordinate is the mesh x coordinate and must be "
				 "non-negative; check the mesh orientation or units.");

	  return info;
}

} // namespace axisym
