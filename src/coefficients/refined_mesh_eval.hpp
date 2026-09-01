// Copyright (c) 2026 T. C. Raymond
// SPDX-License-Identifier: MIT

/**
 * @file refined_mesh_eval.hpp
 * @brief Evaluate a GridFunction through an ElementTransformation that may
 *        belong to a refined copy of the GridFunction's own mesh.
 */

#pragma once

#include "mfem.hpp"

namespace refined_eval {

/**
 * @brief Map an integration point on a refined mesh back to its parent element.
 *
 * Results export builds a uniformly refined copy of the solution mesh
 * (mfem::Mesh(&mesh, ref_factor, ...)) and samples fields through element
 * transformations belonging to that finer mesh. A GridFunction defined on the
 * original mesh cannot be evaluated with such a transformation directly:
 * GridFunction::GetValue indexes its own FE space with T.ElementNo, which
 * overruns once the refined mesh has more elements than the coarse one.
 *
 * MFEM's own GridFunctionCoefficient / VectorGridFunctionCoefficient guard this
 * with an internal RefinedToCoarse helper, but that helper is file-local to
 * fem/coefficient.cpp and is not declared in any public header, so custom
 * coefficients must reproduce it.
 *
 * @param      coarse_mesh Mesh the GridFunction is defined on.
 * @param      T           Transformation, possibly on a refined mesh.
 * @param      ip          Integration point in @a T's reference coordinates.
 * @param[out] mapped_ip   Receives the equivalent point in the parent element's
 *                         reference coordinates, or a copy of @a ip when no
 *                         mapping is needed.
 * @return Transformation to evaluate with: the parent coarse element's
 *         transformation, or @a T itself when @a T already refers to
 *         @a coarse_mesh.
 *
 * @warning The returned transformation is the mesh's shared scratch object, as
 *          in MFEM; it is invalidated by the next GetElementTransformation call
 *          on the same mesh.
 */
inline mfem::ElementTransformation& Resolve(mfem::Mesh& coarse_mesh,
											mfem::ElementTransformation& T,
											const mfem::IntegrationPoint& ip,
											mfem::IntegrationPoint& mapped_ip) {
   if (T.mesh == nullptr || T.mesh->GetNE() == coarse_mesh.GetNE()) {
	  mapped_ip = ip;
	  T.SetIntPoint(&mapped_ip);
	  return T;
   }

   const mfem::Mesh& fine_mesh = *T.mesh;
   const mfem::CoarseFineTransformations& cf =
	  fine_mesh.GetRefinementTransforms();
   const int fine_element = T.ElementNo;
   const int coarse_element = cf.embeddings[fine_element].parent;

   const mfem::Geometry::Type geom = T.GetGeometryType();
   mfem::IntegrationPointTransformation fine_to_coarse;
   mfem::IsoparametricTransformation& emb_tr = fine_to_coarse.Transf;
   emb_tr.SetIdentityTransformation(geom);
   emb_tr.SetPointMat(
	  cf.point_matrices[geom](cf.embeddings[fine_element].matrix));
   fine_to_coarse.Transform(ip, mapped_ip);

   mfem::ElementTransformation* coarse_T =
	  coarse_mesh.GetElementTransformation(coarse_element);
   coarse_T->SetIntPoint(&mapped_ip);
   return *coarse_T;
}

/// Convenience overload resolving the coarse mesh from @a gf itself.
inline mfem::ElementTransformation& Resolve(const mfem::GridFunction& gf,
											mfem::ElementTransformation& T,
											const mfem::IntegrationPoint& ip,
											mfem::IntegrationPoint& mapped_ip) {
   return Resolve(*gf.FESpace()->GetMesh(), T, ip, mapped_ip);
}

} // namespace refined_eval
