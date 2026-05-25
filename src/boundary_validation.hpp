// Copyright (c) 2026 T. C. Raymond
// SPDX-License-Identifier: MIT

#pragma once
#include "mfem.hpp"
#include <vector>
#include <map>
#include <set>
#include <sstream>
#include <tuple>
#include <utility>

/// Validates that boundary conditions don't create conflicts
/// Throws an exception if the same DOF is assigned multiple different values
class BoundaryConditionValidator {
private:
    mfem::Mesh& mesh;
    mfem::FiniteElementSpace& fespace;

public:
    BoundaryConditionValidator(mfem::Mesh& m, mfem::FiniteElementSpace& fes)
        : mesh(m), fespace(fes) {}

    /// Check if boundary conditions create conflicts
    /// @param bcs Vector of (boundary_marker, value) pairs
    /// @param allow_conflicts If true, warn but don't throw on conflicts
    /// @throws std::runtime_error if conflicts detected and !allow_conflicts
    void ValidateBoundaryConditions(
        const std::vector<std::pair<mfem::Array<int>, double>>& bcs,
        bool allow_conflicts = false)
    {
        // Map from DOF index to (value, boundary_attribute, bdr_element_index)
        std::map<int, std::vector<std::tuple<double, int, int>>> dof_values;

        // For each boundary condition
        for (const auto& [marker, value] : bcs) {
            // Process each boundary element
            for (int i = 0; i < fespace.GetNBE(); i++) {
                int attr = mesh.GetBdrAttribute(i);

                // Skip if this BC doesn't apply to this boundary
                if (attr < 1 || attr > marker.Size() || !marker[attr-1]) {
                    continue;
                }

                // Get DOFs for this boundary element
                mfem::Array<int> vdofs;
                fespace.GetBdrElementVDofs(i, vdofs);

                // Record value for each DOF (track which boundary element raised it so we
                // can report geometric coordinates if a conflict is later detected).
                for (int j = 0; j < vdofs.Size(); j++) {
                    int dof = vdofs[j];
                    if (dof < 0) dof = -1 - dof;  // Handle orientation
                    dof_values[dof].push_back(std::make_tuple(value, attr, i));
                }
            }
        }

        // Check for conflicts
        std::vector<std::string> conflicts;
        const double tolerance = 1e-10;

        for (const auto& [dof, values] : dof_values) {
            if (values.size() > 1) {
                // Check if all values are the same
                double first_val = std::get<0>(values[0]);
                bool has_conflict = false;

                for (size_t i = 1; i < values.size(); i++) {
                    if (std::abs(std::get<0>(values[i]) - first_val) > tolerance) {
                        has_conflict = true;
                        break;
                    }
                }

                if (has_conflict) {
                    std::ostringstream oss;
                    oss << "DOF " << dof
                        << " is constrained by multiple boundaries with DIFFERENT values:\n";

                    // Locate the conflicting DOF. In an H1 FE space the first NV DOFs are
                    // vertex DOFs; the next NE * (order-1) DOFs are edge-interior DOFs;
                    // then face/element DOFs. We classify by index range and print the
                    // appropriate geometric location.
                    const int NV = mesh.GetNV();
                    const int NE = mesh.GetNEdges();
                    oss << "    (NV=" << NV << ", NEdges=" << NE
                        << ", total NDofs=" << fespace.GetNDofs() << ")\n";

                    if (dof < NV) {
                        const double* v = mesh.GetVertex(dof);
                        oss << "    Vertex DOF -> vertex " << dof << " at ("
                            << v[0] << ", " << v[1];
                        if (mesh.SpaceDimension() > 2) oss << ", " << v[2];
                        oss << ")\n";
                    } else if (dof < NV + NE) {
                        int edge_idx = dof - NV;
                        mfem::Array<int> ev;
                        mesh.GetEdgeVertices(edge_idx, ev);
                        const double* v0 = mesh.GetVertex(ev[0]);
                        const double* v1 = mesh.GetVertex(ev[1]);
                        double mx = 0.5 * (v0[0] + v1[0]);
                        double my = 0.5 * (v0[1] + v1[1]);
                        oss << "    Edge-interior DOF -> mesh edge " << edge_idx
                            << " between v" << ev[0] << "=(" << v0[0] << ", " << v0[1] << ")"
                            << " and v" << ev[1] << "=(" << v1[0] << ", " << v1[1] << ")"
                            << "  midpoint=(" << mx << ", " << my << ")\n";
                    } else {
                        // Face / interior DOF -- fall back to the first touching bdr element.
                        int be = std::get<2>(values[0]);
                        mfem::Array<int> verts;
                        mesh.GetBdrElementVertices(be, verts);
                        double cx = 0, cy = 0;
                        for (int vi = 0; vi < verts.Size(); ++vi) {
                            const double* v = mesh.GetVertex(verts[vi]);
                            cx += v[0]; cy += v[1];
                        }
                        cx /= verts.Size(); cy /= verts.Size();
                        oss << "    Face/interior DOF, bdr-element centroid: ("
                            << cx << ", " << cy << ")\n";
                    }

                    // Group entries by (attr, value) so we don't flood the output when
                    // many boundary elements of the same attribute touch a shared DOF.
                    std::map<std::pair<int, double>, int> grouped;
                    for (const auto& t : values) {
                        grouped[{std::get<1>(t), std::get<0>(t)}]++;
                    }
                    for (const auto& [key, count] : grouped) {
                        oss << "    Boundary attribute " << key.first
                            << ": value = " << key.second
                            << "  (" << count << " bdr element"
                            << (count == 1 ? "" : "s") << ")\n";
                    }

                    // Dump the offending boundary elements with FULL raw vdofs so we
                    // can see whether DOF X arrived as a vertex DOF (index < NV) or as
                    // an edge DOF (index >= NV) and which mesh edge produced it.
                    oss << "    Offending boundary elements (attribute  vertices  raw vdofs):\n";
                    std::set<int> printed_be;
                    int dumped = 0;
                    for (const auto& t : values) {
                        int be = std::get<2>(t);
                        if (!printed_be.insert(be).second) continue;
                        mfem::Array<int> verts;
                        mesh.GetBdrElementVertices(be, verts);
                        mfem::Array<int> bvd;
                        fespace.GetBdrElementVDofs(be, bvd);
                        oss << "      attr " << std::get<1>(t) << "  bdr#" << be;
                        for (int vi = 0; vi < verts.Size(); ++vi) {
                            const double* v = mesh.GetVertex(verts[vi]);
                            oss << "  v" << verts[vi]
                                << "=(" << v[0] << ", " << v[1] << ")";
                        }
                        oss << "  vdofs=[";
                        for (int k = 0; k < bvd.Size(); ++k) {
                            if (k) oss << ", ";
                            oss << bvd[k];
                        }
                        oss << "]\n";
                        if (++dumped >= 8) {
                            oss << "      ... (more bdr elements omitted)\n";
                            break;
                        }
                    }
                    conflicts.push_back(oss.str());
                }
            }
        }

        if (!conflicts.empty()) {
            std::ostringstream error_msg;
            error_msg << "\n========================================\n";
            if (allow_conflicts) {
                error_msg << "WARNING: Conflicting Boundary Conditions Detected!\n";
            } else {
                error_msg << "ERROR: Conflicting Boundary Conditions Detected!\n";
            }
            error_msg << "========================================\n\n";
            error_msg << "The same DOF is constrained to multiple different values.\n";
            error_msg << "This typically indicates:\n";
            error_msg << "  1. Mesh problem: Boundaries with different BCs share nodes\n";
            error_msg << "  2. Physics problem: Trying to set a point to two voltages simultaneously\n\n";
            error_msg << "Conflicts found (showing first 5):\n\n";

            int count = 0;
            for (const auto& conflict : conflicts) {
                error_msg << conflict << "\n";
                if (++count >= 5) break;
            }
            if (conflicts.size() > 5) {
                error_msg << "... and " << (conflicts.size() - 5) << " more conflicts\n";
            }

            if (allow_conflicts) {
                error_msg << "\n⚠️  PROCEEDING WITH CONFLICTS (values will be averaged)\n";
                error_msg << "This may produce INCORRECT PHYSICS!\n";
                error_msg << "\nRecommended fixes:\n";
                error_msg << "  - Regenerate mesh with non-overlapping boundaries\n";
                error_msg << "  - Check boundary condition specification\n";
                error_msg << "========================================\n";
                mfem::out << error_msg.str() << std::endl;
            } else {
                error_msg << "\nPlease fix your:\n";
                error_msg << "  - Mesh design (separate boundaries shouldn't share nodes)\n";
                error_msg << "  - Boundary condition specification (check your config.json)\n";
                error_msg << "========================================\n";
                throw std::runtime_error(error_msg.str());
            }
        }
    }
};
