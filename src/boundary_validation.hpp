// Copyright (c) 2026 T. C. Raymond
// SPDX-License-Identifier: MIT

#pragma once
#include "mfem.hpp"
#include <vector>
#include <map>
#include <sstream>
#include <tuple>
#include <utility>
#include <unordered_map>
#include <string>
#include <algorithm>

/// Validates that boundary constraints (closures + terminals) are well-posed.
///
/// Two independent checks are performed:
///   1. Closure value conflicts - a DOF pinned by closure conditions (far-field,
///      symmetry, axis) to two DIFFERENT fixed values is ill-posed.
///   2. Terminal ownership/overlap - terminal values are scenario-dependent, so a
///      DOF may be owned by at most ONE terminal and must not also be a closure
///      DOF. Any shared DOF is guaranteed to conflict in some scenario.
class BoundaryConditionValidator {
private:
    mfem::Mesh& mesh;
    mfem::FiniteElementSpace& fespace;

    /// Collect the orientation-normalized DOFs touched by a boundary marker,
    /// remembering one touching boundary element per DOF for diagnostics.
    void CollectMarkerDofs(const mfem::Array<int>& marker,
                           std::map<int, int>& dof_to_be) const {
        for (int i = 0; i < fespace.GetNBE(); ++i) {
            int attr = mesh.GetBdrAttribute(i);
            if (attr < 1 || attr > marker.Size() || !marker[attr - 1]) continue;
            mfem::Array<int> vdofs;
            fespace.GetBdrElementVDofs(i, vdofs);
            for (int j = 0; j < vdofs.Size(); ++j) {
                int dof = vdofs[j];
                if (dof < 0) dof = -1 - dof;  // handle orientation
                dof_to_be.emplace(dof, i);    // first touching bdr element wins
            }
        }
    }

    /// Human-readable geometric description of where a DOF lives. In an H1 space
    /// the first NV DOFs are vertex DOFs, the next NEdges*(order-1) are edge
    /// interiors, then face/interior DOFs (located via a touching bdr element).
    std::string DescribeDof(int dof, int fallback_be) const {
        std::ostringstream oss;
        const int NV = mesh.GetNV();
        const int NE = mesh.GetNEdges();
        if (dof < NV) {
            const double* v = mesh.GetVertex(dof);
            oss << "vertex " << dof << " at (" << v[0] << ", " << v[1];
            if (mesh.SpaceDimension() > 2) oss << ", " << v[2];
            oss << ")";
        } else if (dof < NV + NE) {
            int edge_idx = dof - NV;
            mfem::Array<int> ev;
            mesh.GetEdgeVertices(edge_idx, ev);
            const double* v0 = mesh.GetVertex(ev[0]);
            const double* v1 = mesh.GetVertex(ev[1]);
            oss << "edge-interior on mesh edge " << edge_idx
                << " midpoint=(" << 0.5 * (v0[0] + v1[0]) << ", "
                << 0.5 * (v0[1] + v1[1]) << ")";
        } else if (fallback_be >= 0) {
            mfem::Array<int> verts;
            mesh.GetBdrElementVertices(fallback_be, verts);
            double cx = 0, cy = 0;
            for (int vi = 0; vi < verts.Size(); ++vi) {
                const double* v = mesh.GetVertex(verts[vi]);
                cx += v[0]; cy += v[1];
            }
            cx /= verts.Size(); cy /= verts.Size();
            oss << "face/interior near bdr-element centroid ("
                << cx << ", " << cy << ")";
        } else {
            oss << "index " << dof;
        }
        return oss.str();
    }

public:
    BoundaryConditionValidator(mfem::Mesh& m, mfem::FiniteElementSpace& fes)
        : mesh(m), fespace(fes) {}

    /// @param closures      (boundary_marker, value) pairs for closure conditions.
    /// @param terminals     terminal name -> boundary marker for voltage terminals.
    /// @param allow_overlap If true, warn but don't throw on detected problems.
    /// @throws std::runtime_error if problems are detected and !allow_overlap.
    void ValidateBoundaryConditions(
        const std::vector<std::pair<mfem::Array<int>, double>>& closures,
        const std::unordered_map<std::string, mfem::Array<int>>& terminals,
        bool allow_overlap = false)
    {
        std::vector<std::string> problems;
        const double tolerance = 1e-10;

        // (1) Closure value conflicts: same DOF pinned to different fixed values.
        std::map<int, std::vector<std::tuple<double, int, int>>> dof_values;
        for (const auto& [marker, value] : closures) {
            for (int i = 0; i < fespace.GetNBE(); ++i) {
                int attr = mesh.GetBdrAttribute(i);
                if (attr < 1 || attr > marker.Size() || !marker[attr - 1]) continue;
                mfem::Array<int> vdofs;
                fespace.GetBdrElementVDofs(i, vdofs);
                for (int j = 0; j < vdofs.Size(); ++j) {
                    int dof = vdofs[j];
                    if (dof < 0) dof = -1 - dof;  // handle orientation
                    dof_values[dof].push_back(std::make_tuple(value, attr, i));
                }
            }
        }

        for (const auto& [dof, values] : dof_values) {
            if (values.size() < 2) continue;
            double first_val = std::get<0>(values[0]);
            bool conflict = false;
            for (size_t i = 1; i < values.size(); ++i) {
                if (std::abs(std::get<0>(values[i]) - first_val) > tolerance) {
                    conflict = true;
                    break;
                }
            }
            if (!conflict) continue;

            std::ostringstream oss;
            oss << "DOF " << dof << " (" << DescribeDof(dof, std::get<2>(values[0]))
                << ") is pinned by closures to DIFFERENT values:\n";
            std::map<std::pair<int, double>, int> grouped;
            for (const auto& t : values) grouped[{std::get<1>(t), std::get<0>(t)}]++;
            for (const auto& [key, count] : grouped) {
                oss << "    boundary attribute " << key.first
                    << ": value = " << key.second
                    << "  (" << count << " bdr element" << (count == 1 ? "" : "s") << ")\n";
            }
            problems.push_back(oss.str());
        }

        // (2) Terminal ownership/overlap: each DOF owned by at most one terminal,
        //     and never shared with a closure. Terminal values vary per scenario,
        //     so any shared DOF is guaranteed to conflict in some scenario.
        std::map<int, int> closure_dofs;  // dof -> touching bdr element
        for (const auto& [marker, value] : closures) {
            (void)value;
            CollectMarkerDofs(marker, closure_dofs);
        }

        std::vector<std::string> term_names;
        term_names.reserve(terminals.size());
        for (const auto& kv : terminals) term_names.push_back(kv.first);
        std::sort(term_names.begin(), term_names.end());  // deterministic output

        std::map<int, std::vector<std::string>> dof_terminals;  // dof -> terminal names
        std::map<int, int> terminal_dof_be;                     // dof -> touching bdr element
        for (const auto& name : term_names) {
            std::map<int, int> td;
            CollectMarkerDofs(terminals.at(name), td);
            for (const auto& [dof, be] : td) {
                dof_terminals[dof].push_back(name);
                terminal_dof_be.emplace(dof, be);
            }
        }

        for (const auto& [dof, names] : dof_terminals) {
            const bool multi_terminal = names.size() > 1;
            const bool hits_closure   = closure_dofs.count(dof) > 0;
            if (!multi_terminal && !hits_closure) continue;

            int be = terminal_dof_be.count(dof) ? terminal_dof_be.at(dof) : -1;
            std::ostringstream oss;
            oss << "DOF " << dof << " (" << DescribeDof(dof, be)
                << ") is constrained by overlapping essential sites:\n";
            if (multi_terminal) {
                oss << "    terminals: ";
                for (size_t i = 0; i < names.size(); ++i) {
                    if (i) oss << ", ";
                    oss << names[i];
                }
                oss << "\n";
            } else {
                oss << "    terminal: " << names[0] << "\n";
            }
            if (hits_closure) {
                oss << "    closure boundary condition (also pins this DOF)\n";
            }
            problems.push_back(oss.str());
        }

        if (problems.empty()) return;

        std::ostringstream msg;
        msg << "\n========================================\n";
        msg << (allow_overlap ? "WARNING" : "ERROR")
            << ": Conflicting Boundary Constraints Detected!\n";
        msg << "========================================\n\n";
        msg << "The same DOF is constrained by incompatible essential conditions.\n";
        msg << "This typically indicates:\n";
        msg << "  1. Mesh: boundaries/terminals with different roles share nodes\n";
        msg << "  2. Config: a node is driven to two values at once (now or per-scenario)\n\n";
        msg << "Problems found (showing first 5):\n\n";

        int count = 0;
        for (const auto& p : problems) {
            msg << p << "\n";
            if (++count >= 5) break;
        }
        if (problems.size() > 5) {
            msg << "... and " << (problems.size() - 5) << " more\n";
        }

        if (allow_overlap) {
            msg << "\nPROCEEDING WITH OVERLAPS (values will be averaged at shared DOFs).\n";
            msg << "This may produce INCORRECT PHYSICS!\n";
            msg << "\nRecommended fixes:\n";
            msg << "  - Regenerate the mesh so distinct boundaries/terminals don't share nodes\n";
            msg << "  - Check terminal attribute_ids and boundary specification\n";
            msg << "========================================\n";
            mfem::out << msg.str() << std::endl;
        } else {
            msg << "\nPlease fix your:\n";
            msg << "  - Mesh design (separate boundaries/terminals shouldn't share nodes)\n";
            msg << "  - Terminal/boundary specification (check your config.json)\n";
            msg << "========================================\n";
            throw std::runtime_error(msg.str());
        }
    }
};
