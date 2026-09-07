// Copyright (c) 2026 T. C. Raymond
// SPDX-License-Identifier: MIT
//
// Writes results in Gmsh MSH ASCII format, version 2.2 (default) or 4.1. The
// single output file contains one mesh block followed by zero or more
// $NodeData / $ElementNodeData views. View names are exposed via Gmsh's
// StringTags[0] and form the contract with downstream consumers (e.g.
// TfmrLib's FEMSolution loader).
//
// Elements are emitted as native Gmsh Lagrange elements of the solution order
// (types 2/9/21/23... for triangles, 3/10/36/37... for quads), so Gmsh
// interpolates the field with the matching high-order shape functions. No
// refined/tessellated export copy of the mesh is made.
//
// An $InterpolationScheme block carries the shape functions themselves (a
// monomial exponent matrix plus Lagrange coefficients) and every view names it
// via StringTags[1]. A consumer can therefore evaluate any field at an
// arbitrary point without hardcoding Gmsh node ordering or Lagrange formulas,
// and without changes when the solution order changes.
//
// FIELD REPRESENTATION POLICY: fields are exported exactly as the FE solution
// represents them. Continuous primaries (V, A) go out as $NodeData; derived
// quantities that are genuinely discontinuous across elements (E = -grad V,
// B = curl A, and anything built from them) go out as $ElementNodeData with
// per-element values. Inter-element jumps are real results of the
// discretization and are passed through unsmoothed; how to treat them (average,
// recover, or respect) is the consumer's decision, since it depends on the
// post-processing being done.
//
// FORMAT SEAM: only the MESH sections differ between MSH 2.2 and 4.1, and all
// of that is confined to WriteMeshFormat, WriteMeshBlock22, and
// WriteMeshBlock41 behind the WriteMeshBlock dispatcher. Gmsh's
// post-processing sections are not versioned along with the mesh format, so
// $InterpolationScheme, $NodeData, and $ElementNodeData are emitted
// byte-identically for both versions. The node layout (HoLayout/GetHoLayout),
// the MFEM->Gmsh permutation, and all field sampling are likewise
// format-agnostic.
//
// The one semantic difference worth knowing: 2.2 stores each element's
// attribute on the element line, whereas 4.1 stores it once per model entity
// in $Entities and groups elements under those entities. Both paths therefore
// round-trip the same MFEM attributes.
//
// Format references:
//   2.2: https://gmsh.info/doc/texinfo/gmsh.html#MSH-file-format-version-2-_0028Legacy_0029
//   4.1: https://gmsh.info/doc/texinfo/gmsh.html#MSH-file-format

#pragma once

#include <algorithm>
#include <array>
#include <cerrno>
#include <charconv>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <functional>
#include <ios>
#include <map>
#include <stdexcept>
#include <string>
#include <system_error>
#include <utility>
#include <vector>

#include "mfem.hpp"

namespace gmsh_results {

/// MSH file format version of the emitted mesh sections.
///
/// V2_2 is the default and the format the existing downstream C# consumer
/// (TfmrLib's FEMSolution loader) understands. V4_1 is opt-in.
enum class MshVersion { V2_2, V4_1 };

/// Parse a format string ("2.2" / "4.1") as used in the simulation config.
/// Throws std::runtime_error on an unrecognized value.
inline MshVersion ParseMshVersion(const std::string& text) {
    if (text == "2.2") { return MshVersion::V2_2; }
    if (text == "4.1") { return MshVersion::V4_1; }
    throw std::runtime_error(
        "gmsh_results: unsupported MSH format '" + text
        + "'; expected \"2.2\" or \"4.1\"");
}

/// Description of one Gmsh view to emit alongside the mesh block.
struct View {
    enum class Kind { NodeData, ElementNodeData };

    std::string name;        ///< View name, written as StringTags[0].
    Kind        kind;
    int         num_components; ///< 1 (scalar) or 3 (vector, padded in 2D).

    /// Called once per sample point with the point in reference coordinates;
    /// out has num_components entries. Both view kinds sample by evaluation,
    /// so continuous (NodeData) and discontinuous (ElementNodeData) fields
    /// share one callback shape.
    std::function<void(int elem_id,
                       const mfem::IntegrationPoint& ip,
                       mfem::ElementTransformation& T,
                       double* out)> elem_node_eval;
};

namespace detail {

// Reference-space node layout of a Gmsh Lagrange element of a given order.
//
// Node ORDER is Gmsh's own, defined recursively: all corner nodes, then the
// nodes interior to each edge (edges traversed in element-local order), then
// the nodes interior to the face, which are themselves laid out as a lower
// order element of the same shape. Coordinates are in MFEM's reference domain
// (unit triangle / unit square) so they feed straight into
// ElementTransformation::Transform.
//
// Working in reference COORDINATES rather than a permutation of MFEM DOF
// indices keeps this independent of MFEM's internal nodal layout and of which
// basis (Gauss-Lobatto vs equispaced) the solution happens to use: every
// exported value is obtained by evaluating at a point.
struct HoLayout {
    int gmsh_type = 0;
    std::vector<std::array<double, 2>> ref;
};

using RefPt = std::array<double, 2>;

inline RefPt Lerp(const RefPt& a, const RefPt& b, double t) {
    return { a[0] + (b[0] - a[0]) * t, a[1] + (b[1] - a[1]) * t };
}

// Gmsh element type codes, indexed by order. Index 0 is unused.
inline int TriangleGmshType(int order) {
    static const int kTypes[] = { 0, 2, 9, 21, 23, 25, 42, 43, 44, 45, 46 };
    return (order >= 1 && order <= 10) ? kTypes[order] : 0;
}

inline int QuadGmshType(int order) {
    static const int kTypes[] = { 0, 3, 10, 36, 37, 38, 47, 48, 49, 50, 51 };
    return (order >= 1 && order <= 10) ? kTypes[order] : 0;
}

// Emit an order-p triangular lattice over the triangle (a, b, c) in Gmsh order.
// p == 0 degenerates to the single centroid node, which is how the recursion
// terminates for orders that leave exactly one interior node.
inline void AppendTriangleNodes(int p, const RefPt& a, const RefPt& b,
                                const RefPt& c, std::vector<RefPt>& out) {
    if (p == 0) {
        out.push_back({ (a[0] + b[0] + c[0]) / 3.0, (a[1] + b[1] + c[1]) / 3.0 });
        return;
    }
    out.push_back(a);
    out.push_back(b);
    out.push_back(c);

    const double inv = 1.0 / static_cast<double>(p);
    const RefPt* edges[3][2] = { { &a, &b }, { &b, &c }, { &c, &a } };
    for (auto& e : edges) {
        for (int i = 1; i < p; ++i) {
            out.push_back(Lerp(*e[0], *e[1], i * inv));
        }
    }

    // Interior nodes form a triangle of order p-3, inset by one lattice step.
    if (p < 3) return;
    const RefPt ai = { a[0] + (b[0] - a[0] + c[0] - a[0]) * inv,
                       a[1] + (b[1] - a[1] + c[1] - a[1]) * inv };
    const RefPt bi = { b[0] + (a[0] - b[0] + c[0] - b[0]) * inv,
                       b[1] + (a[1] - b[1] + c[1] - b[1]) * inv };
    const RefPt ci = { c[0] + (a[0] - c[0] + b[0] - c[0]) * inv,
                       c[1] + (a[1] - c[1] + b[1] - c[1]) * inv };
    AppendTriangleNodes(p - 3, ai, bi, ci, out);
}

// Emit an order-p quadrilateral lattice over (a, b, c, d) in Gmsh order.
inline void AppendQuadNodes(int p, const RefPt& a, const RefPt& b,
                            const RefPt& c, const RefPt& d,
                            std::vector<RefPt>& out) {
    if (p == 0) {
        out.push_back({ (a[0] + b[0] + c[0] + d[0]) * 0.25,
                        (a[1] + b[1] + c[1] + d[1]) * 0.25 });
        return;
    }
    out.push_back(a);
    out.push_back(b);
    out.push_back(c);
    out.push_back(d);

    const double inv = 1.0 / static_cast<double>(p);
    const RefPt* edges[4][2] = { { &a, &b }, { &b, &c }, { &c, &d }, { &d, &a } };
    for (auto& e : edges) {
        for (int i = 1; i < p; ++i) {
            out.push_back(Lerp(*e[0], *e[1], i * inv));
        }
    }

    // Interior nodes form a quad of order p-2, inset by one lattice step.
    if (p < 2) return;
    const RefPt ai = { a[0] + (b[0] - a[0] + d[0] - a[0]) * inv,
                       a[1] + (b[1] - a[1] + d[1] - a[1]) * inv };
    const RefPt bi = { b[0] + (a[0] - b[0] + c[0] - b[0]) * inv,
                       b[1] + (a[1] - b[1] + c[1] - b[1]) * inv };
    const RefPt ci = { c[0] + (d[0] - c[0] + b[0] - c[0]) * inv,
                       c[1] + (d[1] - c[1] + b[1] - c[1]) * inv };
    const RefPt di = { d[0] + (c[0] - d[0] + a[0] - d[0]) * inv,
                       d[1] + (c[1] - d[1] + a[1] - d[1]) * inv };
    AppendQuadNodes(p - 2, ai, bi, ci, di, out);
}

inline const char* GeometryName(mfem::Geometry::Type geom) {
    switch (geom) {
        case mfem::Geometry::SEGMENT:     return "segment";
        case mfem::Geometry::TRIANGLE:    return "triangle";
        case mfem::Geometry::SQUARE:      return "quadrilateral";
        case mfem::Geometry::TETRAHEDRON: return "tetrahedron";
        case mfem::Geometry::CUBE:        return "hexahedron";
        case mfem::Geometry::PRISM:       return "prism";
        default:                          return "unknown";
    }
}

// Layout for one (geometry, order) pair. Built once and cached: the recursion
// is cheap but this is called per element.
inline const HoLayout& GetHoLayout(mfem::Geometry::Type geom, int order) {
    static std::map<std::pair<int, int>, HoLayout> cache;
    const auto key = std::make_pair(static_cast<int>(geom), order);
    auto it = cache.find(key);
    if (it != cache.end()) return it->second;

    if (order < 1 || order > 10) {
        throw std::runtime_error(
            "gmsh_results: unsupported export order " + std::to_string(order)
            + " for a " + GeometryName(geom)
            + "; Gmsh Lagrange elements are supported for orders 1-10");
    }

    HoLayout layout;
    switch (geom) {
        case mfem::Geometry::TRIANGLE:
            layout.gmsh_type = TriangleGmshType(order);
            AppendTriangleNodes(order, { 0.0, 0.0 }, { 1.0, 0.0 }, { 0.0, 1.0 },
                                layout.ref);
            break;
        case mfem::Geometry::SQUARE:
            layout.gmsh_type = QuadGmshType(order);
            AppendQuadNodes(order, { 0.0, 0.0 }, { 1.0, 0.0 }, { 1.0, 1.0 },
                            { 0.0, 1.0 }, layout.ref);
            break;
        default:
            throw std::runtime_error(
                std::string("gmsh_results: unsupported element geometry '")
                + GeometryName(geom)
                + "' for MSH export; only triangles and quadrilaterals are "
                  "handled (the export path assumes a 2D mesh)");
    }
    return cache.emplace(key, std::move(layout)).first->second;
}

// Hot-path numeric formatting helpers.
//
// std::ostream::operator<<(double) is ~5-10x slower than std::to_chars on
// MSVC because it queries the stream locale and routes through num_put. For
// MSH export we write millions of doubles, so format into a small stack
// buffer and write raw bytes through filebuf instead.
//
// 9 significant decimal digits is well past single-precision round-trip
// (max_digits10 = 9) and matches what most visualization consumers parse;
// drop the previous setprecision(16) which doubled the per-value byte count
// for no visible benefit.
inline void AppendInt(std::string& s, long long v) {
    char buf[24];
    auto r = std::to_chars(buf, buf + sizeof(buf), v);
    s.append(buf, r.ptr);
}

inline void AppendDouble(std::string& s, double v) {
    char buf[32];
    auto r = std::to_chars(buf, buf + sizeof(buf), v,
                           std::chars_format::scientific, 9);
    s.append(buf, r.ptr);
}

// Per-element map from Gmsh node slot -> MFEM local DOF index.
//
// Derived by matching reference coordinates rather than by hardcoding MFEM's
// nodal ordering, so it stays correct if that ordering ever changes. All
// elements of the same geometry and order share one permutation.
inline std::vector<int> BuildDofPermutation(const mfem::FiniteElement& fe,
                                            const HoLayout& layout) {
    const mfem::IntegrationRule& ir = fe.GetNodes();
    const int n = ir.GetNPoints();
    if (n != static_cast<int>(layout.ref.size())) {
        throw std::runtime_error(
            "gmsh_results: node count mismatch between MFEM element ("
            + std::to_string(n) + ") and Gmsh layout ("
            + std::to_string(layout.ref.size()) + ")");
    }

    std::vector<int> perm(n, -1);
    for (int k = 0; k < n; ++k) {
        for (int j = 0; j < n; ++j) {
            const mfem::IntegrationPoint& ip = ir.IntPoint(j);
            if (std::fabs(ip.x - layout.ref[k][0]) < 1e-10 &&
                std::fabs(ip.y - layout.ref[k][1]) < 1e-10) {
                perm[k] = j;
                break;
            }
        }
        if (perm[k] < 0) {
            throw std::runtime_error(
                "gmsh_results: no MFEM node matches Gmsh node slot "
                + std::to_string(k)
                + "; the export space must use equispaced (ClosedUniform) nodes");
        }
    }
    return perm;
}

// Monomial exponents and Lagrange coefficients defining an element's shape
// functions, in the form Gmsh's $InterpolationScheme block expects.
//
// Gmsh's model is: given exponent matrix E (n_terms x n_vars) and coefficient
// matrix C (n_nodes x n_terms), shape function i is
//
//     phi_i(u, v) = sum_j C[i][j] * u^E[j][0] * v^E[j][1]
//
// and a field is reconstructed as sum_i phi_i(u, v) * value_i, where value_i
// are the nodal values listed for that element (in the same node order as the
// mesh block's connectivity).
//
// Shipping this in the file is what makes the consumer independent of Gmsh's
// node ordering AND of the element order: a reader evaluates the monomials and
// does a matrix-vector product without knowing anything about Lagrange bases.
// This matters most for DISCONTINUOUS ($ElementNodeData) fields such as E and
// B, where nodal values alone cannot be interpolated -- the consumer must
// evaluate the element-local basis to get a value anywhere but a node.
struct InterpScheme {
    std::vector<std::array<int, 2>>  exponents;  // per monomial term
    std::vector<std::vector<double>> coeffs;     // [node][term]
};

// Monomial exponents spanning the polynomial space of an order-p element.
// Triangles use the total-degree space (u^a v^b, a + b <= p); quads use the
// tensor-product space (a <= p, b <= p). These match the node lattices built
// by AppendTriangleNodes / AppendQuadNodes, so the Vandermonde system below
// is square and non-singular.
inline std::vector<std::array<int, 2>> MonomialExponents(
    mfem::Geometry::Type geom, int order) {
    std::vector<std::array<int, 2>> e;
    if (geom == mfem::Geometry::TRIANGLE) {
        for (int d = 0; d <= order; ++d) {
            for (int i = 0; i <= d; ++i) { e.push_back({ d - i, i }); }
        }
    } else {  // SQUARE
        for (int a = 0; a <= order; ++a) {
            for (int b = 0; b <= order; ++b) { e.push_back({ a, b }); }
        }
    }
    return e;
}

// Build the Lagrange coefficient matrix by inverting the Vandermonde system.
//
// Requiring phi_i(node_k) = delta_ik gives V * C^T = I, where
// V[k][j] = monomial_j(node_k). So C^T = V^-1, i.e. C = (V^-1)^T. Solved with
// Gauss-Jordan and partial pivoting; the equispaced lattices used here are
// well enough conditioned at the supported orders (<= 10).
inline InterpScheme BuildInterpScheme(mfem::Geometry::Type geom, int order) {
    InterpScheme scheme;
    scheme.exponents = MonomialExponents(geom, order);
    const HoLayout& layout = GetHoLayout(geom, order);

    const int n = static_cast<int>(layout.ref.size());
    const int m = static_cast<int>(scheme.exponents.size());
    if (n != m) {
        throw std::runtime_error(
            "gmsh_results: interpolation basis size (" + std::to_string(m)
            + ") does not match node count (" + std::to_string(n)
            + ") for a " + GeometryName(geom) + " of order "
            + std::to_string(order));
    }

    // Augmented [V | I]; Gauss-Jordan leaves V^-1 in the right half.
    std::vector<std::vector<double>> a(n, std::vector<double>(2 * n, 0.0));
    for (int k = 0; k < n; ++k) {
        const double u = layout.ref[k][0];
        const double v = layout.ref[k][1];
        for (int j = 0; j < n; ++j) {
            a[k][j] = std::pow(u, scheme.exponents[j][0]) *
                      std::pow(v, scheme.exponents[j][1]);
        }
        a[k][n + k] = 1.0;
    }

    for (int col = 0; col < n; ++col) {
        int piv = col;
        for (int r = col + 1; r < n; ++r) {
            if (std::fabs(a[r][col]) > std::fabs(a[piv][col])) { piv = r; }
        }
        if (std::fabs(a[piv][col]) < 1e-12) {
            throw std::runtime_error(
                "gmsh_results: singular Vandermonde matrix building the "
                "interpolation scheme for a " + std::string(GeometryName(geom))
                + " of order " + std::to_string(order));
        }
        std::swap(a[col], a[piv]);

        const double inv_p = 1.0 / a[col][col];
        for (int j = 0; j < 2 * n; ++j) { a[col][j] *= inv_p; }
        for (int r = 0; r < n; ++r) {
            if (r == col) { continue; }
            const double f = a[r][col];
            if (f == 0.0) { continue; }
            for (int j = 0; j < 2 * n; ++j) { a[r][j] -= f * a[col][j]; }
        }
    }

    // C = (V^-1)^T: row i holds the monomial coefficients of phi_i.
    scheme.coeffs.assign(n, std::vector<double>(n, 0.0));
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) { scheme.coeffs[i][j] = a[j][n + i]; }
    }
    return scheme;
}

// Node numbering and geometry for the exported high-order mesh.
//
// Node ids are the DOF indices of an equispaced H1 space, which are shared
// across element boundaries, so the emitted mesh is continuous. Each node also
// records one (element, reference point) pair so field values can be sampled
// there later without recomputing the layout.
struct ExportNodes {
    int order = 1;
    std::vector<std::array<double, 3>> coord;  // indexed by node id
    std::vector<int>   node_elem;              // representative element
    std::vector<RefPt> node_ref;               // reference point in that element
    std::vector<int>   elem_type;              // Gmsh type code, per element
    std::vector<std::vector<int>> elem_nodes;  // Gmsh-ordered node ids, per element
};

inline ExportNodes BuildExportNodes(mfem::Mesh& mesh,
                                    mfem::FiniteElementSpace& fes,
                                    int order) {
    ExportNodes nodes;
    nodes.order = order;

    const int nd = fes.GetNDofs();
    const int ne = mesh.GetNE();
    nodes.coord.assign(nd, { 0.0, 0.0, 0.0 });
    nodes.node_elem.assign(nd, -1);
    nodes.node_ref.assign(nd, RefPt{ 0.0, 0.0 });
    nodes.elem_type.assign(ne, 0);
    nodes.elem_nodes.assign(ne, {});

    std::map<int, std::vector<int>> perm_cache;  // keyed by geometry
    mfem::Array<int> dofs;
    mfem::Vector phys;

    for (int e = 0; e < ne; ++e) {
        const auto geom = mesh.GetElementBaseGeometry(e);
        const HoLayout& layout = GetHoLayout(geom, order);
        nodes.elem_type[e] = layout.gmsh_type;

        auto pit = perm_cache.find(static_cast<int>(geom));
        if (pit == perm_cache.end()) {
            pit = perm_cache.emplace(
                static_cast<int>(geom),
                BuildDofPermutation(*fes.GetFE(e), layout)).first;
        }
        const std::vector<int>& perm = pit->second;

        fes.GetElementDofs(e, dofs);
        mfem::ElementTransformation* T = mesh.GetElementTransformation(e);

        const int n = static_cast<int>(layout.ref.size());
        nodes.elem_nodes[e].resize(n);
        for (int k = 0; k < n; ++k) {
            const int dof = dofs[perm[k]];
            nodes.elem_nodes[e][k] = dof;

            mfem::IntegrationPoint ip;
            ip.Set2(layout.ref[k][0], layout.ref[k][1]);

            // Transform() gives the true geometric position for curved
            // elements (mesh.GetNodes() populated) and reduces to the expected
            // edge midpoints/lattice for straight-sided ones, so no separate
            // curved-vs-straight handling is needed.
            T->SetIntPoint(&ip);
            T->Transform(ip, phys);

            nodes.coord[dof] = { phys(0),
                                 phys.Size() > 1 ? phys(1) : 0.0,
                                 phys.Size() > 2 ? phys(2) : 0.0 };
            nodes.node_elem[dof] = e;
            nodes.node_ref[dof] = layout.ref[k];
        }
    }
    return nodes;
}

inline void WriteMeshFormat(std::ostream& out, MshVersion version) {
    // Fields are: version, file-type (0 = ASCII), data-size (sizeof(double)).
    out << (version == MshVersion::V4_1
                ? "$MeshFormat\n4.1 0 8\n$EndMeshFormat\n"
                : "$MeshFormat\n2.2 0 8\n$EndMeshFormat\n");
}

inline void WriteMeshBlock22(std::ostream& out, mfem::Mesh& mesh,
                             const ExportNodes& nodes) {
    const int nv = static_cast<int>(nodes.coord.size());
    const int ne = mesh.GetNE();

    // Stage each section into a single std::string and flush with one
    // out.write(). Avoids per-token ostream overhead and lets the OS see
    // bulk I/O.
    std::string blk;
    blk.reserve(static_cast<size_t>(nv) * 48 + 64);

    blk.assign("$Nodes\n");
    AppendInt(blk, nv);
    blk.push_back('\n');
    for (int i = 0; i < nv; ++i) {
        AppendInt(blk, i + 1);
        blk.push_back(' '); AppendDouble(blk, nodes.coord[i][0]);
        blk.push_back(' '); AppendDouble(blk, nodes.coord[i][1]);
        blk.push_back(' '); AppendDouble(blk, nodes.coord[i][2]);
        blk.push_back('\n');
    }
    blk.append("$EndNodes\n");
    out.write(blk.data(), static_cast<std::streamsize>(blk.size()));

    blk.clear();
    blk.reserve(static_cast<size_t>(ne) * 32 + 64);
    blk.append("$Elements\n");
    AppendInt(blk, ne);
    blk.push_back('\n');
    for (int e = 0; e < ne; ++e) {
        const int attr = mesh.GetAttribute(e);
        AppendInt(blk, e + 1);
        blk.push_back(' '); AppendInt(blk, nodes.elem_type[e]);
        blk.append(" 2 "); AppendInt(blk, attr);
        blk.push_back(' '); AppendInt(blk, attr);
        for (int id : nodes.elem_nodes[e]) {
            blk.push_back(' ');
            AppendInt(blk, id + 1);
        }
        blk.push_back('\n');

        if (blk.size() > (1u << 20)) {
            out.write(blk.data(), static_cast<std::streamsize>(blk.size()));
            blk.clear();
        }
    }
    blk.append("$EndElements\n");
    out.write(blk.data(), static_cast<std::streamsize>(blk.size()));
}

// MSH 4.1 mesh sections.
//
// The structural difference from 2.2 is where an element's attribute lives.
// In 2.2 every element line carries its own tags ("2 <attr> <attr>"). In 4.1
// elements carry no tags at all; they are grouped into blocks keyed by the
// model entity they are classified on, and a separate $Entities section maps
// each entity to its physical tags. To preserve exactly the attribute
// semantics the 2.2 path exposes, we synthesize one surface entity per
// distinct element attribute, with entityTag == physicalTag == attribute.
//
// Node tags remain the global 1..nv numbering used by the 2.2 path, so the
// $NodeData / $ElementNodeData sections that follow are byte-identical between
// the two formats.
inline void WriteMeshBlock41(std::ostream& out, mfem::Mesh& mesh,
                             const ExportNodes& nodes) {
    const int nv = static_cast<int>(nodes.coord.size());
    const int ne = mesh.GetNE();

    // Group elements by (attribute, gmsh element type). A single attribute may
    // legitimately contain both triangles and quads, and MSH 4.1 requires one
    // block per element type, so the pair is the block key. std::map keeps the
    // emission order deterministic.
    std::map<std::pair<int, int>, std::vector<int>> blocks;
    for (int e = 0; e < ne; ++e) {
        blocks[{ mesh.GetAttribute(e), nodes.elem_type[e] }].push_back(e);
    }

    // Distinct attributes become the surface entities.
    std::vector<int> attrs;
    for (const auto& kv : blocks) {
        if (std::find(attrs.begin(), attrs.end(), kv.first.first) == attrs.end()) {
            attrs.push_back(kv.first.first);
        }
    }
    std::sort(attrs.begin(), attrs.end());

    std::string s;

    // $Entities: numPoints numCurves numSurfaces numVolumes, then one line per
    // surface: tag, bounding box, physical tags, bounding curves (none, since
    // we do not synthesize a curve topology).
    s.append("$Entities\n0 0 ");
    AppendInt(s, static_cast<long long>(attrs.size()));
    s.append(" 0\n");
    for (int attr : attrs) {
        // Bounding box over the export nodes of every element with this
        // attribute. Gmsh tolerates a loose box; it is used for display only.
        double lo[3] = { 1e300, 1e300, 1e300 };
        double hi[3] = { -1e300, -1e300, -1e300 };
        for (const auto& kv : blocks) {
            if (kv.first.first != attr) { continue; }
            for (int e : kv.second) {
                for (int id : nodes.elem_nodes[e]) {
                    for (int c = 0; c < 3; ++c) {
                        lo[c] = std::min(lo[c], nodes.coord[id][c]);
                        hi[c] = std::max(hi[c], nodes.coord[id][c]);
                    }
                }
            }
        }
        AppendInt(s, attr);
        for (int c = 0; c < 3; ++c) { s.push_back(' '); AppendDouble(s, lo[c]); }
        for (int c = 0; c < 3; ++c) { s.push_back(' '); AppendDouble(s, hi[c]); }
        s.append(" 1 ");        // one physical tag...
        AppendInt(s, attr);     // ...which is the attribute itself
        s.append(" 0\n");       // zero bounding curves
    }
    s.append("$EndEntities\n");
    out.write(s.data(), static_cast<std::streamsize>(s.size()));

    // $Nodes: numEntityBlocks numNodes minNodeTag maxNodeTag, then per block
    // entityDim entityTag parametric numNodesInBlock, all tags, then all
    // coordinates. Export nodes are shared across attributes, and MSH 4.1
    // requires each node tag to appear exactly once, so all nodes go in a
    // single block classified on the first surface entity.
    s.clear();
    s.reserve(static_cast<size_t>(nv) * 48 + 128);
    if (attrs.empty() || nv == 0) {
        s.append("$Nodes\n0 0 0 0\n$EndNodes\n");
    } else {
        s.append("$Nodes\n1 ");
        AppendInt(s, nv); s.append(" 1 "); AppendInt(s, nv); s.push_back('\n');
        s.append("2 ");             // entityDim = 2 (surface)
        AppendInt(s, attrs.front()); // entityTag
        s.append(" 0 ");            // parametric = 0
        AppendInt(s, nv); s.push_back('\n');
        for (int i = 0; i < nv; ++i) {
            AppendInt(s, i + 1);
            s.push_back('\n');
        }
        for (int i = 0; i < nv; ++i) {
            AppendDouble(s, nodes.coord[i][0]);
            s.push_back(' '); AppendDouble(s, nodes.coord[i][1]);
            s.push_back(' '); AppendDouble(s, nodes.coord[i][2]);
            s.push_back('\n');
        }
        s.append("$EndNodes\n");
    }
    out.write(s.data(), static_cast<std::streamsize>(s.size()));

    // $Elements: numEntityBlocks numElements minElementTag maxElementTag, then
    // per block entityDim entityTag elementType numElementsInBlock followed by
    // "elementTag nodeTag..." lines. Element tags keep the global 1..ne
    // numbering so they still line up with $ElementNodeData below.
    s.clear();
    s.reserve(static_cast<size_t>(ne) * 32 + 128);
    s.append("$Elements\n");
    AppendInt(s, static_cast<long long>(blocks.size()));
    s.push_back(' '); AppendInt(s, ne);
    s.append(ne > 0 ? " 1 " : " 0 ");
    AppendInt(s, ne);
    s.push_back('\n');
    for (const auto& kv : blocks) {
        s.append("2 ");                          // entityDim = 2 (surface)
        AppendInt(s, kv.first.first);            // entityTag = attribute
        s.push_back(' '); AppendInt(s, kv.first.second);  // elementType
        s.push_back(' ');
        AppendInt(s, static_cast<long long>(kv.second.size()));
        s.push_back('\n');
        for (int e : kv.second) {
            AppendInt(s, e + 1);
            for (int id : nodes.elem_nodes[e]) {
                s.push_back(' ');
                AppendInt(s, id + 1);
            }
            s.push_back('\n');

            if (s.size() > (1u << 20)) {
                out.write(s.data(), static_cast<std::streamsize>(s.size()));
                s.clear();
            }
        }
    }
    s.append("$EndElements\n");
    out.write(s.data(), static_cast<std::streamsize>(s.size()));
}

inline void WriteMeshBlock(std::ostream& out, mfem::Mesh& mesh,
                           const ExportNodes& nodes, MshVersion version) {
    WriteMeshFormat(out, version);
    if (version == MshVersion::V4_1) {
        WriteMeshBlock41(out, mesh, nodes);
    } else {
        WriteMeshBlock22(out, mesh, nodes);
    }
}

// Emit one $InterpolationScheme block describing the shape functions for every
// element type present in the mesh.
//
// Gmsh associates the scheme with views by name; WriteViewHeader writes the
// same name as the view's StringTags[1], which is how a reader (and Gmsh
// itself) links a view's nodal values to the basis that interpolates them.
//
// Each element topology gets two matrices per Gmsh's format: the coefficient
// matrix then the exponent matrix.
inline void WriteInterpolationScheme(std::ostream& out,
                                     const std::string& scheme_name,
                                     mfem::Mesh& mesh,
                                     int order) {
    // Collect the distinct geometries actually used, preserving a stable order.
    std::vector<mfem::Geometry::Type> geoms;
    for (int e = 0; e < mesh.GetNE(); ++e) {
        const auto g = mesh.GetElementBaseGeometry(e);
        if (std::find(geoms.begin(), geoms.end(), g) == geoms.end()) {
            geoms.push_back(g);
        }
    }
    if (geoms.empty()) { return; }

    std::string s;
    s.append("$InterpolationScheme\n");
    s.append("\"" + scheme_name + "\"\n");
    AppendInt(s, static_cast<long long>(geoms.size()));
    s.push_back('\n');

    for (const auto g : geoms) {
        const HoLayout&    layout = GetHoLayout(g, order);
        const InterpScheme scheme = BuildInterpScheme(g, order);
        const int n = static_cast<int>(scheme.coeffs.size());
        const int m = static_cast<int>(scheme.exponents.size());

        AppendInt(s, layout.gmsh_type);
        s.append("\n2\n");  // two matrices follow: coefficients, then exponents

        AppendInt(s, n); s.push_back(' '); AppendInt(s, m); s.push_back('\n');
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < m; ++j) {
                if (j) { s.push_back(' '); }
                AppendDouble(s, scheme.coeffs[i][j]);
            }
            s.push_back('\n');
        }

        AppendInt(s, m); s.append(" 2\n");
        for (int j = 0; j < m; ++j) {
            AppendInt(s, scheme.exponents[j][0]);
            s.push_back(' ');
            AppendInt(s, scheme.exponents[j][1]);
            s.push_back('\n');
        }
    }

    s.append("$EndInterpolationScheme\n");
    out.write(s.data(), static_cast<std::streamsize>(s.size()));
}

inline void WriteViewHeader(std::ostream& out,
                            const char*   tag,
                            const std::string& name,
                            const std::string& scheme_name,
                            int num_components,
                            int num_entities) {
    std::string s;
    s.reserve(96 + name.size() + scheme_name.size());
    s.push_back('$'); s.append(tag); s.push_back('\n');
    // StringTags: [0] view name, [1] interpolation scheme name. The second tag
    // is what binds this view to the $InterpolationScheme block above.
    s.append("2\n\""); s.append(name); s.append("\"\n");
    s.append("\""); s.append(scheme_name); s.append("\"\n");
    s.append("1\n0.0\n");
    s.append("3\n0\n");
    AppendInt(s, num_components); s.push_back('\n');
    AppendInt(s, num_entities);   s.push_back('\n');
    out.write(s.data(), static_cast<std::streamsize>(s.size()));
}

inline void WriteNodeData(std::ostream& out, mfem::Mesh& mesh,
                          const ExportNodes& nodes, const View& v,
                          const std::string& scheme_name) {
    const int nv = static_cast<int>(nodes.coord.size());
    WriteViewHeader(out, "NodeData", v.name, scheme_name, v.num_components, nv);

    std::vector<double> buf(v.num_components);
    std::string blk;
    blk.reserve(static_cast<size_t>(nv) * (16 + v.num_components * 18));
    for (int i = 0; i < nv; ++i) {
        // Sample at the node's recorded reference point rather than reading a
        // DOF directly: the field may live in a different space or basis than
        // the equispaced one used for node numbering.
        mfem::IntegrationPoint ip;
        ip.Set2(nodes.node_ref[i][0], nodes.node_ref[i][1]);
        mfem::ElementTransformation* T =
            mesh.GetElementTransformation(nodes.node_elem[i]);
        T->SetIntPoint(&ip);
        v.elem_node_eval(nodes.node_elem[i], ip, *T, buf.data());

        AppendInt(blk, i + 1);
        for (int c = 0; c < v.num_components; ++c) {
            blk.push_back(' ');
            AppendDouble(blk, buf[c]);
        }
        blk.push_back('\n');

        if (blk.size() > (1u << 20)) {
            out.write(blk.data(), static_cast<std::streamsize>(blk.size()));
            blk.clear();
        }
    }
    blk.append("$EndNodeData\n");
    out.write(blk.data(), static_cast<std::streamsize>(blk.size()));
}

inline void WriteElementNodeData(std::ostream& out,
                                 mfem::Mesh& mesh,
                                 const ExportNodes& nodes,
                                 const View& v,
                                 const std::string& scheme_name) {
    const int ne = mesh.GetNE();
    WriteViewHeader(out, "ElementNodeData", v.name, scheme_name,
                    v.num_components, ne);

    std::vector<double> buf(v.num_components);
    std::string blk;
    // Rough upper bound; grows on demand if needed.
    blk.reserve(static_cast<size_t>(ne) * (16 + 4 * v.num_components * 18));
    for (int e = 0; e < ne; ++e) {
        // Values must be listed in the same node order as the element's
        // connectivity in the mesh block, so drive the loop from the Gmsh
        // layout rather than from an MFEM space's own nodal ordering.
        const HoLayout& layout =
            GetHoLayout(mesh.GetElementBaseGeometry(e), nodes.order);
        mfem::ElementTransformation* T = mesh.GetElementTransformation(e);

        const int n_local = static_cast<int>(layout.ref.size());
        AppendInt(blk, e + 1);
        blk.push_back(' ');
        AppendInt(blk, n_local);
        for (int k = 0; k < n_local; ++k) {
            mfem::IntegrationPoint ip;
            ip.Set2(layout.ref[k][0], layout.ref[k][1]);
            T->SetIntPoint(&ip);
            v.elem_node_eval(e, ip, *T, buf.data());
            for (int c = 0; c < v.num_components; ++c) {
                blk.push_back(' ');
                AppendDouble(blk, buf[c]);
            }
        }
        blk.push_back('\n');

        // Flush periodically so the staging buffer doesn't grow without
        // bound for very large meshes.
        if (blk.size() > (1u << 20)) {
            out.write(blk.data(), static_cast<std::streamsize>(blk.size()));
            blk.clear();
        }
    }
    blk.append("$EndElementNodeData\n");
    out.write(blk.data(), static_cast<std::streamsize>(blk.size()));
}

} // namespace detail

/// Convenience: NodeData view of a continuous scalar GridFunction. The field is
/// sampled by evaluation, so it may be of any order or basis; it need not match
/// the equispaced space used to number the export nodes.
inline View MakeScalarNodeView(const std::string& name,
                               mfem::GridFunction& gf) {
    View v;
    v.name = name;
    v.kind = View::Kind::NodeData;
    v.num_components = 1;
    v.elem_node_eval = [&gf](int /*elem_id*/,
                             const mfem::IntegrationPoint& ip,
                             mfem::ElementTransformation& T,
                             double* out) {
        out[0] = gf.GetValue(T, ip);
    };
    return v;
}

/// Convenience: ElementNodeData view of a vector Coefficient. Output is always
/// padded to 3 components.
///
/// Derived fields are sampled straight from the coefficient at each export
/// node, so there is no intermediate projection onto an L2 space and no
/// projection error; the emitted values are the coefficient's own.
inline View MakeVectorCoefficientView(const std::string& name,
                                      mfem::VectorCoefficient& vec_coeff) {
    View v;
    v.name = name;
    v.kind = View::Kind::ElementNodeData;
    v.num_components = 3;
    v.elem_node_eval = [&vec_coeff](int /*elem_id*/,
                                    const mfem::IntegrationPoint& ip,
                                    mfem::ElementTransformation& T,
                                    double* out) {
        mfem::Vector val(vec_coeff.GetVDim());
        vec_coeff.Eval(val, T, ip);
        out[0] = val.Size() > 0 ? val(0) : 0.0;
        out[1] = val.Size() > 1 ? val(1) : 0.0;
        out[2] = val.Size() > 2 ? val(2) : 0.0;
    };
    return v;
}

/// Convenience: ElementNodeData scalar view sampling a scalar Coefficient.
inline View MakeScalarCoefficientView(const std::string& name,
                                      mfem::Coefficient& coeff) {
    View v;
    v.name = name;
    v.kind = View::Kind::ElementNodeData;
    v.num_components = 1;
    v.elem_node_eval = [&coeff](int /*elem_id*/,
                                const mfem::IntegrationPoint& ip,
                                mfem::ElementTransformation& T,
                                double* out) {
        out[0] = coeff.Eval(T, ip);
    };
    return v;
}

/// Writes the mesh and all views to @p path in Gmsh MSH ASCII, using native
/// Gmsh Lagrange elements of order @p order.
///
/// @param path     Destination file (will be overwritten).
/// @param mesh     Mesh to embed. Exported at its own resolution; no refined
///                 copy is made.
/// @param order    Lagrange order of the emitted elements and of the node layout
///                 every view is sampled on. Typically the solution order.
/// @param views    Views to emit, in order.
/// @param version  MSH format of the mesh sections. Defaults to 2.2, which is
///                 what the downstream C# consumer reads; 4.1 is opt-in.
inline void WriteGmshResults(const std::string& path,
                             mfem::Mesh& mesh,
                             int order,
                             const std::vector<View>& views,
                             MshVersion version = MshVersion::V2_2) {
    // std::ofstream / fopen will not create missing parent directories on
    // Windows or POSIX; opening "./foo/bar.msh" silently fails with ENOENT
    // when ./foo does not yet exist. Create the parent chain up front so a
    // results path like "<mesh-dir>/<case>.results.msh" works on a fresh
    // working tree.
    namespace fs = std::filesystem;
    fs::path out_path(path);
    if (out_path.has_parent_path()) {
        std::error_code ec;
        fs::create_directories(out_path.parent_path(), ec);
        if (ec) {
            throw std::runtime_error(
                "WriteGmshResults: cannot create directory '"
                + out_path.parent_path().string() + "': " + ec.message());
        }
    }

    // Open in binary mode so '\n' is not translated to "\r\n" on Windows;
    // halves write traffic on text-heavy MSH output. Give the filebuf a
    // 1 MiB buffer (default is 8 KiB) so the staged std::string blocks
    // emitted by WriteMeshBlock / WriteNodeData / WriteElementNodeData
    // hit the OS in a handful of large writes instead of hundreds.
    std::ofstream out;
    static thread_local std::array<char, 1u << 20> io_buf;
    out.rdbuf()->pubsetbuf(io_buf.data(),
                           static_cast<std::streamsize>(io_buf.size()));
    out.open(out_path, std::ios::out | std::ios::binary | std::ios::trunc);
    if (!out.is_open()) {
        std::error_code cwd_ec;
        const auto cwd = fs::current_path(cwd_ec);
        throw std::runtime_error(
            "WriteGmshResults: cannot open '" + path + "': "
            + std::strerror(errno)
            + " (cwd=" + (cwd_ec ? std::string("?") : cwd.string()) + ")");
    }
    // Numeric formatting goes through std::to_chars in detail::AppendDouble
    // (locale-independent, ~5-10x faster than operator<<). No stream-side
    // setprecision / scientific needed.

    // Equispaced nodes so the DOF positions coincide with Gmsh's Lagrange
    // node lattice; the resulting DOF indices double as shared node ids.
    mfem::H1_FECollection fec(order, mesh.Dimension(),
                              mfem::BasisType::ClosedUniform);
    mfem::FiniteElementSpace fes(&mesh, &fec);
    const detail::ExportNodes nodes = detail::BuildExportNodes(mesh, fes, order);

    detail::WriteMeshBlock(out, mesh, nodes, version);

    // One scheme shared by every view: all views are sampled on the same node
    // layout at the same order, so they interpolate with the same basis.
    const std::string scheme_name = "MFEM_Lagrange_P" + std::to_string(order);
    detail::WriteInterpolationScheme(out, scheme_name, mesh, order);

    for (const auto& v : views) {
        switch (v.kind) {
            case View::Kind::NodeData:
                detail::WriteNodeData(out, mesh, nodes, v, scheme_name);
                break;
            case View::Kind::ElementNodeData:
                detail::WriteElementNodeData(out, mesh, nodes, v, scheme_name);
                break;
        }
    }
}

} // namespace gmsh_results
