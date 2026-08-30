// Copyright (c) 2026 T. C. Raymond
// SPDX-License-Identifier: MIT
//
// Writes results in Gmsh MSH 2.2 ASCII format. The single output file contains
// one mesh block followed by zero or more $NodeData / $ElementNodeData views.
// View names are exposed via Gmsh's StringTags[0] and form the contract with
// downstream consumers (e.g. TfmrLib's FEMSolution loader).
//
// Format reference: https://gmsh.info/doc/texinfo/gmsh.html#MSH-file-format-version-2-_0028Legacy_0029

#pragma once

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
#include <stdexcept>
#include <string>
#include <system_error>
#include <vector>

#include "mfem.hpp"

namespace gmsh_results {

/// Description of one Gmsh view to emit alongside the mesh block.
struct View {
    enum class Kind { NodeData, ElementNodeData };

    std::string name;        ///< View name, written as StringTags[0].
    Kind        kind;
    int         num_components; ///< 1 (scalar) or 3 (vector, padded in 2D).

    /// NodeData: out has num_components entries; called once per mesh node.
    std::function<void(int node_id, double* out)> node_eval;

    /// ElementNodeData: called once per (element, local node) with the
    /// integration point in reference coordinates; out has num_components.
    std::function<void(int elem_id,
                       const mfem::IntegrationPoint& ip,
                       mfem::ElementTransformation& T,
                       double* out)> elem_node_eval;
};

namespace detail {

inline int GmshElementType(mfem::Geometry::Type geom) {
    using G = mfem::Geometry;
    switch (geom) {
        case G::SEGMENT:     return 1;
        case G::TRIANGLE:    return 2;
        case G::SQUARE:      return 3;
        case G::TETRAHEDRON: return 4;
        case G::CUBE:        return 5;
        case G::PRISM:       return 6;
        default:
            throw std::runtime_error(
                "gmsh_results: unsupported element geometry for MSH 2.2 export");
    }
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

inline void WriteMeshBlock(std::ostream& out, mfem::Mesh& mesh) {
    out << "$MeshFormat\n2.2 0 8\n$EndMeshFormat\n";

    const int nv   = mesh.GetNV();
    const int ne   = mesh.GetNE();
    const int sdim = mesh.SpaceDimension();

    // Stage each section into a single std::string and flush with one
    // out.write(). Avoids per-token ostream overhead and lets the OS see
    // bulk I/O.
    std::string blk;
    blk.reserve(static_cast<size_t>(nv) * 48 + 64);

    blk.assign("$Nodes\n");
    AppendInt(blk, nv);
    blk.push_back('\n');
    for (int i = 0; i < nv; ++i) {
        const double* v = mesh.GetVertex(i);
        const double x = v[0];
        const double y = (sdim > 1) ? v[1] : 0.0;
        const double z = (sdim > 2) ? v[2] : 0.0;
        AppendInt(blk, i + 1);
        blk.push_back(' '); AppendDouble(blk, x);
        blk.push_back(' '); AppendDouble(blk, y);
        blk.push_back(' '); AppendDouble(blk, z);
        blk.push_back('\n');
    }
    blk.append("$EndNodes\n");
    out.write(blk.data(), static_cast<std::streamsize>(blk.size()));

    blk.clear();
    blk.reserve(static_cast<size_t>(ne) * 32 + 64);
    blk.append("$Elements\n");
    AppendInt(blk, ne);
    blk.push_back('\n');
    mfem::Array<int> verts;
    for (int e = 0; e < ne; ++e) {
        const auto geom = mesh.GetElementBaseGeometry(e);
        const int  type = GmshElementType(geom);
        const int  attr = mesh.GetAttribute(e);
        mesh.GetElementVertices(e, verts);

        AppendInt(blk, e + 1);
        blk.push_back(' '); AppendInt(blk, type);
        blk.append(" 2 "); AppendInt(blk, attr);
        blk.push_back(' '); AppendInt(blk, attr);
        for (int k = 0; k < verts.Size(); ++k) {
            blk.push_back(' ');
            AppendInt(blk, verts[k] + 1);
        }
        blk.push_back('\n');
    }
    blk.append("$EndElements\n");
    out.write(blk.data(), static_cast<std::streamsize>(blk.size()));
}

inline void WriteViewHeader(std::ostream& out,
                            const char*   tag,
                            const std::string& name,
                            int num_components,
                            int num_entities) {
    std::string s;
    s.reserve(64 + name.size());
    s.push_back('$'); s.append(tag); s.push_back('\n');
    s.append("1\n\""); s.append(name); s.append("\"\n");
    s.append("1\n0.0\n");
    s.append("3\n0\n");
    AppendInt(s, num_components); s.push_back('\n');
    AppendInt(s, num_entities);   s.push_back('\n');
    out.write(s.data(), static_cast<std::streamsize>(s.size()));
}

inline void WriteNodeData(std::ostream& out, mfem::Mesh& mesh, const View& v) {
    const int nv = mesh.GetNV();
    WriteViewHeader(out, "NodeData", v.name, v.num_components, nv);

    std::vector<double> buf(v.num_components);
    std::string blk;
    blk.reserve(static_cast<size_t>(nv) * (16 + v.num_components * 18));
    for (int i = 0; i < nv; ++i) {
        v.node_eval(i, buf.data());
        AppendInt(blk, i + 1);
        for (int c = 0; c < v.num_components; ++c) {
            blk.push_back(' ');
            AppendDouble(blk, buf[c]);
        }
        blk.push_back('\n');
    }
    blk.append("$EndNodeData\n");
    out.write(blk.data(), static_cast<std::streamsize>(blk.size()));
}

inline void WriteElementNodeData(std::ostream& out,
                                 mfem::Mesh& mesh,
                                 mfem::FiniteElementSpace& sample_fes,
                                 const View& v) {
    const int ne = mesh.GetNE();
    WriteViewHeader(out, "ElementNodeData", v.name, v.num_components, ne);

    std::vector<double> buf(v.num_components);
    std::string blk;
    // Rough upper bound; grows on demand if needed.
    blk.reserve(static_cast<size_t>(ne) * (16 + 4 * v.num_components * 18));
    for (int e = 0; e < ne; ++e) {
        const mfem::FiniteElement* fe = sample_fes.GetFE(e);
        const mfem::IntegrationRule& ir = fe->GetNodes();
        mfem::ElementTransformation* T = mesh.GetElementTransformation(e);

        const int n_local = ir.GetNPoints();
        AppendInt(blk, e + 1);
        blk.push_back(' ');
        AppendInt(blk, n_local);
        for (int k = 0; k < n_local; ++k) {
            const mfem::IntegrationPoint& ip = ir.IntPoint(k);
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

/// Convenience: NodeData view from a scalar H1 GridFunction whose underlying
/// FE space lives on @p mesh and uses linear nodal H1 (one DOF per vertex).
inline View MakeScalarNodeView(const std::string& name,
                               mfem::GridFunction& gf) {
    View v;
    v.name = name;
    v.kind = View::Kind::NodeData;
    v.num_components = 1;
    v.node_eval = [&gf](int node_id, double* out) {
        // For an order-1 H1 space the i-th vertex DOF equals gf[i].
        out[0] = gf(node_id);
    };
    return v;
}

/// Convenience: ElementNodeData view of a vector GridFunction (typically L2,
/// vdim = SpaceDimension). Output is always padded to 3 components.
inline View MakeVectorElementNodeView(const std::string& name,
                                      mfem::GridFunction& vec_gf) {
    View v;
    v.name = name;
    v.kind = View::Kind::ElementNodeData;
    v.num_components = 3;
    v.elem_node_eval = [&vec_gf](int elem_id,
                                 const mfem::IntegrationPoint& ip,
                                 mfem::ElementTransformation& T,
                                 double* out) {
        mfem::Vector val;
        vec_gf.GetVectorValue(T, ip, val);
        out[0] = val.Size() > 0 ? val(0) : 0.0;
        out[1] = val.Size() > 1 ? val(1) : 0.0;
        out[2] = val.Size() > 2 ? val(2) : 0.0;
    };
    return v;
}

/// Convenience: ElementNodeData scalar view sampling a scalar GridFunction.
/// Used for derived scalar coefficients projected onto the export mesh.
inline View MakeScalarElementNodeView(const std::string& name,
                                      mfem::GridFunction& scalar_gf) {
    View v;
    v.name = name;
    v.kind = View::Kind::ElementNodeData;
    v.num_components = 1;
    v.elem_node_eval = [&scalar_gf](int elem_id,
                                    const mfem::IntegrationPoint& ip,
                                    mfem::ElementTransformation& T,
                                    double* out) {
        out[0] = scalar_gf.GetValue(T, ip);
    };
    return v;
}

/// Convenience: ElementNodeData scalar view containing |vec_gf|.
inline View MakeMagnitudeElementNodeView(const std::string& name,
                                         mfem::GridFunction& vec_gf) {
    View v;
    v.name = name;
    v.kind = View::Kind::ElementNodeData;
    v.num_components = 1;
    v.elem_node_eval = [&vec_gf](int elem_id,
                                 const mfem::IntegrationPoint& ip,
                                 mfem::ElementTransformation& T,
                                 double* out) {
        mfem::Vector val;
        vec_gf.GetVectorValue(T, ip, val);
        double s = 0.0;
        for (int i = 0; i < val.Size(); ++i) s += val(i) * val(i);
        out[0] = std::sqrt(s);
    };
    return v;
}

/// Writes the mesh and all views to @p path in MSH 2.2 ASCII.
///
/// @param path          Destination file (will be overwritten).
/// @param mesh          Mesh to embed (typically a refined export copy).
/// @param sample_fes    FE space on @p mesh whose per-element nodal
///                      IntegrationRule defines the local-node layout used by
///                      every ElementNodeData view. Pass the L2 vector space
///                      backing the "E" view.
/// @param views         Views to emit, in order.
inline void WriteGmshResults(const std::string& path,
                             mfem::Mesh& mesh,
                             mfem::FiniteElementSpace& sample_fes,
                             const std::vector<View>& views) {
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

    detail::WriteMeshBlock(out, mesh);

    for (const auto& v : views) {
        switch (v.kind) {
            case View::Kind::NodeData:
                detail::WriteNodeData(out, mesh, v);
                break;
            case View::Kind::ElementNodeData:
                detail::WriteElementNodeData(out, mesh, sample_fes, v);
                break;
        }
    }
}

} // namespace gmsh_results
