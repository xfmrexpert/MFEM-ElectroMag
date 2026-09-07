// Copyright (c) 2026 T. C. Raymond
// SPDX-License-Identifier: MIT
//
// Covers both directions of the Gmsh MSH format contract:
//
//  INPUT  - MFEM 4.10 added MSH 4.1 reading alongside the long-standing 2.2
//           support. The mesher that produces case.msh emits 2.2 today, so the
//           guarantee we care about is that BOTH versions load and yield the
//           same element attributes, since attributes are what the config's
//           entity_groups/attribute_ids bind to. If 4.1 reading ever silently
//           dropped or renumbered attributes, every region and boundary
//           condition would bind to the wrong elements.
//
//  OUTPUT - gmsh_results_writer emits 2.2 by default and 4.1 on request. The
//           downstream C# consumer reads 2.2, so the default must not drift.
//           For 4.1 the attribute moves out of the element line and into
//           $Entities as a physical tag; these tests pin that mapping.

#include <catch2/catch_test_macros.hpp>

#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "mfem.hpp"
#include "io/gmsh_results_writer.hpp"

namespace {

// Two triangles forming a unit square, split across two physical surfaces
// (attributes 7 and 12, mirroring the winding/conductor ids used in the
// examples). The same mesh is expressed below in MSH 2.2 and MSH 4.1.
const char* kMsh22 = R"MSHFILE(
$MeshFormat
2.2 0 8
$EndMeshFormat
$Nodes
4
1 0 0 0
2 1 0 0
3 1 1 0
4 0 1 0
$EndNodes
$Elements
2
1 2 2 7 7 1 2 3
2 2 2 12 12 1 3 4
$EndElements
)MSHFILE";

const char* kMsh41 = R"MSHFILE(
$MeshFormat
4.1 0 8
$EndMeshFormat
$Entities
0 0 2 0
7 0 0 0 1 1 0 1 7 0
12 0 0 0 1 1 0 1 12 0
$EndEntities
$Nodes
1 4 1 4
2 7 0 4
1
2
3
4
0 0 0
1 0 0
1 1 0
0 1 0
$EndNodes
$Elements
2 2 1 2
2 7 2 1
1 1 2 3
2 12 2 1
2 1 3 4
$EndElements
)MSHFILE";

// MFEM's Gmsh reader works off a stream, so no temp files are needed.
mfem::Mesh LoadMsh(const char* text) {
	// Skip the leading newline the raw literal introduces.
	std::string s(text);
	if (!s.empty() && s.front() == '\n') { s.erase(0, 1); }
	std::istringstream in(s);
	return mfem::Mesh(in, /*generate_edges=*/0, /*refine=*/0);
}

std::vector<int> Attributes(mfem::Mesh& mesh) {
	std::vector<int> attrs;
	for (int e = 0; e < mesh.GetNE(); ++e) {
		attrs.push_back(mesh.GetAttribute(e));
	}
	return attrs;
}

std::string WriteToString(mfem::Mesh& mesh, int order,
						  gmsh_results::MshVersion version) {
	const auto path = std::filesystem::temp_directory_path()
					  / (version == gmsh_results::MshVersion::V4_1
							 ? "mfem_em_fmt_41.msh"
							 : "mfem_em_fmt_22.msh");
	gmsh_results::WriteGmshResults(path.string(), mesh, order, {}, version);

	std::ifstream in(path, std::ios::binary);
	std::ostringstream ss;
	ss << in.rdbuf();
	return ss.str();
}

bool Contains(const std::string& haystack, const std::string& needle) {
	return haystack.find(needle) != std::string::npos;
}

} // namespace

TEST_CASE("MSH 2.2 and 4.1 input meshes load identically", "[gmsh][msh][input]") {
	mfem::Mesh m22 = LoadMsh(kMsh22);
	mfem::Mesh m41 = LoadMsh(kMsh41);

	REQUIRE(m22.GetNE() == 2);
	REQUIRE(m41.GetNE() == m22.GetNE());
	REQUIRE(m41.GetNV() == m22.GetNV());
	REQUIRE(m41.Dimension() == m22.Dimension());

	// The contract that entity_groups/attribute_ids depend on: physical group
	// ids survive as MFEM element attributes, identically in both formats.
	REQUIRE(Attributes(m41) == Attributes(m22));
	REQUIRE(Attributes(m22) == std::vector<int>{ 7, 12 });
}

TEST_CASE("MSH 4.1 input preserves per-attribute element counts", "[gmsh][msh][input]") {
	mfem::Mesh m41 = LoadMsh(kMsh41);

	int n7 = 0, n12 = 0;
	for (int e = 0; e < m41.GetNE(); ++e) {
		if (m41.GetAttribute(e) == 7) { ++n7; }
		if (m41.GetAttribute(e) == 12) { ++n12; }
	}
	REQUIRE(n7 == 1);
	REQUIRE(n12 == 1);
}

TEST_CASE("Results writer defaults to MSH 2.2", "[gmsh][msh][output]") {
	mfem::Mesh mesh = LoadMsh(kMsh22);

	// Called without a version argument, exactly as existing callers do.
	const auto path = std::filesystem::temp_directory_path()
					  / "mfem_em_fmt_default.msh";
	gmsh_results::WriteGmshResults(path.string(), mesh, 1, {});

	std::ifstream in(path, std::ios::binary);
	std::ostringstream ss;
	ss << in.rdbuf();
	const std::string out = ss.str();

	// The downstream C# consumer reads 2.2; this default must not drift.
	REQUIRE(Contains(out, "$MeshFormat\n2.2 0 8\n$EndMeshFormat\n"));
	REQUIRE_FALSE(Contains(out, "$Entities"));
}

TEST_CASE("MSH 4.1 output carries attributes as entity physical tags", "[gmsh][msh][output]") {
	mfem::Mesh mesh = LoadMsh(kMsh22);
	const std::string out =
		WriteToString(mesh, 1, gmsh_results::MshVersion::V4_1);

	REQUIRE(Contains(out, "$MeshFormat\n4.1 0 8\n$EndMeshFormat\n"));

	// One surface entity per distinct attribute, each carrying exactly one
	// physical tag equal to the attribute.
	REQUIRE(Contains(out, "$Entities\n0 0 2 0\n"));
	REQUIRE(Contains(out, " 1 7 0\n"));
	REQUIRE(Contains(out, " 1 12 0\n"));

	// Element blocks are keyed by (entityTag == attribute, elementType).
	// Type 2 is the 3-node triangle, one element per block here.
	REQUIRE(Contains(out, "2 7 2 1\n"));
	REQUIRE(Contains(out, "2 12 2 1\n"));
}

TEST_CASE("MSH 4.1 output round-trips through the MFEM reader", "[gmsh][msh][output]") {
	// The strongest available check that the emitted 4.1 is well formed:
	// read it back and confirm the mesh and its attributes are unchanged.
	mfem::Mesh original = LoadMsh(kMsh22);
	const std::string out =
		WriteToString(original, 1, gmsh_results::MshVersion::V4_1);

	std::istringstream in(out);
	mfem::Mesh reloaded(in, /*generate_edges=*/0, /*refine=*/0);

	REQUIRE(reloaded.GetNE() == original.GetNE());
	REQUIRE(reloaded.GetNV() == original.GetNV());
	REQUIRE(Attributes(reloaded) == Attributes(original));
}
