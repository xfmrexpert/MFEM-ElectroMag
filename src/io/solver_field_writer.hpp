// Copyright (c) 2026 T. C. Raymond
// SPDX-License-Identifier: MIT

#pragma once
#include <algorithm>
#include <filesystem>
#include <list>
#include <optional>
#include <string>
#include <vector>
#include "mfem.hpp"
#include "field_export.hpp"
#include "gmsh_results_writer.hpp"
#include "status_reporter.hpp"

/**
 * @brief Serializes a FieldExportSet to the supported result formats.
 *
 * A dumb sink: it computes nothing and invents no fields. Solvers decide WHAT
 * to export by building a FieldExportSet; this class only decides HOW to write
 * it. Keeping it out of the solver hierarchy means the serialization details
 * (export meshes, projection spaces, file naming) are not mixed in with the
 * physics.
 *
 * @warning The mesh reference must outlive this writer instance.
 */
class SolverFieldWriter {
public:
	// Gmsh export strategy. Linear is the resample-to-linear path used today and
	// by every current consumer; HighOrder is reserved for true high-order MSH
	// output (native Lagrange elements + $InterpolationScheme) and is not wired
	// up yet.
	enum class GmshStrategy { Linear, HighOrder };

	// @p solution_order is the order of the H1 space the primary fields live in;
	// it drives both the ParaView Lagrange cell order and the default Gmsh
	// export refinement.
	SolverFieldWriter(mfem::Mesh& mesh,
					  std::string results_directory,
					  std::string mesh_path,
					  int solution_order,
					  std::optional<int> export_refine)
		: mesh(mesh),
		  results_directory(std::move(results_directory)),
		  mesh_path(std::move(mesh_path)),
		  solution_order(solution_order),
		  export_refine(export_refine) {}

	// ParaView serializer. Primary scalars are registered at native (high) order
	// so ParaView's Lagrange cells render them faithfully; derived coefficients
	// are projected into an L2 space one order below the H1 solution.
	void WriteParaview(const std::string& collection_name,
					   const FieldExportSet& fields) const
	{
		const int dim = mesh.Dimension();
		MFEM_ASSERT(dim == 2, "ParaView export assumes a 2D mesh.");

		mfem::ParaViewDataCollection pv(collection_name, &mesh);

		// SetHighOrderOutput emits Lagrange cells of the solution order rather
		// than subdividing into linear cells.
		pv.SetLevelsOfDetail(solution_order);
		pv.SetHighOrderOutput(solution_order > 1);

		// Derived fields live in L2 spaces one order below the H1 solution.
		// ParaView stores raw pointers, so the backing spaces and grid functions
		// must outlive Save(); a std::list never relocates its elements, so the
		// registered addresses stay valid as we append.
		const int l2_order = std::max(0, solution_order - 1);
		mfem::L2_FECollection fec_l2(l2_order, dim);
		mfem::FiniteElementSpace fes_l2_scalar(&mesh, &fec_l2);
		mfem::FiniteElementSpace fes_l2_vec(&mesh, &fec_l2, dim);
		std::list<mfem::GridFunction> derived;

		for (const auto& f : fields.Fields()) {
			switch (f.kind) {
				case FieldExport::Kind::Primary:
					pv.RegisterField(f.name, f.primary);
					break;
				case FieldExport::Kind::DerivedScalar:
					derived.emplace_back(&fes_l2_scalar);
					derived.back().ProjectCoefficient(*f.scalar);
					pv.RegisterField(f.name, &derived.back());
					break;
				case FieldExport::Kind::DerivedVector:
					derived.emplace_back(&fes_l2_vec);
					derived.back() = 0.0;
					derived.back().ProjectCoefficient(*f.vector);
					pv.RegisterField(f.name, &derived.back());
					break;
			}
		}

		pv.SetCycle(0);
		pv.SetTime(0.0);
		pv.Save();

		Reporter().Diagnostic("Wrote ParaView collection " + collection_name);
	}

	// Gmsh serializer. The Linear strategy resamples every field onto a refined
	// export mesh as linear views (the standard tessellation approach, faithful
	// at ref_factor ~ order).
	void WriteGmsh(const std::string& scenario_name,
				   const FieldExportSet& fields,
				   GmshStrategy strategy = GmshStrategy::Linear) const
	{
		switch (strategy) {
			case GmshStrategy::Linear:
				WriteGmshLinear(scenario_name, fields);
				return;
			case GmshStrategy::HighOrder:
				WriteGmshHighOrder(scenario_name, fields);
				return;
		}
	}

private:
	mfem::Mesh& mesh;
	std::string results_directory;
	std::string mesh_path;
	int solution_order;
	std::optional<int> export_refine;

	StatusReporter& Reporter() const {
		return StatusReporter::Global();
	}

	// Resample-to-linear Gmsh path: a refined export mesh with linear H1 (node)
	// and L2 (element-node) views. Primary scalars are sampled into H1 nodal
	// views; derived scalars/vectors are projected onto the export mesh.
	void WriteGmshLinear(const std::string& scenario_name,
						 const FieldExportSet& fields) const
	{
		namespace fs = std::filesystem;
		const fs::path results_dir = results_directory.empty()
			? fs::path(mesh_path).parent_path()
			: fs::path(results_directory);
		const fs::path out_path = results_dir / (scenario_name + ".results.msh");

		const int ref_factor = export_refine.value_or(std::max(1, solution_order));
		mfem::Mesh export_mesh(&mesh, ref_factor, mfem::BasisType::ClosedUniform);

		const int dim = export_mesh.Dimension();
		const int sdim = export_mesh.SpaceDimension();
		MFEM_ASSERT(dim == 2, "Gmsh export assumes a 2D mesh.");

		// Linear export spaces, shared by every view of each kind.
		mfem::H1_FECollection fec_h1_lin(1, dim);
		mfem::FiniteElementSpace fes_scalar(&export_mesh, &fec_h1_lin);
		mfem::L2_FECollection fec_l2_lin(1, dim);
		mfem::FiniteElementSpace fes_vec(&export_mesh, &fec_l2_lin, sdim);

		// Resampled fields must outlive WriteGmshResults; std::list keeps element
		// addresses stable as views capture them by reference.
		std::list<mfem::GridFunction> resampled;
		std::vector<gmsh_results::View> views;

		for (const auto& f : fields.Fields()) {
			switch (f.kind) {
				case FieldExport::Kind::Primary: {
					// This path resamples through a scalar coefficient, so it only
					// handles scalar primaries. A vector-valued primary needs a
					// vector view instead; fail loudly rather than silently
					// exporting component 0.
					MFEM_VERIFY(f.primary->VectorDim() == 1,
						"Gmsh export of vector-valued primary field '" + f.name +
						"' is not implemented.");
					resampled.emplace_back(&fes_scalar);
					mfem::GridFunctionCoefficient sc(f.primary);
					resampled.back().ProjectCoefficient(sc);
					views.push_back(
						gmsh_results::MakeScalarNodeView(f.name, resampled.back()));
					break;
				}
				case FieldExport::Kind::DerivedScalar: {
					resampled.emplace_back(&fes_scalar);
					resampled.back().ProjectCoefficient(*f.scalar);
					views.push_back(
						gmsh_results::MakeScalarElementNodeView(f.name, resampled.back()));
					break;
				}
				case FieldExport::Kind::DerivedVector: {
					resampled.emplace_back(&fes_vec);
					resampled.back() = 0.0;
					resampled.back().ProjectCoefficient(*f.vector);
					views.push_back(
						gmsh_results::MakeVectorElementNodeView(f.name, resampled.back()));
					break;
				}
			}
		}

		gmsh_results::WriteGmshResults(out_path.string(), export_mesh, fes_vec, views);
		Reporter().Diagnostic("Wrote " + out_path.string());
	}

	// True high-order MSH output (native high-order elements + per-view
	// interpolation scheme). Deferred: needs an MFEM -> Gmsh high-order node
	// ordering map. Stubbed so the strategy switch is complete and the work is
	// a single localized follow-up.
	void WriteGmshHighOrder(const std::string& scenario_name,
							const FieldExportSet& /*fields*/) const
	{
		Reporter().Diagnostic("WriteGmsh(HighOrder) not implemented yet for scenario: "
							  + scenario_name + "; falling back to no output.");
	}
};
