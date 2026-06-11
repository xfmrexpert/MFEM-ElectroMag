// Copyright (c) 2026 T. C. Raymond
// SPDX-License-Identifier: MIT

#pragma once
#include <memory>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <unordered_map>
#include "mfem.hpp"
#include "physics_solver.hpp"
#include "axisymmetric_diffusion_integrator.hpp"
#include "input_parser.hpp"
#include "boundary_validation.hpp"
#include "gmsh_results_writer.hpp"

class ElectrostaticSolver : public PhysicsSolver {
	ModelType type = ModelType::Axisymmetric;

	// Primary Spaces
	std::unique_ptr<mfem::H1_FECollection> fec;
	std::unique_ptr<mfem::FiniteElementSpace> fespace;
	std::unique_ptr<mfem::GridFunction> x; // Electric Potential (V)

	// Physics
	std::unique_ptr<mfem::PWConstCoefficient> epsilon_coeff;
	std::unique_ptr<mfem::LinearForm> b;
	std::unique_ptr<mfem::BilinearForm> a;

	std::unique_ptr<mfem::SparseMatrix> K0;
	std::unique_ptr<mfem::DenseMatrix> C; // Coupling Matrix for terminals

	// Cached constrained system + factorization. The matrix is identical for every
	// solve (same bilinear form and essential DOFs), so we factor once and reuse it
	// for all scenarios / coupling columns. A_op must outlive `umf` because the
	// UMFPack solver keeps a pointer to the matrix it factored.
	mfem::OperatorPtr A_op;
#ifdef MFEM_USE_SUITESPARSE
	std::unique_ptr<mfem::UMFPackSolver> umf;
#endif

	mfem::Array<int> ess_bdr;
	mfem::Array<int> ess_tdof_list;
	std::unordered_map<std::string, mfem::Array<int>> terminal_markers; // Terminal name to boundary marker mapping

public:
	ElectrostaticSolver(mfem::Mesh& m, const json& c) : PhysicsSolver(m, c) {}

	void Setup() override {
		// Config & solver type
		InputParser parser(config_json);
		config = parser.GetProblemConfig();

		int order = config.Order;
		const int dim = mesh.Dimension();

		// Axisymmetric or Planar
		type = config.ModelType;

		// FE spaces
		fec = std::make_unique<mfem::H1_FECollection>(order, dim);
		fespace = std::make_unique<mfem::FiniteElementSpace>(&mesh, fec.get());
		x = std::make_unique<mfem::GridFunction>(fespace.get());

		// Material Properties (Permittivity)
		std::vector<Material> materials = config.Materials;
		mfem::Vector epsilon_values(mesh.attributes.Max());
		epsilon_values = 0.0;

		for (auto& region : config.Regions) {
			for (auto attribute_id : region.AttributeIds) {
				if (attribute_id > 0 && attribute_id <= mesh.attributes.Max()) {
					auto& material = materials[region.Material];
					epsilon_values[attribute_id - 1] = material.RelPermittivity * Constants::EPSILON_0;
				}
			}
		}

		epsilon_coeff = std::make_unique<mfem::PWConstCoefficient>(epsilon_values);

		// Boundary Attributes
		std::vector<std::pair<mfem::Array<int>, double>> bcs;
		for (const auto& bc : config.BoundaryConditions) {
			auto marker = MarkerFromAttrs(bc.AttributeIds);
			bcs.push_back({ marker, bc.Value });
		}

		terminal_markers = std::unordered_map<std::string, mfem::Array<int>>();
		for (const auto& term : config.Terminals) {
			terminal_markers[term.Name] = MarkerFromAttrs(term.AttributeIds);
		}

		BoundaryConditionValidator validator(mesh, *fespace);
		validator.ValidateBoundaryConditions(bcs, terminal_markers, /*allow_overlap=*/false);

		ess_bdr.SetSize(mesh.bdr_attributes.Max());
		ess_bdr = 0;
		for (const auto& [marker, val] : bcs) { // closures
			for (int i = 0; i < marker.Size(); ++i) {
				if (marker[i]) ess_bdr[i] = 1;
			}
		}
		for (const auto& [name, marker] : terminal_markers) { // voltage terminals
			for (int i = 0; i < marker.Size(); ++i) {
				if (marker[i]) ess_bdr[i] = 1;
			}
		}

		// Assemble Stiffness Matrix
		a = std::make_unique<mfem::BilinearForm>(fespace.get());

		if (type == ModelType::Axisymmetric) {
			// Solves: Div( r * eps * Grad(V) ) = 0
			// Integrator handles 'r' weight and 'epsilon'
			a->AddDomainIntegrator(new AxisymmetricDiffusionIntegrator(*epsilon_coeff));
		}
		else {
			// Solves: Div( eps * Grad(V) ) = 0
			// Standard Cartesian Laplacian
			a->AddDomainIntegrator(new mfem::DiffusionIntegrator(*epsilon_coeff));
		}

		a->Assemble();

		K0 = std::make_unique<mfem::SparseMatrix>(a->SpMat());

		// Linear Form (RHS)
		b = std::make_unique<mfem::LinearForm>(fespace.get());
		*b = 0.0;

		fespace->GetEssentialTrueDofs(ess_bdr, ess_tdof_list);
	}

	void ImprintScenario(const Scenario& sc) {
		*x = 0.0; // Reset solution for new scenario
		for (const auto& bc : config.BoundaryConditions) {
			if (bc.Value != 0.0) {
				auto marker = MarkerFromAttrs(bc.AttributeIds);
				mfem::ConstantCoefficient c(bc.Value);
				x->ProjectBdrCoefficient(c, marker);
			}
		}
		for (const auto& term : config.Terminals) {
			if (term.Excitation == Quantity::Voltage) {
				for (const auto& exc : sc.Excitations) {
					if (!exc.Floating) {
						if (term.Name == exc.TerminalName) {
							auto marker = MarkerFromAttrs(term.AttributeIds);
							mfem::ConstantCoefficient c(exc.Value);
							x->ProjectBdrCoefficient(c, marker);
						}
					}
				}
			}
			else
			{
				std::cerr << "Excitation type not supported in ElectrostaticSolver. Skipping terminal " << term.Name << ".\n";
				continue;
			}
		}
	}

	void Run() override {
		if (config.StudyType == StudyType::CouplingMatrix) {
			const int num_terminals = config.Terminals.size();
			C = std::make_unique<mfem::DenseMatrix>(num_terminals, num_terminals);
			*C = 0.0;
			mfem::ConstantCoefficient one(1.0); // Unit drive on terminal i
			mfem::Vector Q(fespace->GetVSize());
			// Axisymmetric integrators omit the global 2*pi (charge is integrated per
			// radian), so scale to physical Coulombs. Planar 2D is per-unit-depth.
			const double charge_scale = (type == ModelType::Axisymmetric) ? Constants::TWO_PI : 1.0;
			// One solve per conductor: a 1 V drive on terminal i yields the full
			// i-th column of C. The constrained matrix is identical for every drive,
			// so SolveSystem() factors it once and only back-substitutes here.
			for (int i = 0; i < num_terminals; ++i) {
				const auto& term_i = config.Terminals[i];
				*x = 0.0;
				x->ProjectBdrCoefficient(one, terminal_markers[term_i.Name]);
				*b = 0.0;
				SolveSystem(); // factor-once + back-substitute; recovers full *x

				K0->Mult(*x, Q); // nodal charges Q = K0 * V (K0 is the pre-BC matrix)

				for (int k = 0; k < num_terminals; ++k) {
					const auto& term_k = config.Terminals[k];
					mfem::Array<int> vdofs_k; // marker array of size ndofs
					fespace->GetEssentialVDofs(terminal_markers[term_k.Name], vdofs_k);
					double Qk = 0.0;
					// Loop through all dofs and sum the charges corresponding to the vdofs of terminal k.
					for (int n = 0; n < vdofs_k.Size(); ++n) {
						if (vdofs_k[n]) Qk += Q(n);
					}
					// Drive is exactly 1 V, so C(k,i) = Q_k / V_i = Q_k (scaled).
					(*C)(k, i) = charge_scale * Qk;
				}
			}
		}
		else {
			for (const auto& sc : config.Scenarios) {
				ImprintScenario(sc);
				SolveSystem();
				SaveScenario(sc.Name);
			}
		}
	}

	void SolveSystem() {
		mfem::Vector B, X;
		a->FormLinearSystem(ess_tdof_list, *x, *b, A_op, X, B);
#ifdef MFEM_USE_SUITESPARSE
		if (!umf) {
			umf = std::make_unique<mfem::UMFPackSolver>();
			umf->Control[UMFPACK_PRL] = 1;
			umf->SetOperator(*A_op); // factor once; reused for every subsequent RHS
		}
		umf->Mult(B, X);
#else
		auto* sp = dynamic_cast<mfem::SparseMatrix*>(A_op.Ptr());
		MFEM_ASSERT(sp, "Expected SparseMatrix operator from FormLinearSystem.");
		mfem::GSSmoother M(*sp);
		mfem::PCG(*A_op, M, B, X, config.SolverPrintLevel,
			config.SolverMaxIter, config.SolverTolerance, 0.0);
#endif
		a->RecoverFEMSolution(X, *b, *x);
	}

	void SaveScenario(const std::string& scenario_name) override {
		if (config.OutputParaview) {
			WriteParaviewResultsFile(scenario_name);
		}
		if (config.OutputGmsh) {
			WriteGmshResultsFile(scenario_name);
		}
	}

	void SaveStudy() override {
		if (config.StudyType == StudyType::CouplingMatrix) {
			WriteCouplingMatrix();
		}
	}

private:
	// Writes the computed Maxwell (short-circuit) capacitance matrix.
	// C(k,i) is the charge induced on conductor k when conductor i is held at
	// 1 V and all other conductors at 0 V. The matrix is symmetric; diagonals
	// are positive and off-diagonals negative. Each row sums to that
	// conductor's capacitance to the grounded/closure boundary. 
	void WriteCouplingMatrix() {
		if (!C) {
			std::cerr << "WriteCouplingMatrix: coupling matrix not computed.\n";
			return;
		}

		const int n = C->Height();
		const auto& terminals = config.Terminals;

		// Console summary.
		std::cout << "\n=== Capacitance Matrix [F] ===\n";
		std::cout << std::setw(18) << " ";
		for (int j = 0; j < n; ++j) {
			std::cout << std::setw(18) << terminals[j].Name;
		}
		std::cout << "\n";
		for (int k = 0; k < n; ++k) {
			std::cout << std::setw(18) << terminals[k].Name;
			for (int i = 0; i < n; ++i) {
				std::cout << std::setw(18) << std::scientific << std::setprecision(6) << (*C)(k, i);
			}
			std::cout << "\n";
		}

		// CSV next to the mesh (same path convention as the Gmsh/ParaView writers).
		namespace fs = std::filesystem;
		fs::path mesh_dir = fs::path(config.MeshPath).parent_path();
		fs::path out_path = mesh_dir / "capacitance_matrix.csv";

		std::ofstream ofs(out_path);
		if (!ofs) {
			std::cerr << "WriteCouplingMatrix: failed to open " << out_path.string() << "\n";
			return;
		}

		ofs << "Terminal";
		for (int j = 0; j < n; ++j) { ofs << "," << terminals[j].Name; }
		ofs << "\n";
		ofs << std::scientific << std::setprecision(12);
		for (int k = 0; k < n; ++k) {
			ofs << terminals[k].Name;
			for (int i = 0; i < n; ++i) { ofs << "," << (*C)(k, i); }
			ofs << "\n";
		}

		std::cout << "Wrote " << out_path.string() << std::endl;
	}

	void WriteParaviewResultsFile(const std::string& scenario_name) {
		// Paraview output (should we make an option?)
		mfem::ParaViewDataCollection paraview(("results_electrostatic_" + scenario_name).c_str(), &mesh);
		paraview.SetLevelsOfDetail(1);
		paraview.RegisterField("V", x.get());

		// Electric Field: E = -Grad(V)
		mfem::L2_FECollection fec_l2(fec->GetOrder() - 1, mesh.Dimension());

		// Electric Field Vector Space
		mfem::FiniteElementSpace fespace_l2_vec(&mesh, &fec_l2, mesh.Dimension());
		mfem::GridFunction E(&fespace_l2_vec);

		mfem::GradientGridFunctionCoefficient grad_x(x.get());
		E.ProjectCoefficient(grad_x);
		E *= -1.0;  // E = -Grad(V)

		paraview.RegisterField("E", &E);

		// Permittivity
		mfem::FiniteElementSpace fespace_l2_scalar(&mesh, &fec_l2);
		mfem::GridFunction eps_gf(&fespace_l2_scalar);
		eps_gf.ProjectCoefficient(*epsilon_coeff);
		paraview.RegisterField("Permittivity", &eps_gf);

		paraview.SetCycle(0);
		paraview.SetTime(0.0);
		paraview.Save();
	}

	void WriteGmshResultsFile(const std::string& scenario_name) {
		namespace fs = std::filesystem;

		// Resolve output path: explicit override > derive from mesh path.
		// A relative override is resolved against the mesh's directory rather
		// than the process CWD, so the results land next to the input mesh
		// regardless of where the executable was launched from (VS launches
		// from out\build\<cfg>\, CTest from the build dir, etc.).
		fs::path mesh_dir = fs::path(config.MeshPath).parent_path();
		fs::path out_path = mesh_dir / (scenario_name + ".results.msh");

		// Refinement factor: explicit override > Order (>=1).
		int ref_factor = (config.ExportRefine > 0) ? config.ExportRefine
			: std::max(1, config.Order);

		// Build refined export mesh. The ClosedUniform basis matches the
		// ParaView export path and gives uniformly distributed sub-vertices.
		mfem::Mesh export_mesh(&mesh, ref_factor, mfem::BasisType::ClosedUniform);

		const int dim = export_mesh.Dimension();
		const int sdim = export_mesh.SpaceDimension();

		// Linear H1 for V on the refined mesh (one DOF per vertex).
		mfem::H1_FECollection fec_h1_lin(1, dim);
		mfem::FiniteElementSpace fes_V(&export_mesh, &fec_h1_lin);
		mfem::GridFunction V_h(&fes_V);
		{
			mfem::GridFunctionCoefficient phi_coeff(x.get());
			V_h.ProjectCoefficient(phi_coeff);
		}

		// Linear L2 vector space for E.
		mfem::L2_FECollection fec_l2_lin(1, dim);
		mfem::FiniteElementSpace fes_E(&export_mesh, &fec_l2_lin, sdim);
		mfem::GridFunction E_h(&fes_E);
		{
			mfem::GradientGridFunctionCoefficient grad_phi(x.get());
			mfem::ScalarVectorProductCoefficient minus_grad(-1.0, grad_phi);
			E_h.ProjectCoefficient(minus_grad);
		}

		std::vector<gmsh_results::View> views;
		views.push_back(gmsh_results::MakeScalarNodeView("V", V_h));
		views.push_back(gmsh_results::MakeVectorElementNodeView("E", E_h));
		views.push_back(gmsh_results::MakeMagnitudeElementNodeView("|E|", E_h));

		// TODO: per-surface tangential field views "Et_<surface>" from
		// creep-surface boundary attributes (FaceElementTransformations +
		// E - (E.n) n). The C# loader simply ignores missing keys.

		gmsh_results::WriteGmshResults(out_path.string(), export_mesh, fes_E, views);

		std::cout << "Wrote " << out_path.string() << std::endl;
	}
};