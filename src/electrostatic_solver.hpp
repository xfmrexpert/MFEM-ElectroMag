// Copyright (c) 2026 T. C. Raymond
// SPDX-License-Identifier: MIT

#pragma once
#include <memory>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <unordered_map>
#include <vector>
#include "mfem.hpp"
#include "physics_solver.hpp"
#include "axisymmetric_diffusion_integrator.hpp"
#include "input_parser.hpp"
#include "boundary_validation.hpp"
#include "gmsh_results_writer.hpp"
#include "amr_support.hpp"

class ElectrostaticSolver : public PhysicsSolver {
	GeometryType geometry = GeometryType::Axisymmetric;

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

	// Cached constrained system for the CURRENT mesh. The matrix is identical
	// for every solve on a given mesh (same bilinear form and essential DOFs),
	// so it is assembled once per mesh in BuildOperators() and reused for all
	// scenarios / coupling columns. AMR refinement rebuilds it via
	// BuildOperators().
	mfem::OperatorPtr A_op;

	mfem::Array<int> ess_bdr;
	mfem::Array<int> ess_tdof_list;
	std::unordered_map<std::string, mfem::Array<int>> terminal_markers; // Terminal name to boundary marker mapping

public:
	// One AMR iteration's diagnostics, recorded during RunAdaptive(). Exposed for
	// tests (convergence / conformity assertions) and console logging.
	struct AmrIterationInfo {
		long   true_dofs;     // global true DOFs on that iteration's mesh
		double global_error;  // sqrt(sum_k eta_k^2), combined over scenarios
		double peak_absE;     // max |E| sampled over the mesh
	};

private:
	std::vector<AmrIterationInfo> amr_history;

public:
	ElectrostaticSolver(mfem::Mesh& m, const json& c) : PhysicsSolver(m, c) {}

	void Setup() override {
		// Config & solver type
		InputParser parser(config_json);
		config = parser.GetProblemConfig();

		int order = config.Order;
		const int dim = mesh.Dimension();

		// Axisymmetric or Planar
		geometry = config.GeometryType;

		// FE collection. Depends only on (order, dim), NOT on the mesh, so it is
		// refinement-invariant and built once here; BuildOperators() reuses it.
		fec = std::make_unique<mfem::H1_FECollection>(order, dim);

		// Material Properties (Permittivity). epsilon_coeff is a PWConstCoefficient
		// keyed by mesh DOMAIN attribute. AMR refinement subdivides elements but
		// preserves their attributes, so this mapping is refinement-invariant.
		std::vector<Material> materials = config.Materials;
		mfem::Vector epsilon_values(mesh.attributes.Max());
		epsilon_values = 0.0;

		for (auto& region : config.Regions) {
			const std::string& group_name = region.EntityGroupName;
			const EntityGroup& group = config.EntityGroups.at(group_name);
			for (auto attribute_id : group.AttributeIds) {
				if (attribute_id > 0 && attribute_id <= mesh.attributes.Max()) {
					auto& material = materials[region.Material];
					epsilon_values[attribute_id - 1] = material.RelPermittivity * Constants::EPSILON_0;
				}
			}
		}

		epsilon_coeff = std::make_unique<mfem::PWConstCoefficient>(epsilon_values);

		// Boundary Attributes. Markers are keyed by mesh BOUNDARY attribute, which
		// refinement preserves, so bcs / terminal_markers / ess_bdr are all
		// refinement-invariant and built once here.
		std::vector<std::pair<mfem::Array<int>, double>> bcs;
		for (const auto& bc : config.BoundaryConditions) {
			const std::string& group_name = bc.EntityGroupName;
			const EntityGroup& group = config.EntityGroups.at(group_name);
			auto marker = MarkerFromAttrs(group.AttributeIds);
			bcs.push_back({ marker, bc.Value });
		}

		terminal_markers = std::unordered_map<std::string, mfem::Array<int>>();
		for (const auto& [term_name, term] : config.Terminals) {
			const std::string& group_name = term.EntityGroupName;
			const EntityGroup& group = config.EntityGroups.at(group_name);
			terminal_markers[term_name] = MarkerFromAttrs(group.AttributeIds);
		}

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

		// Build the FE space and everything bound to it for the starting mesh.
		BuildOperators();

		// Validate boundary conditions once. The check is over the (refinement-
		// invariant) mesh topology / attributes; it merely needs an FE space for
		// DOF queries, so running it after the first BuildOperators() is correct.
		BoundaryConditionValidator validator(mesh, *fespace);
		validator.ValidateBoundaryConditions(bcs, terminal_markers, /*allow_overlap=*/false);
	}

	// (Re)build the FE space and every object bound to it for the CURRENT mesh.
	// Called once from Setup() and again after each AMR refinement. The
	// refinement-invariant data (fec, epsilon_coeff, ess_bdr, terminal_markers)
	// persists across calls and is reused.
	//
	// The constrained system matrix is assembled here; within a single mesh that
	// matrix is reused for every scenario / coupling column (the bilinear form
	// and essential-DOF set do not change between solves). AMR refinement
	// invalidates the mesh, so this is re-run to rebuild on the new mesh.
	void BuildOperators() {
		fespace = std::make_unique<mfem::FiniteElementSpace>(&mesh, fec.get());
		x = std::make_unique<mfem::GridFunction>(fespace.get());
		*x = 0.0;

		// Assemble Stiffness Matrix
		a = std::make_unique<mfem::BilinearForm>(fespace.get());
		a->AddDomainIntegrator(MakeDiffusionIntegrator()); // a takes ownership
		a->Assemble();

		// Snapshot the UNCONSTRAINED stiffness matrix (used for charge Q = K0*x)
		// before FormSystemMatrix eliminates the essential DOFs from a's SpMat.
		K0 = std::make_unique<mfem::SparseMatrix>(a->SpMat());

		// Linear Form (RHS)
		b = std::make_unique<mfem::LinearForm>(fespace.get());

		fespace->GetEssentialTrueDofs(ess_bdr, ess_tdof_list);

		// Form the constrained system operator. The eliminated-column part
		// (mat_e, used to build each scenario's RHS) is bound to A_op, which is
		// reused for every scenario's solve on this mesh.
		a->FormSystemMatrix(ess_tdof_list, A_op);
	}

	// Create the domain diffusion integrator matching the active geometry. Used
	// both by the solve (owned by `a`) and the AMR error estimator (a separate,
	// independently-owned instance), so the estimated error is consistent with
	// the assembled operator - axisymmetric (2*pi*r, eps) or planar (eps).
	mfem::BilinearFormIntegrator* MakeDiffusionIntegrator() const {
		if (geometry == GeometryType::Axisymmetric) {
			// Solves: Div( r * eps * Grad(V) ) = 0; integrator handles 'r' and 'eps'.
			return new AxisymmetricDiffusionIntegrator(*epsilon_coeff);
		}
		// Solves: Div( eps * Grad(V) ) = 0; standard Cartesian Laplacian.
		return new mfem::DiffusionIntegrator(*epsilon_coeff);
	}

	// Estimate per-element error on the CURRENT mesh, folding every scenario /
	// coupling column into a single indicator via the element-wise maximum
	//     eta_k = max over solves s of eta_k(s),
	// so one shared mesh is refined for all scenarios (spec: identical
	// $Nodes/$Elements across every <scenario>.results.msh).
	//
	// Uses the serial recovery-based ZienkiewiczZhuEstimator (the L2 variant is
	// MPI-only). A dedicated integrator instance (separate from `a`'s) and an H1
	// vector flux space are constructed here per call; SetFluxAveraging(1) keeps
	// the recovered flux from smoothing across material-attribute interfaces so
	// per-region permittivity discontinuities are respected. SetWithCoeff(true)
	// makes the flux eps*grad(V), consistent with the integrator's energy norm.
	//
	// @param combined  Output: per-element max error indicator (sized to NE).
	// @return Global error sqrt(sum_k combined_k^2).
	double EstimateCombinedError(mfem::Vector& combined) {
		const int ne = mesh.GetNE();
		combined.SetSize(ne);
		combined = 0.0;

		const int sdim = mesh.SpaceDimension();
		std::unique_ptr<mfem::BilinearFormIntegrator> flux_integ(MakeDiffusionIntegrator());
		mfem::FiniteElementSpace flux_fes(&mesh, fec.get(), sdim);
		mfem::ZienkiewiczZhuEstimator estimator(*flux_integ, *x, flux_fes);
		estimator.SetWithCoeff(true);     // flux = eps * grad(V)
		estimator.SetFluxAveraging(1);    // do not average across attribute interfaces

		auto fold_current_solution = [&]() {
			estimator.Reset(); // force recompute: same mesh, new solution in *x
			const mfem::Vector& errs = estimator.GetLocalErrors();
			for (int k = 0; k < ne; ++k) {
				if (errs(k) > combined(k)) { combined(k) = errs(k); }
			}
		};

		if (config.AnalysisType == AnalysisType::CouplingMatrix) {
			// One unit excitation per terminal (drive i = 1 V, rest = 0).
			mfem::ConstantCoefficient one(1.0);
			for (const auto& [term_name, marker] : terminal_markers) {
				*x = 0.0;
				x->ProjectBdrCoefficient(one, terminal_markers.at(term_name));
				*b = 0.0;
				SolveSystem();
				fold_current_solution();
			}
		}
		else {
			for (const auto& [sc_name, sc] : config.Scenarios) {
				ImprintScenario(sc);
				SolveSystem();
				fold_current_solution();
			}
		}

		double sum_sq = 0.0;
		for (int k = 0; k < ne; ++k) { sum_sq += combined(k) * combined(k); }
		return std::sqrt(sum_sq);
	}

	// Peak field magnitude |E| = |grad(V)| over the current solution *x, sampled
	// at element nodes. Used as an AMR convergence diagnostic (peak |E| at a
	// conductor corner should settle as refinement resolves the singularity).
	// Reflects whichever solution currently lives in *x (the last one solved).
	double ComputePeakFieldMagnitude() const {
		if (!x) { return 0.0; }
		double peak = 0.0;
		mfem::Vector grad;
		for (int e = 0; e < fespace->GetNE(); ++e) {
			const mfem::FiniteElement* fe = fespace->GetFE(e);
			mfem::ElementTransformation* T = fespace->GetElementTransformation(e);
			const mfem::IntegrationRule& nodes = fe->GetNodes();
			for (int i = 0; i < nodes.GetNPoints(); ++i) {
				const mfem::IntegrationPoint& ip = nodes.IntPoint(i);
				T->SetIntPoint(&ip);
				x->GetGradient(*T, grad); // grad(V); |E| = |grad(V)|
				const double mag = grad.Norml2();
				if (mag > peak) { peak = mag; }
			}
		}
		return peak;
	}

	// AMR per-iteration diagnostics from the most recent RunAdaptive(). Empty when
	// AMR is disabled. Consumed by the regression tests and useful for logging.
	const std::vector<AmrIterationInfo>& GetAmrHistory() const { return amr_history; }

	void ImprintScenario(const Scenario& sc) {
		*x = 0.0; // Reset solution for new scenario
		*b = 0.0; // Reset RHS for new scenario
		for (const auto& bc : config.BoundaryConditions) {
			if (bc.Value != 0.0) {
				const std::string& group_name = bc.EntityGroupName;
				const EntityGroup& group = config.EntityGroups.at(group_name);
				auto marker = MarkerFromAttrs(group.AttributeIds);
				mfem::ConstantCoefficient c(bc.Value);
				x->ProjectBdrCoefficient(c, marker);
			}
		}
		for (const auto& [term_name, term] : config.Terminals) {
			if (term.Excitation == Quantity::Voltage) {
				for (const auto& exc : sc.Excitations) {
					if (!exc.Floating) {
						if (term_name == exc.TerminalName) {
							const std::string& group_name = term.EntityGroupName;
							const EntityGroup& group = config.EntityGroups.at(group_name);
							auto marker = MarkerFromAttrs(group.AttributeIds);
							mfem::ConstantCoefficient c(exc.Value);
							x->ProjectBdrCoefficient(c, marker);
						}
					}
				}
			}
			else
			{
				std::cerr << "Excitation type not supported in ElectrostaticSolver. Skipping terminal " << term_name << ".\n";
				continue;
			}
		}
	}

	void Run() override {
		if (config.Amr.Enabled) {
			RunAdaptive();
		}
		else {
			RunFixed();
		}
	}

	// Solve + save on the CURRENT mesh/operators. This is the legacy Run() body:
	// CouplingMatrix assembles C from unit terminal excitations; Field solves and
	// saves every authored scenario. AMR calls this once more on the final mesh.
	void RunFixed() {
		if (config.AnalysisType == AnalysisType::CouplingMatrix) {
			const int num_terminals = static_cast<int>(config.Terminals.size());

			// Positional view over the (name-sorted) terminal map. C is inherently
			// positional, so snapshot the order once and index that.
			std::vector<const std::pair<const std::string, Terminal>*> terms;
			terms.reserve(config.Terminals.size());
			for (const auto& kv : config.Terminals) terms.push_back(&kv);

			C = std::make_unique<mfem::DenseMatrix>(num_terminals, num_terminals);
			*C = 0.0;
			mfem::ConstantCoefficient one(1.0);
			mfem::Vector Q(fespace->GetVSize());
			const double charge_scale = (geometry == GeometryType::Axisymmetric) ? Constants::TWO_PI : 1.0;

			for (int i = 0; i < num_terminals; ++i) {
				const std::string& term_name = terms[i]->first;
				*x = 0.0;
				x->ProjectBdrCoefficient(one, terminal_markers[term_name]);
				*b = 0.0;
				SolveSystem();

				K0->Mult(*x, Q);

				for (int k = 0; k < num_terminals; ++k) {
					const std::string& term_name_k = terms[k]->first;
					mfem::Array<int> vdofs_k;
					fespace->GetEssentialVDofs(terminal_markers[term_name_k], vdofs_k);
					double Qk = 0.0;
					for (int n = 0; n < vdofs_k.Size(); ++n) {
						if (vdofs_k[n]) Qk += Q(n);
					}
					(*C)(k, i) = charge_scale * Qk;
				}
			}
		}
		else {
			for (const auto& [sc_name, sc] : config.Scenarios) {
				ImprintScenario(sc);
				SolveSystem();
				SaveScenario(sc_name);
			}
		}
	}

	// h-adaptive loop: estimate combined (max-over-scenarios) error on the current
	// mesh, stop on iteration / DOF / error-tolerance caps, otherwise bulk-mark
	// (Dorfler) and refine conformingly, rebuild operators, and repeat. A final
	// RunFixed() on the converged mesh produces the saved fields / C matrix so the
	// exported results correspond exactly to the exported (refined) mesh.
	void RunAdaptive() {
		const AmrSettings& amr = config.Amr;
		amr_history.clear();

		const int max_it = std::max(1, amr.MaxIterations);
		for (int it = 0; it < max_it; ++it) {
			const long cdofs = fespace->GetTrueVSize();

			mfem::Vector errors;
			const double global_err = EstimateCombinedError(errors);
			const double peak_absE = ComputePeakFieldMagnitude();

			amr_history.push_back({ cdofs, global_err, peak_absE });

			std::cout << "AMR iteration " << it
				<< ": elements=" << mesh.GetNE()
				<< ", true_dofs=" << cdofs
				<< ", global_error=" << std::scientific << std::setprecision(6) << global_err
				<< ", peak|E|=" << peak_absE << std::endl;

			// Stopping criteria (any one stops): error tolerance, DOF budget, or
			// this being the last permitted iteration.
			if (amr.ErrorTolerance > 0.0 && global_err < amr.ErrorTolerance) {
				std::cout << "AMR: global error below tolerance. Stop." << std::endl;
				break;
			}
			if (amr.MaxDofs > 0 && cdofs > amr.MaxDofs) {
				std::cout << "AMR: reached the maximum number of DOFs. Stop." << std::endl;
				break;
			}
			if (it + 1 >= max_it) {
				std::cout << "AMR: reached the maximum number of iterations. Stop." << std::endl;
				break;
			}

			// Mark and refine conformingly (throws if the mesh cannot refine
			// without hanging nodes), then rebuild the FE space / operators.
			mfem::Array<int> marked;
			amr::MarkElementsDorfler(errors, amr.ErrorFraction, marked);
			if (marked.Size() == 0) {
				std::cout << "AMR: no elements marked for refinement. Stop." << std::endl;
				break;
			}
			amr::RefineConforming(mesh, marked);
			BuildOperators();
		}

		// Authoritative final solve on the final mesh: produces the saved fields
		// (Field) or the C matrix (CouplingMatrix) on the exported mesh.
		RunFixed();
	}

	void SolveSystem() {
		// The system matrix and its essential-DOF elimination (mat_e) were built for
		// the current mesh in BuildOperators(). FormLinearSystem here re-derives ONLY
		// this scenario's eliminated RHS from the freshly imprinted x/b (b was zeroed
		// in ImprintScenario, so no previous scenario's load survives in it); it
		// reuses the same already-eliminated operator. Each scenario is therefore
		// solved independently, and AMR re-runs BuildOperators() after refinement to
		// rebuild on the new mesh.
		mfem::Vector B, X;
		a->FormLinearSystem(ess_tdof_list, *x, *b, A_op, X, B);
		auto* sp = dynamic_cast<mfem::SparseMatrix*>(A_op.Ptr());
		MFEM_ASSERT(sp, "Expected SparseMatrix operator from FormLinearSystem.");
		mfem::GSSmoother M(*sp);
		mfem::PCG(*A_op, M, B, X, config.SolverPrintLevel,
			config.SolverMaxIter, config.SolverTolerance, 0.0);
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

	void SaveAnalysis() override {
		if (config.AnalysisType == AnalysisType::CouplingMatrix) {
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

		std::vector<const std::pair<const std::string, Terminal>*> terms;
		terms.reserve(config.Terminals.size());
		for (const auto& kv : config.Terminals) terms.push_back(&kv);

		// Console summary.
		std::cout << "\n=== Capacitance Matrix [F] ===\n";
		std::cout << std::setw(18) << " ";
		for (int j = 0; j < n; ++j) {
			std::cout << std::setw(18) << terms[j]->first;
		}
		std::cout << "\n";
		for (int k = 0; k < n; ++k) {
			std::cout << std::setw(18) << terms[k]->first;
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
		for (int j = 0; j < n; ++j) { ofs << "," << terms[j]->first; }
		ofs << "\n";
		ofs << std::scientific << std::setprecision(12);
		for (int k = 0; k < n; ++k) {
			ofs << terms[k]->first;
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
		//
		// This is a VISUALIZATION sub-sampling factor (p-detail): each base element
		// is split ref_factor times per dimension so a high-Order (curved) solution
		// renders smoothly on the linear export mesh. It depends only on the
		// (refinement-invariant) polynomial Order, NOT on the mesh, so it composes
		// orthogonally with AMR: AMR performs h-refinement (more base elements) while
		// this performs p-sub-sampling (sub-vertices per element). AMR does not
		// compound ref_factor, so the export element count is the AMR element count
		// times ref_factor^dim - i.e. for the common Order==1 case ref_factor==1 and
		// the export mesh equals the AMR mesh exactly.
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