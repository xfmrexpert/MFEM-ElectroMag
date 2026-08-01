// Copyright (c) 2026 T. C. Raymond
// SPDX-License-Identifier: MIT

#pragma once
#include <memory>
#include <iostream>
#include <limits>
#include <unordered_map>
#include <vector>
#include "mfem.hpp"
#include "physics_solver.hpp"
#include "axisymmetric_diffusion_integrator.hpp"
#include "boundary_validation.hpp"
#include "gmsh_results_writer.hpp"
#include "amr_support.hpp"

class ElectrostaticSolver : public PhysicsSolver {
	// Primary Spaces
	std::unique_ptr<mfem::GridFunction> x; // Electric Potential (V)

	// Physics
	std::unique_ptr<mfem::PWConstCoefficient> epsilon_coeff;

	std::unique_ptr<mfem::LinearForm> b;
	mfem::Vector neumann_rhs;
	std::unique_ptr<mfem::BilinearForm> a;

	std::unique_ptr<mfem::SparseMatrix> K0;
	std::unique_ptr<mfem::DenseMatrix> C; // Coupling Matrix for terminals

	// Cached constrained system for the CURRENT mesh. The matrix is identical
	// for every solve on a given mesh (same bilinear form and essential DOFs),
	// so it is assembled once per mesh in BuildOperators() and reused for all
	// scenarios / coupling columns. AMR refinement rebuilds it via BuildOperators().
	mfem::OperatorPtr A_op;

	std::unordered_map<std::string, mfem::Array<int>> terminal_markers; // Terminal name to boundary marker mapping

public:
	ElectrostaticSolver(mfem::Mesh& m, const ProblemConfig& c) : PhysicsSolver(m, c) {}

	void Setup() override {
		int order = config.Order;
		const int dim = mesh.Dimension();

		// Axisymmetric or Planar
		geometry = config.GeometryType;

		// Reject negative radii. Electrostatics needs no axis condition: the
		// natural condition dV/dr = 0 on r = 0 is already the correct one.
		ValidateAxisymmetricGeometry();

		// FE collection
		fec = std::make_unique<mfem::H1_FECollection>(order, dim);

		// Material Properties (Permittivity). epsilon_coeff is a PWConstCoefficient
		// keyed by mesh DOMAIN attribute. AMR refinement subdivides elements but
		// preserves their attributes, so this mapping is refinement-invariant.
	
		epsilon_coeff = MaterialCoefficient(0.0, [](const Material& m) {
			return m.RelPermittivity * Constants::EPSILON_0; });

		// Closures + voltage terminals are all essential (Dirichlet).
		closure_bcs = BuildClosureBcs();

		terminal_markers.clear();
		for (const auto& [term_name, term] : config.Terminals)
			terminal_markers[term_name] = MarkerFromGroup(term.EntityGroupName);

		std::vector<mfem::Array<int>> ess_markers;
		for (const auto& marker : DirichletClosureMarkers(closure_bcs)) {
			ess_markers.push_back(marker);
		}
		for (const auto& [name, marker] : terminal_markers) ess_markers.push_back(marker);
		ess_bdr = EssentialBdrFrom(ess_markers);

		// Build the FE space and everything bound to it for the starting mesh.
		BuildOperators();

		// Validate boundary conditions once. The check is over the (refinement-
		// invariant) mesh topology / attributes; it merely needs an FE space for
		// DOF queries, so running it after the first BuildOperators() is correct.
		BoundaryConditionValidator validator(mesh, *fespace);
		validator.ValidateBoundaryConditions(closure_bcs, terminal_markers, /*allow_overlap=*/false);
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
	void BuildOperators() override {
		fespace = std::make_unique<mfem::FiniteElementSpace>(&mesh, fec.get());
		
		x = std::make_unique<mfem::GridFunction>(fespace.get());
		*x = 0.0;

		// Assemble Stiffness Matrix
		a = std::make_unique<mfem::BilinearForm>(fespace.get());
		a->AddDomainIntegrator(MakeStiffnessIntegrator()); // a takes ownership
		a->Assemble();

		// Snapshot the UNCONSTRAINED stiffness matrix (used for charge Q = K0*x)
		// before FormSystemMatrix eliminates the essential DOFs from a's SpMat.
		K0 = std::make_unique<mfem::SparseMatrix>(a->SpMat());

		// Linear Form (RHS)
		b = std::make_unique<mfem::LinearForm>(fespace.get());
		neumann_rhs = AssembleNeumannBoundaryLoad();

		fespace->GetEssentialTrueDofs(ess_bdr, ess_tdof_list);

		// Form the constrained system operator. The eliminated-column part
		// (mat_e, used to build each scenario's RHS) is bound to A_op, which is
		// reused for every scenario's solve on this mesh.
		a->FormSystemMatrix(ess_tdof_list, A_op);
	}

	// Create the domain diffusion integrator matching the active geometry. Used
	// both by the solve (owned by 'a') and the AMR error estimator (a separate,
	// independently-owned instance), so the estimated error is consistent with
	// the assembled operator - axisymmetric (2*pi*r, eps) or planar (eps).
	mfem::BilinearFormIntegrator* MakeStiffnessIntegrator() const {
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
	// MPI-only). A dedicated integrator instance (separate from a's) and an H1
	// vector flux space are constructed here per call; SetFluxAveraging(1) keeps
	// the recovered flux from smoothing across material-attribute interfaces so
	// per-region permittivity discontinuities are respected. SetWithCoeff(true)
	// makes the flux eps*grad(V), consistent with the integrator's energy norm.
	//
	// @param combined  Output: per-element max error indicator (sized to NE).
	// @return Global error sqrt(sum_k combined_k^2).
	double EstimateCombinedError(mfem::Vector& combined) override {
		const int sdim = mesh.SpaceDimension();
		std::unique_ptr<mfem::BilinearFormIntegrator> flux_integ(MakeStiffnessIntegrator());
		mfem::FiniteElementSpace flux_fes(&mesh, fec.get(), sdim);
		mfem::ZienkiewiczZhuEstimator estimator(*flux_integ, *x, flux_fes);
		estimator.SetWithCoeff(true);     // flux = eps * grad(V)
		estimator.SetFluxAveraging(1);    // do not average across attribute interfaces

		return EstimateScenarioMaximumError(combined,
			[this](const Scenario& scenario) {
				ImprintScenario(scenario);
				SolveSystem();
			},
			[&estimator](mfem::Vector& current) {
				estimator.Reset();
				current = estimator.GetLocalErrors();
			});
	}

	// Peak field magnitude |E| = |grad(V)| over the current solution *x, sampled
	// at element nodes. Used as an AMR convergence diagnostic (peak |E| at a
	// conductor corner should settle as refinement resolves the singularity).
	// Reflects whichever solution currently lives in *x (the last one solved).
	double ComputePeakFieldMagnitude() const override {
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

	void ImprintScenario(const Scenario& sc) {
		*x = 0.0; // Reset solution for new scenario
		*b = neumann_rhs;
		
		for (const auto& bc : config.BoundaryConditions) {
			if (bc.Type == BoundaryConditionType::Dirichlet && bc.Value != 0.0) {
				auto marker = MarkerFromGroup(bc.EntityGroupName);
				mfem::ConstantCoefficient c(bc.Value);
				x->ProjectBdrCoefficient(c, marker);
			}
		}
		for (const auto& [term_name, term] : config.Terminals) {
			if (term.Excitation == Quantity::Voltage) {
				for (const auto& exc : sc.Excitations) {
					if (!exc.Floating) {
						if (term_name == exc.TerminalName) {
							auto marker = MarkerFromGroup(term.EntityGroupName);
							mfem::ConstantCoefficient c(exc.Value);
							x->ProjectBdrCoefficient(c, marker);
						}
					}
				}
			}
			else
			{
				Reporter().Warning("Excitation type not supported in ElectrostaticSolver. "
					"Skipping terminal " + term_name + ".");
				continue;
			}
		}
	}

	// Solve + save on the CURRENT mesh/operators. Both analysis types flow through
	// ONE imprint -> solve loop over BuildSolveScenarios(); only the intrinsic
	// post-solve action differs: CouplingMatrix gathers an induced-charge column
	// into C, Field saves the scenario's fields. AMR calls this once more on the
	// converged mesh so the exported C / fields match the exported mesh.
	void RunOnCurrentMesh() override {
		const bool coupling = (config.AnalysisType == AnalysisType::CouplingMatrix);
		if (coupling) {
			const int n = static_cast<int>(config.Terminals.size());
			C = std::make_unique<mfem::DenseMatrix>(n, n);
			*C = 0.0;
		}

		int col = 0;
		for (const auto& [name, sc] : BuildSolveScenarios()) {
			auto operation = Reporter().Start("scenario '" + name + "'");
			ImprintScenario(sc);
			SolveSystem();
			if (coupling) { GatherChargeColumn(col++); }
			else          { SaveScenario(name); }
		}
	}

	void SolveSystem() {
		auto operation = Reporter().Start("linear system solve");
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
		mfem::PCG(*A_op, M, B, X, Reporter().SolverPrintLevel(config.SolverPrintLevel),
			config.SolverMaxIter, config.SolverTolerance, 0.0);
		a->RecoverFEMSolution(X, *b, *x);
	}

	// Post-solve field recovery: the potential V, the field E = -grad(V), and the
	// per-region permittivity. Serialization is handled by the base class.
	FieldExportSet CollectExportFields() const override {
		FieldExportSet fields;
		fields.AddPrimaryScalar("V", *x);

		auto grad = std::make_unique<mfem::GradientGridFunctionCoefficient>(x.get());
		auto& grad_ref = fields.Own(std::move(grad));
		fields.AddVector("E",
			std::make_unique<mfem::ScalarVectorProductCoefficient>(-1.0, grad_ref));

		fields.AddScalar("Permittivity", *epsilon_coeff);
		return fields;
	}

	void SaveAnalysis() override {
		if (config.AnalysisType == AnalysisType::CouplingMatrix) {
			WriteCouplingMatrix();
		}
	}

private:
	// CouplingMatrix post-solve action for one column: with the just-solved
	// potential in *x (terminal `col` driven at 1 V, the rest grounded by the
	// synthesized scenario), gather the induced charge Q = K0 * x onto every
	// conductor's boundary DOFs and write column `col` of C. Off-diagonals are
	// negative, diagonals positive; 2*pi scaling accounts for the axisymmetric
	// measure omitted in the integrators. Terminal order matches
	// BuildSolveScenarios() / WriteCouplingMatrix() (config.Terminals order).
	void GatherChargeColumn(int col) {
		std::vector<const std::pair<const std::string, Terminal>*> terms;
		terms.reserve(config.Terminals.size());
		for (const auto& kv : config.Terminals) terms.push_back(&kv);

		mfem::Vector Q(fespace->GetVSize());
		K0->Mult(*x, Q);

		const double charge_scale = (geometry == GeometryType::Axisymmetric) ? Constants::TWO_PI : 1.0;
		for (int k = 0; k < static_cast<int>(terms.size()); ++k) {
			mfem::Array<int> vdofs_k;
			fespace->GetEssentialVDofs(terminal_markers[terms[k]->first], vdofs_k);
			double Qk = 0.0;
			for (int n = 0; n < vdofs_k.Size(); ++n) {
				if (vdofs_k[n]) Qk += Q(n);
			}
			(*C)(k, col) = charge_scale * Qk;
		}
	}

	// Writes the computed Maxwell (short-circuit) capacitance matrix.
	// C(k,i) is the charge induced on conductor k when conductor i is held at
	// 1 V and all other conductors at 0 V. The matrix is symmetric; diagonals
	// are positive and off-diagonals negative. Each row sums to that
	// conductor's capacitance to the grounded/closure boundary. 
	void WriteCouplingMatrix() {
		if (!C) {
			Reporter().Warning("WriteCouplingMatrix: coupling matrix not computed.");
			return;
		}

		SaveCouplingMatrix(*C, "Capacitance Matrix [F]", "capacitance_matrix.csv");
	}
};