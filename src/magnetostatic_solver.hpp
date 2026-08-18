// Copyright (c) 2026 T. C. Raymond
// SPDX-License-Identifier: MIT

#pragma once

#include <memory>
#include <cmath>
#include <algorithm>
#include <sstream>

#include "mfem.hpp"
#include "physics_solver.hpp"
#include "axisymmetric_curl_curl_integrator.hpp"
#include "axisymmetric_lf_integrator.hpp"
#include "magnetic_field_coefficient.hpp"
#include "boundary_validation.hpp"
#include "constants.hpp"
#include "gmsh_results_writer.hpp"

class MagnetostaticSolver : public PhysicsSolver
{
private:
	// Resources (order of declaration = order of destruction)
	std::unique_ptr<mfem::GridFunction> A; // A_phi (axisym) or A_z (planar scalar)

	std::unique_ptr<mfem::PWConstCoefficient> nu_coeff; // nu=1/mu (reluctivity)
	std::unique_ptr<mfem::PWConstCoefficient> j_coeff; // J_phi (axisym) or J (planar scalar src)

	std::unique_ptr<mfem::LinearForm> b;
	std::unique_ptr<mfem::BilinearForm> a;
	mfem::Vector neumann_rhs;

	// Cached constrained system for the CURRENT mesh. The matrix is identical
	// for every solve on a given mesh (same bilinear form and essential DOFs),
	// so it is assembled once per mesh in BuildOperators() and reused for all
	// scenarios / coupling columns. AMR refinement rebuilds it via BuildOperators().
	mfem::OperatorHandle A_op;

	std::unique_ptr<mfem::DenseMatrix> L; // Inductance matrix (coupling matrix) for the current mesh

public:
	MagnetostaticSolver(mfem::Mesh& m, const ProblemConfig& c) : PhysicsSolver(m, c) {}

	void Setup() override
	{
		int order = config.Order;
		const int dim = mesh.Dimension();

		// Axisymmetric or Planar
		geometry = config.GeometryType;

		// Reject negative radii and record whether the domain reaches r = 0.
		ValidateAxisymmetricGeometry();

		// FE collection
		fec = std::make_unique<mfem::H1_FECollection>(order, dim);

		// Material Properties (Reluctivity nu = 1/mu), keyed by mesh DOMAIN attribute.
		nu_coeff = MaterialCoefficient(1.0 / Constants::MU_0, [](const Material& m) {
			return 1.0 / (Constants::MU_0 * m.RelPermeability); });

		closure_bcs = BuildClosureBcs();
		ValidateMagneticAxisBoundaryValues();

		std::vector<mfem::Array<int>> ess_markers =
			DirichletClosureMarkers(closure_bcs);
		ess_bdr = EssentialBdrFrom(ess_markers);

		// Axis regularity: the setup-time mesh inspection identifies a dedicated
		// r=0 boundary attribute. Merge it into the essential marker so A_phi = 0.
		if (geometry == GeometryType::Axisymmetric)
		{
			MFEM_VERIFY(ess_bdr.Size() == axisymmetric_mesh.axis_boundary.Size(),
				"Axis boundary marker does not match the mesh boundary attributes.");
			for (int i = 0; i < ess_bdr.Size(); ++i)
			{
				ess_bdr[i] = ess_bdr[i] || axisymmetric_mesh.axis_boundary[i];
			}
		}

		// Build the FE space and everything bound to it for the starting mesh.
		BuildOperators();

		BoundaryConditionValidator validator(mesh, *fespace);
		validator.ValidateBoundaryConditions(
			closure_bcs, /*terminals=*/{}, /*allow_overlap=*/false);

	}

	void BuildOperators() override {
		fespace = std::make_unique<mfem::FiniteElementSpace>(&mesh, fec.get());

		A = std::make_unique<mfem::GridFunction>(fespace.get());
		*A = 0.0;
		neumann_rhs = AssembleNeumannBoundaryLoad();

		a = std::make_unique<mfem::BilinearForm>(fespace.get());
		a->AddDomainIntegrator(MakeStiffnessIntegrator()); // a takes ownership
		a->Assemble();

		fespace->GetEssentialTrueDofs(ess_bdr, ess_tdof_list);

		// Form the constrained system operator. The eliminated-column part
		// (mat_e, used to build each scenario's RHS) is bound to A_op, which is
		// reused for every scenario's solve on this mesh.
		a->FormSystemMatrix(ess_tdof_list, A_op);
	}

	// Build the domain stiffness integrator matching the active geometry:
	// axisymmetric curl-curl (curl (nu curl A) with the 2*pi*r measure / A * nu/r^2 term) or
	// planar diffusion (nu * grad A * grad nu). A fresh instance is returned each call so the
	// solve (owned by 'a') and AMR error estimator can own separate copies.
	mfem::BilinearFormIntegrator* MakeStiffnessIntegrator() const {
		if (geometry == GeometryType::Axisymmetric)
		{
			auto* integ = new AxisymmetricCurlCurlIntegrator(*nu_coeff);
			return integ;
		}
		else
		{
			return new mfem::DiffusionIntegrator(*nu_coeff);
		}
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
	// per-region reluctivity discontinuities are respected. SetWithCoeff(true)
	// makes the flux nu*grad(A) (planar) / the curl-curl flux (axisym),
	// consistent with the integrator's energy norm.
	//
	// @param combined  Output: per-element max error indicator (sized to NE).
	// @return Global error sqrt(sum_k combined_k^2).
	double EstimateCombinedError(mfem::Vector& combined) override {
		const int sdim = mesh.SpaceDimension();
		std::unique_ptr<mfem::BilinearFormIntegrator> flux_integ(MakeStiffnessIntegrator());
		mfem::FiniteElementSpace flux_fes(&mesh, fec.get(), sdim);
		mfem::ZienkiewiczZhuEstimator estimator(*flux_integ, *A, flux_fes);
		estimator.SetWithCoeff(true);     // flux = nu * grad(A)
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

	// Peak flux density |B| over the current solution *A, sampled at element
	// nodes. AMR convergence diagnostic: the peak near a conductor corner or
	// high-permeability edge should settle as refinement resolves it.
	//
	// Planar: B = (dA/dy, -dA/dx), so |B| == |grad(A)| exactly.
	// Axisymmetric: B_r=-dA/dz, B_z=A/r+dA/dr, so the A/r term matters; reuse
	// MagneticFieldCoefficient (same B reconstruction as the exporters, incl.
	// the r->0 limit). Reflects whichever solution currently lives in *A.
	double ComputePeakFieldMagnitude() const override {
		if (!A) { return 0.0; }

		MagneticFieldCoefficient B_axi(A.get());
		const bool axi = (geometry == GeometryType::Axisymmetric);

		double peak = 0.0;
		mfem::Vector B;
		for (int e = 0; e < fespace->GetNE(); ++e) {
			const mfem::FiniteElement* fe = fespace->GetFE(e);
			mfem::ElementTransformation* T = fespace->GetElementTransformation(e);
			const mfem::IntegrationRule& nodes = fe->GetNodes();
			for (int i = 0; i < nodes.GetNPoints(); ++i) {
				const mfem::IntegrationPoint& ip = nodes.IntPoint(i);
				T->SetIntPoint(&ip);
				if (axi) { B_axi.Eval(B, *T, ip); }  // true |B| incl. A/r term
				else { A->GetGradient(*T, B); }   // |B| == |grad(A)| (planar)
				const double mag = B.Norml2();
				if (mag > peak) { peak = mag; }
			}
		}
		return peak;
	}

	void ImprintScenario(const Scenario& sc) {
		*A = 0.0; // Reset solution for new scenario

		// Re-apply non-zero essential (closure) BC values on this mesh's A.
		// ess_tdof values are lifted into the RHS by FormLinearSystem at solve time,
		// so they must be set AFTER the *A = 0.0 reset, every scenario.
		for (const auto& bc : config.BoundaryConditions) {
			if (bc.Type == BoundaryConditionType::Dirichlet && bc.Value != 0.0) {
				const EntityGroup& group = config.EntityGroups.at(bc.EntityGroupName);
				auto marker = MarkerFromAttrs(group.AttributeIds);
				mfem::ConstantCoefficient c(bc.Value);
				A->ProjectBdrCoefficient(c, marker);
			}
		}
		// Axis stays at A=0 (already zero from the reset; no projection needed).

		auto j_src = BuildCurrentDensity(sc);
		j_coeff = std::make_unique<mfem::PWConstCoefficient>(j_src);
		// RHS
		b = std::make_unique<mfem::LinearForm>(fespace.get());

		if (geometry == GeometryType::Axisymmetric)
		{
			// Integrates J * v * r  (global 2π omitted consistently)
			b->AddDomainIntegrator(new AxisymmetricLFIntegrator(*j_coeff));
		}
		else
		{
			b->AddDomainIntegrator(new mfem::DomainLFIntegrator(*j_coeff));
		}
		b->Assemble();
		*b += neumann_rhs;
	}

	// Solve + save on the CURRENT mesh/operators. Both analysis types flow through
	// ONE imprint -> solve -> save loop over BuildSolveScenarios() (authored
	// scenarios for Field; synthetic per-terminal unit-current drives for
	// CouplingMatrix). AMR calls this once more on the converged mesh so the
	// exported fields match the exported mesh.
	void RunOnCurrentMesh() override
	{
		if (config.AnalysisType == AnalysisType::CouplingMatrix) {
			L = std::make_unique<mfem::DenseMatrix>(config.Terminals.size());
			// CouplingMatrix synthesizes a unit-current scenario per terminal, so
			// every terminal must be current-driven for the drive to be meaningful.
			for (const auto& [term_name, term] : config.Terminals) {
				MFEM_VERIFY(term.ExcitationType == Quantity::Current,
					"CouplingMatrix terminal '" + term_name +
					"' must be a Current terminal for the magnetostatic solver.");

			}
		}

		for (const auto& [sc_name, sc] : BuildSolveScenarios()) {
			auto operation = Reporter().Start("scenario '" + sc_name + "'");
			ImprintScenario(sc);
			SolveSystem();
			SaveScenario(sc_name);
			if (config.AnalysisType == AnalysisType::CouplingMatrix) {
				// Each scenario is a unit-current drive on one terminal, so the
				// solution's flux linkage / inductance is the corresponding column
				// of the coupling matrix.
				const int col = std::distance(config.Terminals.begin(),
					config.Terminals.find(sc.Excitations[0].TerminalName));
				for (int row = 0; row < L->Height(); ++row) {
					const std::string& row_term = std::next(config.Terminals.begin(), row)->first;
					(*L)(row, col) = ComputeFluxLinkage(row_term);
				}
			}
		}
		if (config.AnalysisType == AnalysisType::CouplingMatrix) {
			WriteCouplingMatrix();
		}
	}

	void SolveSystem() {
		auto operation = Reporter().Start("linear system solve");
		// Form and solve
		mfem::OperatorPtr Aop;
		mfem::Vector X, B;

		a->FormLinearSystem(ess_tdof_list, *A, *b, Aop, X, B);

		if (B.Norml2() < 1e-12 && X.Norml2() < 1e-12)
		{
			mfem::out << "WARNING: Linear system RHS is ~zero. "
				<< "Check that 'sources' in JSON match mesh attributes.\n";
		}

		auto* sp = dynamic_cast<mfem::SparseMatrix*>(Aop.Ptr());
		MFEM_ASSERT(sp, "Expected SparseMatrix operator from FormLinearSystem.");

		mfem::GSSmoother M(*sp);
		mfem::PCG(*sp, M, B, X,
			Reporter().SolverPrintLevel(config.SolverPrintLevel),
			config.SolverMaxIter,
			config.SolverTolerance,
			0.0);

		a->RecoverFEMSolution(X, *b, *A);

		std::ostringstream statistics;
		statistics << "=== A Statistics ===\n"
			<< "  A min:     " << A->Min() << "\n"
			<< "  A max:     " << A->Max() << "\n"
			<< "  A L2 norm: " << A->Norml2();
		Reporter().Diagnostic(statistics.str());
	}

	// Post-solve field recovery: the vector potential A and the flux density B.
	// B = curl(A): the axisymmetric form (B_r = -dA/dz, B_z = dA/dr + A/r) needs
	// MagneticFieldCoefficient; the planar form is the standard curl of A_z.
	// Serialization is handled by the base class.
	FieldExportSet CollectExportFields() const override
	{
		FieldExportSet fields;
		fields.AddPrimaryScalar("A", *A);

		if (geometry == GeometryType::Axisymmetric) {
			fields.AddVector("B", std::make_unique<MagneticFieldCoefficient>(A.get()));
		}
		else {
			fields.AddVector("B", std::make_unique<PlanarMagneticFieldCoefficient>(A.get()));
		}
		return fields;
	}

	void SaveAnalysis() override
	{
		if (config.AnalysisType == AnalysisType::CouplingMatrix) {
			WriteCouplingMatrix();
		}
	}

	void WriteCouplingMatrix() {
		if (!L) {
			Reporter().Warning("WriteCouplingMatrix: coupling matrix not computed.");
			return;
		}

		SaveCouplingMatrix(*L, "Inductance Matrix [H]", "inductance_matrix.csv");
	}

private:

	// Flux linkage of terminal `terminal_name` for the solution currently in *A.
	//
	//   lambda_k = integral over terminal k of (N_k/area_k) * A dV
	//
	// i.e. the winding functional of the MEASURED terminal applied to the field,
	// NOT the source of the driving scenario. With a unit-current drive on
	// terminal i, lambda_k is directly L(k,i) - the mutual inductance for k != i.
	// Mirrors MagnetoquasistaticSolver::ComputeStrandedFluxLinkage.
	double ComputeFluxLinkage(const std::string& terminal_name) const
	{
		mfem::Vector unit_density = BuildTerminalCurrentDensity(terminal_name, 1.0);
		mfem::PWConstCoefficient unit_density_coeff(unit_density);

		mfem::LinearForm winding_functional(fespace.get());
		if (geometry == GeometryType::Axisymmetric) {
			winding_functional.AddDomainIntegrator(
				new AxisymmetricLFIntegrator(unit_density_coeff));
		}
		else {
			winding_functional.AddDomainIntegrator(
				new mfem::DomainLFIntegrator(unit_density_coeff));
		}
		winding_functional.Assemble();

		// Both integrators carry the full geometric measure, so this is webers.
		return winding_functional * *A;
	}

	// Uniform current density I/area over the terminal's domain attributes,
	// laid out per mesh attribute for a PWConstCoefficient.
	mfem::Vector BuildTerminalCurrentDensity(
		const std::string& terminal_name, double current) const {
		const Terminal& term = config.Terminals.at(terminal_name);
		const EntityGroup& group = config.EntityGroups.at(term.EntityGroupName);
		const double area = CalculateRegionArea(group.AttributeIds);
		MFEM_VERIFY(area > 0.0,
			"Current terminal '" + terminal_name + "' has zero cross-section.");

		mfem::Vector current_density(mesh.attributes.Max());
		current_density = 0.0;
		const double density = current / area;
		for (int attr : group.AttributeIds) {
			if (attr > 0 && attr <= current_density.Size()) {
				current_density[attr - 1] = density;
			}
		}
		return current_density;
	}

	mfem::Vector BuildCurrentDensity(const Scenario& sc) const {
		mfem::Vector j_src(mesh.attributes.Max());
		j_src = 0.0;

		for (const auto& [term_name, term] : config.Terminals) {
			if (term.ExcitationType == Quantity::Current) {
				double I = 0.0;
				for (const auto& exc : sc.Excitations) {
					if (exc.TerminalName == term_name) {
						I = exc.Value;
					}
				}

				if (I == 0.0) continue;
				j_src += BuildTerminalCurrentDensity(term_name, I);
			}
		}
		return j_src;
	}
};
