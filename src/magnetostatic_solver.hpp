// Copyright (c) 2026 T. C. Raymond
// SPDX-License-Identifier: MIT

#pragma once

#include <memory>
#include <cmath>
#include <algorithm>
#include <iomanip>
#include <sstream>

#include "mfem.hpp"
#include "physics_solver.hpp"
#include "axisymmetric_curl_curl_integrator.hpp"
#include "axisymmetric_lf_integrator.hpp"
#include "magnetic_field_coefficient.hpp"
#include "input_parser.hpp"
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

	// Cached constrained system for the CURRENT mesh. The matrix is identical
	// for every solve on a given mesh (same bilinear form and essential DOFs),
	// so it is assembled once per mesh in BuildOperators() and reused for all
	// scenarios / coupling columns. AMR refinement rebuilds it via
	// BuildOperators().
	mfem::OperatorPtr A_op;

	std::unique_ptr<mfem::DenseMatrix> L; // Inductance matrix (coupling matrix) for the current mesh

public:
	// One AMR iteration's diagnostics, recorded during RunAdaptive(). Exposed for
	// tests (convergence / conformity assertions) and console logging.
	struct AmrIterationInfo {
		long   true_dofs;     // global true DOFs on that iteration's mesh
		double global_error;  // sqrt(sum_k eta_k^2), combined over scenarios
		double peak_absB;     // max |B| sampled over the mesh
	};
private:
	std::vector<AmrIterationInfo> amr_history;

public:
	MagnetostaticSolver(mfem::Mesh& m, const json& c) : PhysicsSolver(m, c) {}

	void Setup() override
	{
		// Config & solver type
		InputParser parser(config_json);
		config = parser.GetProblemConfig();

		int order = config.Order;
		const int dim = mesh.Dimension();

		// Axisymmetric or Planar
		geometry = config.GeometryType;

		// FE collection
		fec = std::make_unique<mfem::H1_FECollection>(order, dim);

		// Material Properties (Reluctivity nu = 1/mu), keyed by mesh DOMAIN attribute.
		nu_coeff = MaterialCoefficient(1.0 / Constants::MU_0, [](const Material& m) {
			return 1.0 / (Constants::MU_0 * m.RelPermeability); });

		// All closures essential; values lifted from *A by FormLinearSystem.
		auto bcs = BuildClosureBcs();

		std::vector<mfem::Array<int>> ess_markers;
		for (const auto& [marker, val] : bcs) ess_markers.push_back(marker);
		ess_bdr = EssentialBdrFrom(ess_markers);

		// Axis regularity: enforce A_phi = 0 on r=0 as ESSENTIAL.
		// Best practice: mark the axis as an essential boundary via boundary attributes if your mesh has it tagged.
		// If you *don't* have the axis tagged as a boundary attribute, do a geometric fallback:
		if (geometry == GeometryType::Axisymmetric)
		{
			// Geometric fallback: force A=0 on axis boundary vertices by marking the boundary attributes
			// that lie on r=0. This requires detecting boundary elements on the axis and marking their attribute.
			// If your mesh already has an "axis" boundary attribute, prefer using that instead.
			MarkAxisBoundaryAttributesGeometric();
		}

		// Build the FE space and everything bound to it for the starting mesh.
		BuildOperators();

		BoundaryConditionValidator validator(mesh, *fespace);
		validator.ValidateBoundaryConditions(bcs, /*terminals=*/{}, /*allow_overlap=*/false);

	}

	void BuildOperators() {
		fespace = std::make_unique<mfem::FiniteElementSpace>(&mesh, fec.get());

		A = std::make_unique<mfem::GridFunction>(fespace.get());
		*A = 0.0;

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
			return new AxisymmetricCurlCurlIntegrator(*nu_coeff);
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
	double EstimateCombinedError(mfem::Vector& combined) {
		const int ne = mesh.GetNE();
		combined.SetSize(ne);
		combined = 0.0;

		const int sdim = mesh.SpaceDimension();
		std::unique_ptr<mfem::BilinearFormIntegrator> flux_integ(MakeStiffnessIntegrator());
		mfem::FiniteElementSpace flux_fes(&mesh, fec.get(), sdim);
		mfem::ZienkiewiczZhuEstimator estimator(*flux_integ, *A, flux_fes);
		estimator.SetWithCoeff(true);     // flux = nu * grad(A)
		estimator.SetFluxAveraging(1);    // do not average across attribute interfaces

		auto fold_current_solution = [&]() {
			estimator.Reset(); // force recompute: same mesh, new solution in *A
			const mfem::Vector& errs = estimator.GetLocalErrors();
			for (int k = 0; k < ne; ++k) {
				if (errs(k) > combined(k)) { combined(k) = errs(k); }
			}
			};

		// Same scenario set as the real solves (authored scenarios for Field, or
		// the synthetic per-terminal unit drives for CouplingMatrix), so the AMR
		// indicator reflects exactly what will be exported.
		for (const auto& [sc_name, sc] : BuildSolveScenarios()) {
			ImprintScenario(sc);
			SolveSystem();
			fold_current_solution();
		}

		double sum_sq = 0.0;
		for (int k = 0; k < ne; ++k) { sum_sq += combined(k) * combined(k); }
		return std::sqrt(sum_sq);
	}

	// Peak flux density |B| over the current solution *A, sampled at element
	// nodes. AMR convergence diagnostic: the peak near a conductor corner or
	// high-permeability edge should settle as refinement resolves it.
	//
	// Planar: B = (dA/dy, -dA/dx), so |B| == |grad(A)| exactly.
	// Axisymmetric: B_r=-dA/dz, B_z=A/r+dA/dr, so the A/r term matters; reuse
	// MagneticFieldCoefficient (same B reconstruction as the exporters, incl.
	// the r->0 limit). Reflects whichever solution currently lives in *A.
	double ComputePeakFieldMagnitude() const {
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

	// AMR per-iteration diagnostics from the most recent RunAdaptive(). Empty when
	// AMR is disabled. Consumed by the regression tests and useful for logging.
	const std::vector<AmrIterationInfo>& GetAmrHistory() const { return amr_history; }

	void ImprintScenario(const Scenario& sc) {
		*A = 0.0; // Reset solution for new scenario

		// Re-apply non-zero essential (closure) BC values on this mesh's A.
		// ess_tdof values are lifted into the RHS by FormLinearSystem at solve time,
		// so they must be set AFTER the *A = 0.0 reset, every scenario.
		for (const auto& bc : config.BoundaryConditions) {
			if (bc.Value != 0.0) {
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
	}

	void Run() override {
		if (config.Amr.Enabled) {
			RunAdaptive();
		}
		else {
			RunFixed();
		}
	}

	// Solve + save on the CURRENT mesh/operators. Both analysis types flow through
	// ONE imprint -> solve -> save loop over BuildSolveScenarios() (authored
	// scenarios for Field; synthetic per-terminal unit-current drives for
	// CouplingMatrix). AMR calls this once more on the converged mesh so the
	// exported fields match the exported mesh.
	void RunFixed()
	{
		if (config.AnalysisType == AnalysisType::CouplingMatrix) {
			L = std::make_unique<mfem::DenseMatrix>(config.Terminals.size());
			// CouplingMatrix synthesizes a unit-current scenario per terminal, so
			// every terminal must be current-driven for the drive to be meaningful.
			for (const auto& [term_name, term] : config.Terminals) {
				MFEM_VERIFY(term.Excitation == Quantity::Current,
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
					(*L)(row, col) = ComputeFluxLinkage(row_term, sc);
				}
			}
		}
		if (config.AnalysisType == AnalysisType::CouplingMatrix) {
			WriteCouplingMatrix();
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
			const double peak_absB = ComputePeakFieldMagnitude();

			amr_history.push_back({ cdofs, global_err, peak_absB });

			std::ostringstream diagnostic;
			diagnostic << "AMR iteration " << it
				<< ": elements=" << mesh.GetNE()
				<< ", true_dofs=" << cdofs
				<< ", global_error=" << std::scientific << std::setprecision(6) << global_err
				<< ", peak|B|=" << peak_absB;
			Reporter().Diagnostic(diagnostic.str());

			// Stopping criteria (any one stops): error tolerance, DOF budget, or
			// this being the last permitted iteration.
			if (amr.ErrorTolerance > 0.0 && global_err < amr.ErrorTolerance) {
				Reporter().Diagnostic("AMR: global error below tolerance. Stop.");
				break;
			}
			if (amr.MaxDofs > 0 && cdofs > amr.MaxDofs) {
				Reporter().Diagnostic("AMR: reached the maximum number of DOFs. Stop.");
				break;
			}
			if (it + 1 >= max_it) {
				Reporter().Diagnostic("AMR: reached the maximum number of iterations. Stop.");
				break;
			}

			// Mark and refine conformingly (throws if the mesh cannot refine
			// without hanging nodes), then rebuild the FE space / operators.
			mfem::Array<int> marked;
			amr::MarkElementsDorfler(errors, amr.ErrorFraction, marked);
			if (marked.Size() == 0) {
				Reporter().Diagnostic("AMR: no elements marked for refinement. Stop.");
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
			std::cerr << "WriteCouplingMatrix: coupling matrix not computed.\n";
			return;
		}

		SaveCouplingMatrix(*L, "Inductance Matrix [H]", "inductance_matrix.csv");
	}

private:

	double ComputeFluxLinkage(const std::string& terminal_name, const Scenario& sc) const
	{
		const mfem::Vector unit_source = BuildCurrentDensity(sc);

		const double geometry_scale =
			geometry == GeometryType::Axisymmetric
			? Constants::TWO_PI
			: 1.0;

		return geometry_scale * (unit_source * *A);
	}

	mfem::Vector BuildCurrentDensity(const Scenario& sc) const {
		mfem::Vector j_src(mesh.attributes.Max());
		j_src = 0.0;

		for (const auto& [term_name, term] : config.Terminals) {
			if (term.Excitation == Quantity::Current) {
				double I = 0.0;
				for (const auto& exc : sc.Excitations) {
					if (exc.TerminalName == term_name) {
						I = exc.Value;
					}
				}

				if (I == 0.0) continue;
				const std::string& group_name = term.EntityGroupName;
				const EntityGroup& group = config.EntityGroups.at(group_name);
				const double A = CalculateRegionArea(group.AttributeIds);
				MFEM_VERIFY(A > 0.0, "Current terminal '" + term_name + "' has zero cross-section.");
				const double J = I / A; // Current density = current / area

				for (int attr : group.AttributeIds) {
					if (attr > 0 && attr <= mesh.attributes.Max()) {
						j_src[attr - 1] = J;
					}
				}
			}
		}
		return j_src;
	}

	void AssembleCurrentLoad(const Scenario& sc)
	{

	}

	void AssembleCurrentLoad(const std::string& terminal_name)
	{
		
	}
};
