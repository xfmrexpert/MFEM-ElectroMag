// Copyright (c) 2026 T. C. Raymond
// SPDX-License-Identifier: MIT

#pragma once

#include <cctype>
#include <iomanip>
#include <limits>
#include <memory> // Required for smart pointers
#include <set>
#include <sstream>
#include "mfem.hpp"
#include "physics_solver.hpp"
#include "axisymmetric_curl_curl_integrator.hpp"
#include "axisymmetric_mass_integrator.hpp"
#include "axisymmetric_lf_integrator.hpp"
#include "magnetic_field_coefficient.hpp"
#include "complex_vector_magnitude_coefficient.hpp"
#include "constants.hpp"
#include "boundary_validation.hpp"
#include "problem_config.hpp"
#include "mqs_massive_port_operator.hpp"
#include "complex_block_layout.hpp"
#include "axisymmetric_conductance_coefficient.hpp"
#include "gmsh_results_writer.hpp"
#include "amr_support.hpp"

class MagnetoquasistaticSolver : public PhysicsSolver {

    double frequency = 0.0;
    mfem::real_t omega = 0.0;
    
    // Complex system objects. S_AA owns the real/imag field matrices referenced
    // by port_operator, so it is declared first and outlives that operator.
	std::unique_ptr<mfem::SesquilinearForm> S_AA;
	std::unique_ptr<mfem::ComplexGridFunction> A; // Complex field solution
	std::unique_ptr<mfem::Vector> Re_port_values; // Real part of port voltages
	std::unique_ptr<mfem::Vector> Im_port_values; // Imag part of port voltages
    std::unique_ptr<MqsMassivePortOperator> port_operator;
	std::unique_ptr<mfem::Vector> x_combined;
	std::unique_ptr<mfem::Vector> b_combined;

    // Coefficients
    std::unique_ptr<mfem::PWConstCoefficient> nu_coeff;
    std::unique_ptr<mfem::PWConstCoefficient> sigma_coeff;
    std::unique_ptr<mfem::PWConstCoefficient> j_coeff;     
    mfem::Vector conductivity_values;
    std::vector<mfem::real_t> port_conductances;
    
    mfem::Array<int> ess_mesh_tdofs;   // scalar-space essential mesh true DOFs
    std::unordered_map<std::string, mfem::Array<int>> terminal_markers; // Terminal name to boundary marker mapping

    struct CouplingResult {
        std::string Name;
        double Frequency;
        std::unique_ptr<mfem::DenseMatrix> Resistance;
        std::unique_ptr<mfem::DenseMatrix> Inductance;
    };
    struct MassivePortDefinition {
        std::string Name;
        std::vector<int> AttributeIds;
        double Conductivity;
    };
    std::vector<CouplingResult> coupling_results;
    mfem::DenseMatrix* resistance_matrix = nullptr;
    mfem::DenseMatrix* inductance_matrix = nullptr;

    // Material property pickers for MaterialVector, named instead of inlined as
    // lambdas so the Setup() coefficient construction reads at a glance.
    static double Reluctivity(const Material& m) {
        return 1.0 / (Constants::MU_0 * m.RelPermeability);
    }
    static double Conductivity(const Material& m) {
        return m.Conductivity;
    }

    // Function to build the port vector for a specific port attribute
    std::unique_ptr<mfem::Vector> BuildPortVector(mfem::FiniteElementSpace* fespace,
                            std::vector<int> port_attributes,
                            double sigma)
    {
        // Restrict integration to this specific port's attribute
        mfem::Array<int> port_marker(fespace->GetMesh()->attributes.Max());
        port_marker = 0;
        for (auto port_attribute : port_attributes) {
            port_marker[port_attribute - 1] = 1; // MFEM attributes are 1-indexed
        }

        // Constant coefficient restricted to this port (axisymmetric and planar
        // alike); referenced only while we assemble below.
        mfem::ConstantCoefficient coeff(sigma);
        mfem::RestrictedCoefficient restricted_coeff(coeff, port_marker);

        // Assemble the LinearForm using a scalar domain integrator
        mfem::LinearForm port_lf(fespace);
        port_lf.AddDomainIntegrator(new mfem::DomainLFIntegrator(restricted_coeff));
        port_lf.Assemble();

        // Extract and return as a standalone Vector
        auto port_vector = std::make_unique<mfem::Vector>(port_lf.Size());
        *port_vector = port_lf;
        return port_vector;
    }

    // Smallest physical radius attained by the elements carrying any of
    // @p attribute_ids, sampled through the element transformations so curved
    // geometry is respected.
    double MinRadiusOverAttributes(const std::vector<int>& attribute_ids) const
    {
        std::set<int> attrs(attribute_ids.begin(), attribute_ids.end());
        double min_r = std::numeric_limits<double>::max();
        mfem::Vector pos(mesh.SpaceDimension());

        for (int e = 0; e < mesh.GetNE(); ++e) {
            if (!attrs.count(mesh.GetAttribute(e))) { continue; }
            mfem::ElementTransformation* T = mesh.GetElementTransformation(e);
            const mfem::IntegrationRule& nodes = fespace->GetFE(e)->GetNodes();
            for (int i = 0; i < nodes.GetNPoints(); ++i) {
                const mfem::IntegrationPoint& ip = nodes.IntPoint(i);
                T->SetIntPoint(&ip);
                T->Transform(ip, pos);
                min_r = std::min(min_r, pos(0));
            }
        }
        return min_r;
    }

    // Function to compute G_dc for a specific port
    double ComputePortConductance(std::vector<int> port_attributes, double sigma)
    {
        // Setup an L2 space of order 0 for pure volumetric integration
        mfem::L2_FECollection l2_fec(0, mesh.Dimension());
        mfem::FiniteElementSpace l2_fes(&mesh, &l2_fec);

        // Create the restriction array for the port attribute
        mfem::Array<int> port_marker(mesh.attributes.Max());
        port_marker = 0;
        for (auto port_attribute : port_attributes) {
            port_marker[port_attribute - 1] = 1; // MFEM attributes are 1-indexed
        }

        // Define the appropriate coefficient
        mfem::Coefficient* base_coeff = nullptr;
        if (geometry == GeometryType::Axisymmetric)
        {
            base_coeff = new AxisymmetricConductanceCoeff(sigma);
        }
        else
        {
            base_coeff = new mfem::ConstantCoefficient(sigma); // For Planar, unit depth assumed? Or maybe per unit length.
            // If planar 2D, G usually per unit depth.
        }

        mfem::RestrictedCoefficient restricted_coeff(*base_coeff, port_marker);

        // Assemble the LinearForm to perform the spatial integration
        mfem::LinearForm g_form(&l2_fes);
        g_form.AddDomainIntegrator(new mfem::DomainLFIntegrator(restricted_coeff));
        g_form.Assemble();

        // The total integral is the sum of the piecewise constant values
        double G_dc = g_form.Sum();

        delete base_coeff;
        return G_dc;
    }

    // The complex block system is solved as a single real vector laid out
    // [Re_Mesh, Re_Port, Im_Mesh, Im_Port]; ComplexPortVectorView and
    // ConstComplexPortVectorView (complex_block_layout.hpp) name the four slots
    // so callers never compute packed indices by hand.

    // Lift the essential (closure) values currently projected into *A onto the
    // matching slots of the packed solution vector. FormLinearSystem constrains
    // essential DOFs to whatever it finds there, so without this the projected
    // non-zero closures would be forced back to zero.
    void LiftEssentialInto(mfem::Vector& x_packed) const {
        auto x = port_operator->View(x_packed);
        for (int k = 0; k < ess_mesh_tdofs.Size(); ++k) {
            const int d = ess_mesh_tdofs[k];
            x.ReMesh(d) = A->real()(d);
            x.ImMesh(d) = A->imag()(d);
        }
    }

	void RecoverSolvedUnknowns(mfem::Vector& x_packed) {
		auto x = port_operator->View(x_packed);
        // Copy field DOFs back into the complex grid function
        for (int i = 0; i < port_operator->Layout().NDofs(); ++i) {
            A->real()(i) = x.ReMesh(i);
            A->imag()(i) = x.ImMesh(i);
        }
        // Copy solved port voltages back into real/imaginary port vectors
		Re_port_values = std::make_unique<mfem::Vector>(port_operator->Layout().NPorts());
		Im_port_values = std::make_unique<mfem::Vector>(port_operator->Layout().NPorts());
		for (int p = 0; p < port_operator->Layout().NPorts(); ++p) {
			// Recover the solved port values from the packed vector
            (*Re_port_values)(p) = x.RePort(p);
            (*Im_port_values)(p) = x.ImPort(p);
		}
	}

public:
    // Constructor deals only with initialization, no manual nullptr assignment needed
    MagnetoquasistaticSolver(mfem::Mesh &m, const ProblemConfig &c) : PhysicsSolver(m, c) {}

    void Setup() override {
        int order = config.Order;
        const int dim = mesh.Dimension();

        MFEM_VERIFY(!config.Scenarios.empty(),
            "Magnetoquasistatic simulations require at least one frequency scenario.");
        frequency = config.Scenarios.front().second.Frequency;
        omega = Constants::TWO_PI * frequency;

        // Axisymmetric or Planar
        geometry = config.GeometryType;

        // Reject negative radii and record whether the domain reaches r = 0.
        ValidateAxisymmetricGeometry();

        // FE spaces
        fec = std::make_unique<mfem::H1_FECollection>(order, mesh.Dimension());
        
        // Materials
        // Real part: reluctivity nu = 1/mu.
        nu_coeff = MaterialCoefficient(1.0 / Constants::MU_0, Reluctivity);

        // Assemble conductivity without frequency scaling so the mass matrix can
        // be reused at every sweep point.
        conductivity_values = MaterialVector(0.0, Conductivity);
        sigma_coeff = std::make_unique<mfem::PWConstCoefficient>(conductivity_values);

        auto bcs = BuildClosureBcs();

        // All closures are essential (Dirichlet). Non-zero values are lifted into
        // x_combined in ImprintScenario() so FormLinearSystem constrains to them.
        std::vector<mfem::Array<int>> ess_markers;
        for (const auto& [marker, val] : bcs)
            ess_markers.push_back(marker);
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

        // Validate that BCs don't create physical conflicts
        BoundaryConditionValidator validator(mesh, *fespace);
        validator.ValidateBoundaryConditions(bcs, /*terminals=*/{}, false);  // Strict mode - reject conflicts
    }

    void BuildOperators() override {
		// Build the FE space and everything bound to it for the starting mesh.
		fespace = std::make_unique<mfem::FiniteElementSpace>(&mesh, fec.get());

        // Setup Complex Billinear Form
        S_AA = std::make_unique<mfem::SesquilinearForm>(fespace.get(), mfem::ComplexOperator::HERMITIAN);
        S_AA->AddDomainIntegrator(MakeStiffnessIntegrator(), nullptr);
        S_AA->AddDomainIntegrator(nullptr, MakeMassIntegrator());

        S_AA->Assemble();
        S_AA->Finalize();

        mfem::BilinearForm& K = S_AA->real();
        mfem::BilinearForm& M_sigma = S_AA->imag();

        const int n_dofs = fespace->GetTrueVSize();

        std::vector<MassivePortDefinition> massive_ports;
        massive_ports.reserve(config.Terminals.size() + config.Regions.size());
        for (const auto& [term_name, term] : config.Terminals) {
            if (term.Conductor != ConductorType::Massive) continue;
            const EntityGroup& group = config.EntityGroups.at(term.EntityGroupName);
            massive_ports.push_back(
                { term_name, group.AttributeIds, TerminalConductivity(term) });
        }
        for (const Region& region : config.Regions) {
            if (region.CurrentConstraint != RegionCurrentConstraint::Open) continue;
            const EntityGroup& group = config.EntityGroups.at(region.EntityGroupName);
            massive_ports.push_back({ region.EntityGroupName, group.AttributeIds,
                config.Materials[region.Material].Conductivity });
        }

        // Explicit massive terminals come first so their solved voltage indices
        // retain terminal order. Passive open-current regions follow and receive
        // an identically zero current RHS in every scenario.
        std::vector<std::unique_ptr<mfem::Vector>> port_loads;
        port_loads.reserve(massive_ports.size());
        port_conductances.clear();
        port_conductances.reserve(massive_ports.size());
        for (const MassivePortDefinition& port : massive_ports) {
            if (geometry == GeometryType::Axisymmetric) {
                // G_dc = integral sigma/(2*pi*r) diverges for a toroidal massive
                // conductor whose cross-section reaches the symmetry axis.
                MFEM_VERIFY(
                    MinRadiusOverAttributes(port.AttributeIds) > axisymmetric_mesh.tolerance,
                    "Massive port '" + port.Name + "' touches the symmetry axis. "
                    "Its DC conductance integral sigma/(2*pi*r) is divergent; "
                    "model it as a stranded conductor or move it off the axis.");
            }
            port_loads.push_back(
                BuildPortVector(fespace.get(), port.AttributeIds, port.Conductivity));
            double G_dc = ComputePortConductance(port.AttributeIds, port.Conductivity);
            MFEM_VERIFY(G_dc > 0.0,
                "Massive port '" + port.Name + "' has zero conductance.");
            port_conductances.push_back(G_dc);
        }

        // Hand the field matrices and port data to the block-system owner, which
        // assembles the complex saddle-point operator and owns all the wiring.
        port_operator = std::make_unique<MqsMassivePortOperator>(
            n_dofs, K.SpMat(), M_sigma.SpMat(), std::move(port_loads),
            port_conductances, omega);

        fespace->GetEssentialTrueDofs(ess_bdr, ess_mesh_tdofs);   // indices in [0, N_DOFs)

        // Each scalar essential DOF constrains both its real and imaginary copy
        // in the packed [Re|Im] layout (half-size = N_DOFs + N_Ports).
        ess_tdof_list = port_operator->MakeEssentialTDofs(ess_mesh_tdofs);

        // Grid Function (for solution recovery later)
        A = std::make_unique<mfem::ComplexGridFunction>(fespace.get());
	}

	mfem::BilinearFormIntegrator* MakeStiffnessIntegrator() {
		if (geometry == GeometryType::Axisymmetric) {
			auto* integ = new AxisymmetricCurlCurlIntegrator(*nu_coeff);
			return integ;
		}
		else {
			return new mfem::DiffusionIntegrator(*nu_coeff);
		}
	}

	mfem::BilinearFormIntegrator* MakeMassIntegrator() {
		if (geometry == GeometryType::Axisymmetric) {
            return new AxisymmetricMassIntegrator(*sigma_coeff);
		}
		else {
            return new mfem::MassIntegrator(*sigma_coeff);
		}
	}

    void ActivateFrequency(const Scenario& sc) {
        MFEM_VERIFY(std::isfinite(sc.Frequency) && sc.Frequency > 0.0,
            "MQS scenario frequency must be finite and positive.");
        frequency = sc.Frequency;
        omega = Constants::TWO_PI * frequency;
        port_operator->SetOmega(omega);
    }

    void ImprintScenario(const Scenario& sc) {
        ActivateFrequency(sc);
        *A = 0.0;

        // Re-apply non-zero essential (closure) BC values on this mesh's A.
        // ess_tdof values are lifted into the RHS by FormLinearSystem at solve time,
        // so they must be set AFTER the *A = 0.0 reset, every scenario.
        for (const auto& bc : config.BoundaryConditions) {
            if (bc.Value != 0.0) {
                const EntityGroup& group = config.EntityGroups.at(bc.EntityGroupName);
                auto marker = MarkerFromAttrs(group.AttributeIds);
                mfem::ConstantCoefficient c_re(bc.Value);
				mfem::ConstantCoefficient c_im(0.0);
                A->ProjectBdrCoefficient(c_re, c_im, marker);
            }
        }
        // Axis stays at A=0 (already zero from the reset; no projection needed).

        // The monolithic real/imag solver vectors (size 2*(N_DOFs+N_Ports)) are
        // laid out [Re_Mesh, Re_Port, Im_Mesh, Im_Port]; assemble the RHS
        // directly into that layout through a typed view rather than raw indices.
        b_combined = std::make_unique<mfem::Vector>(port_operator->Layout().FullSize());
        x_combined = std::make_unique<mfem::Vector>(port_operator->Layout().FullSize());
        *b_combined = 0.0;
        auto b = port_operator->View(*b_combined);

        // Source
        auto j_src = BuildCurrentDensity(sc);
        j_coeff = std::make_unique<mfem::PWConstCoefficient>(j_src);

        // Assemble the source term (J is assumed real) into the Re_Mesh block.
        mfem::LinearForm b_source(fespace.get());
        if (geometry == GeometryType::Axisymmetric) {
            b_source.AddDomainIntegrator(new AxisymmetricLFIntegrator(*j_coeff));
        }
        else {
            b_source.AddDomainIntegrator(new mfem::DomainLFIntegrator(*j_coeff));
        }
        b_source.Assemble();
        const mfem::real_t* b_source_data = b_source.GetData();   // bypass LinearForm::operator()
        for (int d = 0; d < port_operator->Layout().NDofs(); ++d) {
            b.ReMesh(d) += b_source_data[d];
        }

        // Drive the active port(s) via the imaginary port block Im_Port.
        int p = 0;
        for (const auto& [term_name, term] : config.Terminals) {
            if (term.Conductor != ConductorType::Massive) continue;   // keep p aligned
            for (const auto& exc : sc.Excitations)
                if (exc.TerminalName == term_name)
                    b.ImPort(p) = -exc.Value / omega;
            ++p;
        }

        *x_combined = 0.0; // Initial guess
        LiftEssentialInto(*x_combined);
    }

    // Solve + save on the CURRENT mesh/operators. Field analysis performs one
    // solve per concrete scenario. Coupling analysis uses each authored scenario
    // as a frequency point and synthesizes one unit-current solve per terminal.
    void RunOnCurrentMesh() override {
        PrepareAnalysis();

        if (config.AnalysisType == AnalysisType::Field) {
            for (const auto& [name, scenario] : config.Scenarios) {
                auto operation = Reporter().Start("scenario '" + name + "'");
                ImprintScenario(scenario);
                SolveSystem();
                SaveScenario(name);
            }
            return;
        }

        for (const auto& [point_name, point] : config.Scenarios) {
            BeginCouplingPoint(point_name, point.Frequency);
            for (const auto& [term_name, term] : config.Terminals) {
                Scenario column;
                column.Frequency = point.Frequency;
                column.Excitations.push_back({ term_name, 1.0 });
                auto operation = Reporter().Start(
                    "scenario '" + point_name + "', terminal '" + term_name + "'");
                ImprintScenario(column);
                SolveSystem();
                GatherCouplingColumn(column);
            }
        }
    }

    void SolveSystem() {
        auto operation = Reporter().Start("linear system solve");
        // Solve
        mfem::OperatorHandle A_op;
        mfem::Vector B_vec, X_vec;

        mfem::Operator* A_op_ptr;

        mfem::ComplexOperator& complex_system = port_operator->Operator();
        complex_system.FormLinearSystem(ess_tdof_list, *x_combined, *b_combined, A_op_ptr, X_vec, B_vec);
        bool own_A = (A_op_ptr != &complex_system);
        A_op.Reset(A_op_ptr, own_A);

        // Iterative Complex Solver
        mfem::GMRESSolver gmres;
        gmres.SetOperator(*A_op.Ptr());
        gmres.SetPrintLevel(Reporter().SolverPrintLevel(config.SolverPrintLevel));
        gmres.SetRelTol(config.SolverTolerance);
        gmres.SetMaxIter(config.SolverMaxIter);
		gmres.Mult(B_vec, X_vec);

		// X_vec is laid out [Re_Mesh, Re_Port, Im_Mesh, Im_Port]; copy the mesh
		// (field) DOFs back into the complex grid function, dropping the ports.
		RecoverSolvedUnknowns(X_vec);
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
    // per-region reluctivity discontinuities are respected. The real and
    // imaginary indicators are combined as a complex magnitude before the base
    // class folds them across scenarios.
    //
    // @param combined  Output: per-element max error indicator (sized to NE).
    // @return Global error sqrt(sum_k combined_k^2).
    double EstimateCombinedError(mfem::Vector& combined) override {
        const int sdim = mesh.SpaceDimension();
        std::unique_ptr<mfem::BilinearFormIntegrator> flux_integ(MakeStiffnessIntegrator());
        mfem::FiniteElementSpace flux_fes(&mesh, fec.get(), sdim);
        mfem::ZienkiewiczZhuEstimator estimator_re(*flux_integ, A->real(), flux_fes);
        estimator_re.SetWithCoeff(false);     // flux = nu * grad(A)
        estimator_re.SetFluxAveraging(1);    // do not average across attribute interfaces

        mfem::ZienkiewiczZhuEstimator estimator_im(*flux_integ, A->imag(), flux_fes);
        estimator_im.SetWithCoeff(false);     // flux = nu * grad(A)
        estimator_im.SetFluxAveraging(1);    // do not average across attribute interfaces

        const auto amr_scenarios = BuildAmrSolveScenarios();
        return EstimateScenarioMaximumErrorOver(amr_scenarios, combined,
            [this](const Scenario& scenario) {
                ImprintScenario(scenario);
                SolveSystem();
            },
            [&estimator_re, &estimator_im](mfem::Vector& current) {
                estimator_re.Reset();
                const mfem::Vector& errs_re = estimator_re.GetLocalErrors();
                estimator_im.Reset();
                const mfem::Vector& errs_im = estimator_im.GetLocalErrors();
                current.SetSize(errs_re.Size());
                for (int k = 0; k < current.Size(); ++k) {
                    current(k) = std::hypot(errs_re(k), errs_im(k));
                }
            });
    }

    std::vector<std::pair<std::string, Scenario>> BuildAmrSolveScenarios() const {
        if (config.AnalysisType == AnalysisType::Field) {
            return config.Scenarios;
        }

        std::vector<std::pair<std::string, Scenario>> scenarios;
        scenarios.reserve(config.Scenarios.size() * config.Terminals.size());
        for (const auto& [point_name, point] : config.Scenarios) {
            for (const auto& [term_name, term] : config.Terminals) {
                Scenario column;
                column.Frequency = point.Frequency;
                column.Excitations.push_back({ term_name, 1.0 });
                scenarios.emplace_back(point_name + "_" + term_name, std::move(column));
            }
        }
        return scenarios;
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

        MagneticFieldCoefficient B_axi_re(&A->real());
        MagneticFieldCoefficient B_axi_im(&A->imag());
        const bool axi = (geometry == GeometryType::Axisymmetric);

        double peak = 0.0;
        mfem::Vector B_re;
        mfem::Vector B_im;
        for (int e = 0; e < fespace->GetNE(); ++e) {
            const mfem::FiniteElement* fe = fespace->GetFE(e);
            mfem::ElementTransformation* T = fespace->GetElementTransformation(e);
            const mfem::IntegrationRule& nodes = fe->GetNodes();
            for (int i = 0; i < nodes.GetNPoints(); ++i) {
                const mfem::IntegrationPoint& ip = nodes.IntPoint(i);
                T->SetIntPoint(&ip);
                if (axi) {
                    B_axi_re.Eval(B_re, *T, ip);
                    B_axi_im.Eval(B_im, *T, ip);
                }  // true |B| incl. A/r term
                else {
                    A->real().GetGradient(*T, B_re);
                    A->imag().GetGradient(*T, B_im);
                }   // |B| == |grad(A)| (planar)
                const double mag_re = B_re.Norml2();
                const double mag_im = B_im.Norml2();
                const double mag = std::hypot(mag_re, mag_im);
                if (mag > peak) { peak = mag; }
            }
        }
        return peak;
    }


	// Post-solve field recovery for the complex solution: the real/imaginary
	// vector-potential parts (primaries), the real/imaginary flux densities
	// B = curl(A) (derived vectors, geometry-dependent), and the complex flux
	// magnitude |B| = sqrt(|Re B|^2 + |Im B|^2). Serialization is handled by the
	// base class.
	FieldExportSet CollectExportFields() const override {
		FieldExportSet fields;
		fields.AddPrimaryScalar("A_Real", A->real());
		fields.AddPrimaryScalar("A_Imag", A->imag());

		mfem::VectorCoefficient* b_re;
		mfem::VectorCoefficient* b_im;
		if (geometry == GeometryType::Axisymmetric) {
			// Axisymmetric B = Curl(A_phi) = (-dA/dz, 1/r*d(rA)/dr)
			b_re = &fields.AddVector("B_Real",
                std::make_unique<MagneticFieldCoefficient>(&A->real()));
			b_im = &fields.AddVector("B_Imag",
                std::make_unique<MagneticFieldCoefficient>(&A->imag()));
		}
		else {
			// Planar B = Curl(A_z) = (dA/dy, -dA/dx)
			b_re = &fields.AddVector("B_Real",
                std::make_unique<PlanarMagneticFieldCoefficient>(&A->real()));
			b_im = &fields.AddVector("B_Imag",
                std::make_unique<PlanarMagneticFieldCoefficient>(&A->imag()));
		}

		fields.AddScalar("B_Magnitude",
			std::make_unique<ComplexVectorMagnitudeCoefficient>(*b_re, *b_im));
		return fields;
	}

	void GatherCouplingColumn(const Scenario& sc) {
        MFEM_VERIFY(sc.Excitations.size() == 1,
            "MQS coupling scenarios must drive exactly one terminal.");
        const auto driven = config.Terminals.find(sc.Excitations.front().TerminalName);
        MFEM_VERIFY(driven != config.Terminals.end(),
            "MQS coupling scenario references an unknown terminal.");
        const int column = static_cast<int>(std::distance(config.Terminals.begin(), driven));

        int massive_port = 0;
        int row = 0;
        for (const auto& [term_name, term] : config.Terminals) {
            if (term.Conductor == ConductorType::Massive) {
                (*resistance_matrix)(row, column) = (*Re_port_values)(massive_port);
                (*inductance_matrix)(row, column) =
                    (*Im_port_values)(massive_port) / omega;
                ++massive_port;
            }
            else {
                const auto [flux_re, flux_im] =
                    ComputeStrandedFluxLinkage(term_name);
                (*resistance_matrix)(row, column) = -omega * flux_im;
                (*inductance_matrix)(row, column) = flux_re;
            }
            ++row;
        }
	}

	void PrepareAnalysis() {
		if (config.AnalysisType == AnalysisType::CouplingMatrix) {
            for (const auto& [term_name, term] : config.Terminals) {
                MFEM_VERIFY(term.Excitation == Quantity::Current,
                    "MQS CouplingMatrix terminal '" + term_name +
                    "' must be a Current terminal.");
            }

            coupling_results.clear();
            coupling_results.reserve(config.Scenarios.size());
            resistance_matrix = nullptr;
            inductance_matrix = nullptr;
		}
	}

    void BeginCouplingPoint(const std::string& name, double point_frequency) {
        const int num_terminals = static_cast<int>(config.Terminals.size());
        CouplingResult result;
        result.Name = name;
        result.Frequency = point_frequency;
        result.Resistance = std::make_unique<mfem::DenseMatrix>(num_terminals, num_terminals);
        result.Inductance = std::make_unique<mfem::DenseMatrix>(num_terminals, num_terminals);
        *result.Resistance = 0.0;
        *result.Inductance = 0.0;
        coupling_results.push_back(std::move(result));
        resistance_matrix = coupling_results.back().Resistance.get();
        inductance_matrix = coupling_results.back().Inductance.get();
    }

    void SaveAnalysis() override
	{
		if (config.AnalysisType == AnalysisType::CouplingMatrix) {
			WriteCouplingMatrix();
		}
	}

	void WriteCouplingMatrix() {
        if (coupling_results.empty()) {
            Reporter().Warning("WriteCouplingMatrix: MQS coupling matrices not computed.");
            return;
        }

        for (const CouplingResult& result : coupling_results) {
            const std::string frequency_label = FrequencyOutputToken(result.Frequency) + "Hz";
            std::string output_tag = SafeOutputToken(result.Name);
            if (output_tag.size() < 2 ||
                output_tag.compare(output_tag.size() - 2, 2, "Hz") != 0) {
                output_tag += "_" + frequency_label;
            }
            SaveCouplingMatrix(*result.Inductance,
                "Inductance Matrix at " + frequency_label + " [H]",
                "inductance_matrix_" + output_tag + ".csv");
            SaveCouplingMatrix(*result.Resistance,
                "Resistance Matrix at " + frequency_label + " [Ohm]",
                "resistance_matrix_" + output_tag + ".csv");
        }
	}

    static std::string SafeOutputToken(std::string value) {
        for (char& c : value) {
            const unsigned char uc = static_cast<unsigned char>(c);
            if (!std::isalnum(uc) && c != '-' && c != '_') c = '_';
        }
        return value;
    }

    static std::string FrequencyOutputToken(double value) {
        std::ostringstream stream;
        stream << std::setprecision(12) << value;
        return SafeOutputToken(stream.str());
    }

	double TerminalConductivity(const Terminal& term) const {
		// First domain attribute of the terminal's group that a region claims
		// wins (preserves the prior "first match" behavior); 0 => non-conductive.
		const EntityGroup& group = config.EntityGroups.at(term.EntityGroupName);
		for (int attr : group.AttributeIds) {
			if (const Material* mat = MaterialForAttr(attr))
				return mat->Conductivity;
		}
		return 0.0;
	}

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

    std::pair<double, double> ComputeStrandedFluxLinkage(
        const std::string& terminal_name) const {
        mfem::Vector unit_density =
            BuildTerminalCurrentDensity(terminal_name, 1.0);
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

        const double geometry_scale =
            geometry == GeometryType::Axisymmetric ? Constants::TWO_PI : 1.0;
        return {
            geometry_scale * (winding_functional * A->real()),
            geometry_scale * (winding_functional * A->imag())
        };
    }

    mfem::Vector BuildCurrentDensity(const Scenario& sc) const {
        mfem::Vector j_src(mesh.attributes.Max());
        j_src = 0.0;

        for (const auto& [term_name, term] : config.Terminals) {
            if (term.Excitation == Quantity::Current && term.Conductor == ConductorType::Stranded) {
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