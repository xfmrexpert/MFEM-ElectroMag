// Copyright (c) 2026 T. C. Raymond
// SPDX-License-Identifier: MIT

#pragma once

#include <memory> // Required for smart pointers
#include "mfem.hpp"
#include "physics_solver.hpp"
#include "axisymmetric_curl_curl_integrator.hpp"
#include "axisymmetric_mass_integrator.hpp"
#include "axisymmetric_lf_integrator.hpp"
#include "magnetic_field_coefficient.hpp"
#include "complex_vector_magnitude_coefficient.hpp"
#include "input_parser.hpp"
#include "constants.hpp"
#include "boundary_validation.hpp"
#include "problem_config.hpp"
#include "port_coupled_system.hpp"
#include "complex_block_layout.hpp"
#include "axisymmetric_conductance_coefficient.hpp"
#include "gmsh_results_writer.hpp"

class MagnetoquasistaticSolver : public PhysicsSolver {
    double frequency = 60.0;
    mfem::real_t omega = Constants::TWO_PI * frequency;
    
	// Complex System objects. S_AA owns the real/imag field matrices that the
	// port-coupled block system references, so it is declared before 'system'
	// and outlives it (members destroy in reverse declaration order).
	std::unique_ptr<mfem::SesquilinearForm> S_AA;
	std::unique_ptr<mfem::ComplexGridFunction> A;
	std::unique_ptr<PortCoupledComplexSystem> system;
	std::unique_ptr<mfem::Vector> x_combined;
	std::unique_ptr<mfem::Vector> b_combined;

    int N_DOFs;
    int N_Ports;

    // Coefficients
    std::unique_ptr<mfem::PWConstCoefficient> nu_coeff;
    std::unique_ptr<mfem::PWConstCoefficient> omega_sigma_coeff;
    std::unique_ptr<mfem::PWConstCoefficient> j_coeff;     
    
    mfem::Array<int> ess_mesh_tdofs;   // scalar-space essential mesh true DOFs, in [0, N_DOFs)
    std::unordered_map<std::string, mfem::Array<int>> terminal_markers; // Terminal name to boundary marker mapping

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
        ComplexPortVectorView x(x_packed, N_DOFs, N_Ports);
        for (int k = 0; k < ess_mesh_tdofs.Size(); ++k) {
            const int d = ess_mesh_tdofs[k];
            x.ReMesh(d) = A->real()(d);
            x.ImMesh(d) = A->imag()(d);
        }
    }

    // Copy the mesh (field) DOFs out of the solved monolithic vector back into
    // the complex grid function, discarding the port unknowns.
    void UnpackComplexSolution(const mfem::Vector& x_packed) {
        ConstComplexPortVectorView x(x_packed, N_DOFs, N_Ports);
        for (int i = 0; i < N_DOFs; ++i) {
            A->real()(i) = x.ReMesh(i);
            A->imag()(i) = x.ImMesh(i);
        }
    }

public:
    // Constructor deals only with initialization, no manual nullptr assignment needed
    MagnetoquasistaticSolver(mfem::Mesh &m, const json &c) : PhysicsSolver(m, c) {}

    void Setup() override {
        // Config & solver type
        InputParser parser(config_json);
        config = parser.GetProblemConfig();
        
        int order = config.Order;
        const int dim = mesh.Dimension();

        frequency = config.Frequency;
        omega = Constants::TWO_PI * frequency;

        // Axisymmetric or Planar
        geometry = config.GeometryType;

        // FE spaces
        fec = std::make_unique<mfem::H1_FECollection>(order, mesh.Dimension());
        
        // Materials
        // Real part: reluctivity nu = 1/mu.
        nu_coeff = MaterialCoefficient(1.0 / Constants::MU_0, Reluctivity);

        // Imag part: omega * sigma (raw conductivity scaled by omega before wrapping).
        mfem::Vector omega_sigma = MaterialVector(0.0, Conductivity);
        omega_sigma *= omega;
        omega_sigma_coeff = std::make_unique<mfem::PWConstCoefficient>(omega_sigma);

        auto bcs = BuildClosureBcs();

        // All closures are essential (Dirichlet). Non-zero values are lifted into
        // x_combined in ImprintScenario() so FormLinearSystem constrains to them.
        std::vector<mfem::Array<int>> ess_markers;
        for (const auto& [marker, val] : bcs)
            ess_markers.push_back(marker);
        ess_bdr = EssentialBdrFrom(ess_markers);

        // Axis regularity: enforce A_phi = 0 on r=0 as ESSENTIAL.
        // Best practice: mark the axis as an essential boundary via boundary attributes if your mesh has it tagged.
        // If you *don't* have the axis tagged as a boundary attribute, do a geometric fallback:
        if (geometry == GeometryType::Axisymmetric)
        {
            // Geometric fallback: force A=0 on axis boundary vertices by marking the boundary attributes
            // that lie on r=0. This requires detecting boundary elements on the axis and marking their attribute.
            // If your mesh already has an "axis" boundary attribute, prefer using that in InputParser instead.
            MarkAxisBoundaryAttributesGeometric();
        }

        // Build the FE space and everything bound to it for the starting mesh.
        BuildOperators();

        // Validate that BCs don't create physical conflicts
        BoundaryConditionValidator validator(mesh, *fespace);
        validator.ValidateBoundaryConditions(bcs, /*terminals=*/{}, false);  // Strict mode - reject conflicts
    }

	void BuildOperators() {
		// Build the FE space and everything bound to it for the starting mesh.
		fespace = std::make_unique<mfem::FiniteElementSpace>(&mesh, fec.get());

        // Setup Complex Billinear Form
        S_AA = std::make_unique<mfem::SesquilinearForm>(fespace.get(), mfem::ComplexOperator::HERMITIAN);

        if (geometry == GeometryType::Axisymmetric) {
            // Real Part: Curl-Curl (1/mu) with r weight
            S_AA->AddDomainIntegrator(new AxisymmetricCurlCurlIntegrator(*nu_coeff), nullptr);

            // Imag Part: Mass (sigma * omega) with r weight
            S_AA->AddDomainIntegrator(nullptr, new AxisymmetricMassIntegrator(*omega_sigma_coeff));
        }
        else {
            // Planar 2D
            // Real Part: Diffusion (similar to Magnetostatic Planar)
            // Solves Div(nu Grad A)
            S_AA->AddDomainIntegrator(new mfem::DiffusionIntegrator(*nu_coeff), nullptr);

            // Imag Part: Standard Mass (sigma * omega)
            S_AA->AddDomainIntegrator(nullptr, new mfem::MassIntegrator(*omega_sigma_coeff));
        }

        S_AA->Assemble();
        S_AA->Finalize();

        mfem::BilinearForm& K = S_AA->real();
        mfem::BilinearForm& M_sigma = S_AA->imag();

        N_DOFs = fespace->GetTrueVSize();

        // N_Ports = massive terminals only
        N_Ports = 0;
        for (const auto& [name, term] : config.Terminals)
            if (term.Conductor == ConductorType::Massive) ++N_Ports;

        // One load vector + self-admittance diagonal entry per massive port.
        std::vector<std::unique_ptr<mfem::Vector>> port_loads;
        std::vector<mfem::real_t> g_scaled_diag;
        port_loads.reserve(N_Ports);
        g_scaled_diag.reserve(N_Ports);
        for (const auto& [term_name, term] : config.Terminals) {
            if (term.Conductor != ConductorType::Massive) continue;   // skip stranded
            const EntityGroup& group = config.EntityGroups.at(term.EntityGroupName);
            std::vector<int> attribute_ids = group.AttributeIds;
            double conductivity = TerminalConductivity(term);
            port_loads.push_back(BuildPortVector(fespace.get(), attribute_ids, conductivity));
            double G_dc = ComputePortConductance(attribute_ids, conductivity);
            MFEM_VERIFY(G_dc > 0.0, "Massive port '" + term_name + "' has zero conductance.");
            g_scaled_diag.push_back(-1.0 / (omega * G_dc));
        }

        // Hand the field matrices and port data to the block-system owner, which
        // assembles the complex saddle-point operator and owns all the wiring.
        system = std::make_unique<PortCoupledComplexSystem>(N_DOFs, K.SpMat(), M_sigma.SpMat(), std::move(port_loads), g_scaled_diag);

        fespace->GetEssentialTrueDofs(ess_bdr, ess_mesh_tdofs);   // indices in [0, N_DOFs)

        // Each scalar essential DOF constrains both its real and imaginary copy
        // in the packed [Re|Im] layout (half-size = N_DOFs + N_Ports).
        ess_tdof_list = ComplexEssentialTDofs(ess_mesh_tdofs, N_DOFs + N_Ports);

        // Grid Function (for solution recovery later)
        A = std::make_unique<mfem::ComplexGridFunction>(fespace.get());
	}

    void ImprintScenario(const Scenario& sc) {
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
        b_combined = std::make_unique<mfem::Vector>(system->FullSize());
        x_combined = std::make_unique<mfem::Vector>(system->FullSize());
        *b_combined = 0.0;
        ComplexPortVectorView b(*b_combined, N_DOFs, N_Ports);

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
        for (int d = 0; d < N_DOFs; ++d) { b.ReMesh(d) += b_source_data[d]; }

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

    // Solve + save on the CURRENT mesh/operators. Both analysis types flow through
    // ONE imprint -> solve -> save loop over BuildSolveScenarios() (authored
    // scenarios for Field; synthetic per-terminal unit-current drives for
    // CouplingMatrix). MQS has no AMR path, so this is the whole run.
    void Run() override {
        if (config.AnalysisType == AnalysisType::CouplingMatrix) {
            // CouplingMatrix synthesizes a unit-current scenario per terminal, so
            // every terminal must be current-driven for the drive to be meaningful.
            for (const auto& [term_name, term] : config.Terminals) {
                MFEM_VERIFY(term.Excitation == Quantity::Current,
                    "CouplingMatrix terminal '" + term_name +
                    "' must be a Current terminal for the magnetoquasistatic solver.");
            }
        }

        for (const auto& [sc_name, sc] : BuildSolveScenarios()) {
            ImprintScenario(sc);
            SolveSystem();
            SaveScenario(sc_name);
        }
    }

    void SolveSystem() {
        // Solve
        mfem::OperatorHandle A_op;
        mfem::Vector B_vec, X_vec;

        mfem::Operator* A_op_ptr;

        mfem::ComplexOperator& complex_system = system->GetOperator();
        complex_system.FormLinearSystem(ess_tdof_list, *x_combined, *b_combined, A_op_ptr, X_vec, B_vec);
        bool own_A = (A_op_ptr != &complex_system);
        A_op.Reset(A_op_ptr, own_A);

        // Iterative Complex Solver
        mfem::GMRESSolver gmres;
        gmres.SetOperator(*A_op.Ptr());
        gmres.SetPrintLevel(config.SolverPrintLevel);
        gmres.SetRelTol(config.SolverTolerance);
        gmres.SetMaxIter(config.SolverMaxIter);
		gmres.Mult(B_vec, X_vec);

		// X_vec is laid out [Re_Mesh, Re_Port, Im_Mesh, Im_Port]; copy the mesh
		// (field) DOFs back into the complex grid function, dropping the ports.
		UnpackComplexSolution(X_vec);
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
				std::make_unique<mfem::CurlGridFunctionCoefficient>(&A->real()));
			b_im = &fields.AddVector("B_Imag",
				std::make_unique<mfem::CurlGridFunctionCoefficient>(&A->imag()));
		}

		fields.AddScalar("B_Magnitude",
			std::make_unique<ComplexVectorMagnitudeCoefficient>(*b_re, *b_im));
		return fields;
	}

	void SaveAnalysis() override
	{
		if (config.AnalysisType == AnalysisType::CouplingMatrix) {
			WriteCouplingMatrix();
		}
	}

	void WriteCouplingMatrix() {
		// Placeholder: once the admittance matrix is assembled, write it via the
		// shared helper, e.g.:
		//   SaveCouplingMatrix(Y, "Admittance Matrix [S]", "admittance_matrix.csv");
		mfem::out << "WriteCouplingMatrix() not implemented yet.\n";
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
};