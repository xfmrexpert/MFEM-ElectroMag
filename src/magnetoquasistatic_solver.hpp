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
#include "input_parser.hpp"
#include "constants.hpp"
#include "boundary_validation.hpp"
#include "problem_config.hpp"
#include "port_coupling_operator.hpp"
#include "axisymmetric_conductance_coefficient.hpp"

class MagnetoquasistaticSolver : public PhysicsSolver {
    GeometryType geometry = GeometryType::Axisymmetric; // Default initialization

    double frequency = 60.0;
    mfem::real_t omega = Constants::TWO_PI * frequency;

    std::vector<Terminal> terminals;

    // Smart Pointers for MFEM Objects
    // Order matters for destruction: GridFunctions depend on Spaces, Spaces depend on Collections.
    std::unique_ptr<mfem::H1_FECollection> fec;
    std::unique_ptr<mfem::FiniteElementSpace> fespace;
    
    // Complex System objects
    // Declared before the BlockOperators/BlockVectors below: they hold a
    // non-owning reference to this offsets array, so it must outlive them
    // (members destroy in reverse declaration order).
    std::unique_ptr<mfem::Array<int>> block_offsets;
    std::unique_ptr<mfem::SesquilinearForm> S_AA; 
    std::unique_ptr<mfem::ComplexGridFunction> A; 
    std::unique_ptr<mfem::ComplexLinearForm> b;
    std::unique_ptr<mfem::ComplexOperator> global_complex_system;
    std::unique_ptr<mfem::BlockVector> B_Re;
	std::unique_ptr<mfem::BlockVector> B_Im;
    std::unique_ptr<mfem::Vector> x_combined;
	std::unique_ptr<mfem::Vector> b_combined;
    std::unique_ptr<mfem::Array<mfem::Vector*>> port_forms;
    std::unique_ptr<mfem::BlockOperator> M_Re;
    std::unique_ptr<mfem::BlockOperator> M_Im;
    std::unique_ptr<PortCouplingOperator> C_op;
    std::unique_ptr<mfem::ScaledOperator> neg_C_op;
    std::unique_ptr<mfem::TransposeOperator> neg_C_T_op;
    std::unique_ptr<mfem::DenseMatrix> G_scaled;

    int N_DOFs;
    int N_Ports;

    // Coefficients
    std::unique_ptr<mfem::PWConstCoefficient> nu_coeff;
    std::unique_ptr<mfem::PWConstCoefficient> omega_sigma_coeff;
    std::unique_ptr<mfem::PWConstCoefficient> j_coeff;     
    
    mfem::Array<int> ess_bdr;
    mfem::Array<int> ess_tdof_list;

    // Function to build the port vector for a specific port attribute
    mfem::Vector* BuildPortVector(mfem::FiniteElementSpace* fespace, 
                            std::vector<int> port_attributes, 
                            double sigma) 
    {
        // Restrict integration to this specific port's attribute
        mfem::Array<int> port_marker(fespace->GetMesh()->attributes.Max());
        port_marker = 0;
        for (auto port_attribute : port_attributes) {
            port_marker[port_attribute - 1] = 1; // MFEM attributes are 1-indexed
        }

        // Define the appropriate coefficient
        mfem::Coefficient* coeff = nullptr;
        // For Axisymmetric, see constant coefficient reasoning
        coeff = new mfem::ConstantCoefficient(sigma);
        
        // Apply the restriction marker so it only evaluates on the port
        mfem::RestrictedCoefficient restricted_coeff(*coeff, port_marker);

        // Assemble the LinearForm using a scalar domain integrator
        mfem::LinearForm port_lf(fespace);
        port_lf.AddDomainIntegrator(new mfem::DomainLFIntegrator(restricted_coeff));
        port_lf.Assemble();

        // Extract and return as a standalone Vector
        mfem::Vector* port_vector = new mfem::Vector(port_lf.Size());
        *port_vector = port_lf;

        delete coeff;
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

    // Geometric fallback: find boundary attributes whose boundary elements lie on r=0 and mark them essential.
   // This is intentionally conservative. Best practice is to tag the axis in your mesh and handle it in InputParser.
   void MarkAxisBoundaryAttributesGeometric()
   {
      const double tol = Constants::AXIS_TOLERANCE;
      if (!mesh.bdr_attributes.Size()) { return; }

      mfem::Array<int> axis_attr(mesh.bdr_attributes.Max());
      axis_attr = 0;

      for (int be = 0; be < mesh.GetNBE(); be++)
      {
         mfem::Element *bEl = mesh.GetBdrElement(be);
         const int attr = bEl->GetAttribute();
         mfem::Array<int> v;
         bEl->GetVertices(v);

         bool on_axis = true;
         for (int i = 0; i < v.Size(); i++)
         {
            const double *vx = mesh.GetVertex(v[i]);
            if (std::abs(vx[0]) > tol) { on_axis = false; break; }
         }

         if (on_axis)
         {
            axis_attr[attr - 1] = 1; // attributes are 1-based
         }
      }

      // Merge axis boundary attributes into ess_bdr
      for (int i = 0; i < axis_attr.Size(); i++)
      {
         if (axis_attr[i]) { ess_bdr[i] = 1; }
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
        fespace = std::make_unique<mfem::FiniteElementSpace>(&mesh, fec.get());
        
        // Materials
        std::vector<Material> materials = config.Materials;
        // Real Part: Reluctivity
        mfem::Vector nu_vec(mesh.attributes.Max());
        mfem::Vector sigma_vec(mesh.attributes.Max());
        nu_vec = 0.0; sigma_vec = 0.0;
        
        for (auto& region : config.Regions) {
            const std::string& group_name = region.EntityGroupName;
            const EntityGroup& group = config.EntityGroups.at(group_name);

            for (auto attribute_id : group.AttributeIds) {
                if (attribute_id > 0 && attribute_id <= mesh.attributes.Max()) {
                    auto& material = materials[region.Material];
                    nu_vec[attribute_id - 1] = 1.0 / (Constants::MU_0 * material.RelPermeability);
                    sigma_vec[attribute_id - 1] = material.Conductivity;
                }
            }
        }

        nu_coeff = std::make_unique<mfem::PWConstCoefficient>(nu_vec);

        // Imag Part: Omega * Sigma
        mfem::Vector omega_sigma_vec = sigma_vec;
        omega_sigma_vec *= omega;
        omega_sigma_coeff = std::make_unique<mfem::PWConstCoefficient>(omega_sigma_vec);

        // Boundary Attributes
        std::vector<std::pair<mfem::Array<int>, double>> bcs;
        for (const auto& bc : config.BoundaryConditions) {
			const std::string& group_name = bc.EntityGroupName;
			const EntityGroup& group = config.EntityGroups.at(group_name);
            mfem::Array<int> marker(mesh.bdr_attributes.Max());
            marker = MarkerFromAttrs(group.AttributeIds);
            bcs.push_back({ marker, bc.Value });
        }

        // Validate that BCs don't create physical conflicts
        BoundaryConditionValidator validator(mesh, *fespace);
        validator.ValidateBoundaryConditions(bcs, /*terminals=*/{}, false);  // Strict mode - reject conflicts

        ess_bdr.SetSize(mesh.bdr_attributes.Max());
        ess_bdr = 0;
        // Mark essential boundaries (only those with zero BC in this formulation)
        for (const auto& [marker, val] : bcs) {
            if (val == 0.0) {
                for (int i = 0; i < marker.Size(); i++) if (marker[i]) ess_bdr[i] = 1;
            }
        }

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
        N_Ports = config.Terminals.size();

        block_offsets = std::make_unique<mfem::Array<int>>(3);
        (*block_offsets)[0] = 0;
        (*block_offsets)[1] = N_DOFs;
        (*block_offsets)[2] = N_DOFs + N_Ports;

        port_forms = std::make_unique<mfem::Array<mfem::Vector*>>(N_Ports);

        G_scaled = std::make_unique<mfem::DenseMatrix>(N_Ports, N_Ports);
        *G_scaled = 0.0;

		// Generate the forms
		int i = 0;
		for (const auto& [term_name, term] : config.Terminals)
		{
			const std::string& group_name = term.EntityGroupName;
			const EntityGroup& group = config.EntityGroups.at(group_name);
			std::vector<int> attribute_ids = group.AttributeIds;
			//Material material = config.Materials[region.Material];
			//double conductivity = material.Conductivity;
			double conductivity = 0.0;
			// BuildPortVector is the function defined in the previous step
			(*port_forms)[i] = BuildPortVector(fespace.get(), attribute_ids, conductivity);
			(*G_scaled)(i, i) = -1.0 / (omega * ComputePortConductance(attribute_ids, conductivity));
			++i;
		}

        C_op = std::make_unique<PortCouplingOperator>(N_DOFs, N_Ports, (*port_forms));

        neg_C_op = std::make_unique<mfem::ScaledOperator>(C_op.get(), -1.0);
        neg_C_T_op = std::make_unique<mfem::TransposeOperator>(neg_C_op.get());

        M_Re = std::make_unique<mfem::BlockOperator>(*block_offsets);
        M_Re->SetBlock(0, 0, &K.SpMat());
        M_Re->SetBlock(0, 1, neg_C_op.get());
        M_Re->SetBlock(1, 0, neg_C_T_op.get());
        // M_Re Block (1, 1) is 0

        M_Im = std::make_unique<mfem::BlockOperator>(*block_offsets);
        M_Im->SetBlock(0, 0, &M_sigma.SpMat());
        M_Im->SetBlock(1, 1, G_scaled.get());

        global_complex_system = std::make_unique<mfem::ComplexOperator>(M_Re.get(), M_Im.get(), false, false);

        // Setup the complex RHS blocks
        B_Re = std::make_unique<mfem::BlockVector>(*block_offsets);
        B_Im = std::make_unique<mfem::BlockVector>(*block_offsets);

        // Grid Function (for solution recovery later)
        A = std::make_unique<mfem::ComplexGridFunction>(fespace.get());
    }

    void ImprintScenario(const Scenario& sc) {
        *A = 0.0;
        *B_Re = 0.0;
        *B_Im = 0.0;

        // Source
        mfem::Vector j_src(mesh.attributes.Max());
        j_src = 0.0;

        for (const auto& term : config.Terminals) {
            //for (int attr : src.Markers) {
            //    if (attr > 0 && attr <= mesh.attributes.Max()) {
            //        j_src[attr - 1] = src.CurrentDensity;
            //    }
            //}
        }
        j_coeff = std::make_unique<mfem::PWConstCoefficient>(j_src);

        // Assemble the source term (J is assumed real)
        mfem::LinearForm* b_source = new mfem::LinearForm(fespace.get());
        if (geometry == GeometryType::Axisymmetric) {
            b_source->AddDomainIntegrator(new AxisymmetricLFIntegrator(*j_coeff));
        }
        else {
            b_source->AddDomainIntegrator(new mfem::DomainLFIntegrator(*j_coeff));
        }
        b_source->Assemble();

        // Add source to Real part of Mesh RHS (Block 0)
        B_Re->GetBlock(0) += *b_source;
        delete b_source;

        // Set the imaginary RHS for the active port 'k'
        if (N_Ports > 0)
        {
            B_Im->GetBlock(1)(0) = -1.0 / omega;
        }

        // Construct the full RHS vector and Initial Guess
        // System size is 2 * (N_DOFs + N_Ports)
        int total_size = 2 * (N_DOFs + N_Ports);
        b_combined = std::make_unique<mfem::Vector>(total_size);
        x_combined = std::make_unique<mfem::Vector>(total_size);

        *x_combined = 0.0; // Initial guess

        // Copy parts to combined RHS
        // convention: [Re_Mesh, Re_Port, Im_Mesh, Im_Port]
        // This assumes ComplexOperator expects [Re_Block, Im_Block] where Block is the full coupled vector
        for (int i = 0; i < B_Re->Size(); i++) (*b_combined)(i) = (*B_Re)(i);
        for (int i = 0; i < B_Im->Size(); i++) (*b_combined)(B_Re->Size() + i) = (*B_Im)(i);

    }

    void Run() override {
        if (config.AnalysisType == AnalysisType::CouplingMatrix) {
            // For coupling matrix, we solve one scenario per terminal with a unit drive
            for (const auto& [term_name, term] : config.Terminals) {
                //*x = 0.0; // Reset solution for new scenario
                //auto marker = MarkerFromAttrs(term.AttributeIds);
                //mfem::ConstantCoefficient c(1.0); // Unit drive
                //x->ProjectBdrCoefficient(c, marker);
                SolveSystem();
                SaveScenario("CouplingMatrix_" + term_name);
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

    void SolveSystem() {
        // Solve
        mfem::OperatorHandle A_op;
        mfem::Vector B_vec, X_vec;

        mfem::Operator* A_op_ptr;

        global_complex_system->FormLinearSystem(ess_tdof_list, *x_combined, *b_combined, A_op_ptr, X_vec, B_vec);
        bool own_A = (A_op_ptr != global_complex_system.get());
        A_op.Reset(A_op_ptr, own_A);

        // Iterative Complex Solver
        mfem::GMRESSolver gmres;
        gmres.SetOperator(*A_op.Ptr());
        gmres.SetPrintLevel(config.SolverPrintLevel);
        gmres.SetRelTol(config.SolverTolerance);
        gmres.SetMaxIter(config.SolverMaxIter);
        gmres.Mult(B_vec, X_vec);

        // Extract Solution manually
        // System solution X_vec is ordered: [Re_Mesh, Re_Port, Im_Mesh, Im_Port]
        // Copy mesh parts to GridFunction A

        // Real Part
        for (int i = 0; i < N_DOFs; i++) {
            A->real()(i) = X_vec(i);
        }

        // Imag Part (starts after Re_Mesh + Re_Port_Re)
        // Wait, Block offsets for Re part: [0, N_DOFs, N_DOFs+N_Ports]
        // So Im part starts at (N_DOFs + N_Ports)
		int offset_imag = N_DOFs + N_Ports;
		for (int i = 0; i < N_DOFs; i++) {
			A->imag()(i) = X_vec(offset_imag + i);
		}
	}

	void SaveScenario(const std::string& scenario_name) override {
		if (config.OutputParaview) {
			WriteParaviewResultsFile(scenario_name);
		}
	}

	void SaveAnalysis() override {}

	double TerminalConductivity(const Terminal& term) {
		// For simplicity, use the conductivity of the first material region associated with the terminal's attributes
		const std::string& group_name = term.EntityGroupName;
		const EntityGroup& group = config.EntityGroups.at(group_name);
		for (int attr : group.AttributeIds) {
			if (attr > 0 && attr <= mesh.attributes.Max()) {
				int region_id = -1;
				for (size_t i = 0; i < config.Regions.size(); i++) {
					auto eg_it = config.EntityGroups.find(config.Regions[i].EntityGroupName);
					if (eg_it == config.EntityGroups.end()) continue;
					const auto& region_attrs = eg_it->second.AttributeIds;
					if (std::find(region_attrs.begin(), region_attrs.end(), attr) != region_attrs.end()) {
						region_id = i;
						break;
					}
				}
				if (region_id != -1) {
					int material_id = config.Regions[region_id].Material;
					return config.Materials[material_id].Conductivity;
				}
			}
		}
		return 0.0; // Default to non-conductive if no match found
	}

    void WriteParaviewResultsFile(const std::string& scenario_name) {
        mfem::ParaViewDataCollection paraview("results_mqs_" + scenario_name, &mesh);
        
        paraview.RegisterField("A_Real", &A->real());
        paraview.RegisterField("A_Imag", &A->imag());

        // Derived B-Fields
        mfem::FiniteElementSpace fespace_vec(&mesh, fec.get(), mesh.Dimension());
        mfem::GridFunction B_real(&fespace_vec);
        mfem::GridFunction B_imag(&fespace_vec); 
        
        if (geometry == GeometryType::Axisymmetric) {
            // Axisymmetric B = Curl(A_phi) = (-dA/dz, 1/r*d(rA)/dr)
            MagneticFieldCoefficient B_real_coeff(&A->real());
            MagneticFieldCoefficient B_imag_coeff(&A->imag());
            B_real.ProjectCoefficient(B_real_coeff);
            B_imag.ProjectCoefficient(B_imag_coeff);
        }
        else {
            // Planar B = Curl(A_z) = (dA/dy, -dA/dx)
            mfem::CurlGridFunctionCoefficient B_real_coeff(&A->real());
            mfem::CurlGridFunctionCoefficient B_imag_coeff(&A->imag());
            B_real.ProjectCoefficient(B_real_coeff);
            B_imag.ProjectCoefficient(B_imag_coeff);
        }
        
        paraview.RegisterField("B_Real", &B_real);
        paraview.RegisterField("B_Imag", &B_imag);
        
        // Compute magnitude
        mfem::GridFunction B_mag(fespace.get()); 
        int ndofs = fespace->GetNDofs();
        int v_dim = fespace_vec.GetVDim();
        
        for (int i = 0; i < ndofs; i++) {
            double Br_re = B_real(i);
            double Bz_re = B_real(i + ndofs);
            
            double Br_im = B_imag(i);
            double Bz_im = B_imag(i + ndofs);
            
            double mag_sq = (Br_re * Br_re + Bz_re * Bz_re) + 
                            (Br_im * Br_im + Bz_im * Bz_im);

            B_mag(i) = std::sqrt(mag_sq);
        }
        paraview.RegisterField("B_Magnitude", &B_mag);
        
        paraview.Save();
    }
};