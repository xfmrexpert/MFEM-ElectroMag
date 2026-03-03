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
    ModelType type = ModelType::Axisymmetric; // Default initialization

    double frequency = 60.0;
    mfem::real_t omega = Constants::TWO_PI * frequency;

    std::vector<Port> ports;

    // Smart Pointers for MFEM Objects
    // Order matters for destruction: GridFunctions depend on Spaces, Spaces depend on Collections.
    std::unique_ptr<mfem::H1_FECollection> fec;
    std::unique_ptr<mfem::FiniteElementSpace> fespace;
    
    // Complex System objects
    std::unique_ptr<mfem::SesquilinearForm> S_AA; 
    std::unique_ptr<mfem::ComplexGridFunction> A; 
    std::unique_ptr<mfem::ComplexLinearForm> b;

    // Coefficients
    std::unique_ptr<mfem::PWConstCoefficient> nu_coeff;
    std::unique_ptr<mfem::PWConstCoefficient> omega_sigma_coeff;
    std::unique_ptr<mfem::PWConstCoefficient> j_coeff;     
    
    mfem::Array<int> ess_bdr;

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
        if (type == ModelType::Axisymmetric)
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

public:
    // Constructor deals only with initialization, no manual nullptr assignment needed
    MagnetoquasistaticSolver(mfem::Mesh &m, const json &c) : PhysicsSolver(m, c) {}

    void Setup() override {
        // Config
        InputParser parser(config_json);
        config = parser.GetProblemConfig();
        
        int order = config.Order;
        frequency = config.Frequency;
        omega = Constants::TWO_PI * frequency;

        type = config.ModelType;

        // Spaces (make_unique)
        fec = std::make_unique<mfem::H1_FECollection>(order, mesh.Dimension());
        fespace = std::make_unique<mfem::FiniteElementSpace>(&mesh, fec.get());
        
        // Materials
        std::vector<Material> materials = config.Materials;
        // Real Part: Reluctivity
        mfem::Vector nu_vec(mesh.attributes.Max());
        mfem::Vector sigma_vec(mesh.attributes.Max());
        nu_vec = 0.0; sigma_vec = 0.0;
        
        //parser.SetupReluctivity(mesh, nu_vec);
        for (auto& region : config.Regions) {
            for (auto attribute_id : region.AttributeIds) {
                if (attribute_id > 0 && attribute_id <= mesh.attributes.Max()) {
                    auto& material = materials[region.Material];
                    double nu = 1.0 / (Constants::MU_0 * material.RelPermeability);
                    nu_vec[attribute_id - 1] = nu;
                    sigma_vec[attribute_id - 1] = material.Conductivity;
                }
            }
        }
        nu_coeff = std::make_unique<mfem::PWConstCoefficient>(nu_vec);

        // Imag Part: Omega * Sigma
        mfem::Vector omega_sigma_vec = sigma_vec;
        omega_sigma_vec *= omega;
        omega_sigma_coeff = std::make_unique<mfem::PWConstCoefficient>(omega_sigma_vec);

        // Ports
        // ports member variable is redundant if config maps ports, but kept for compatibility
        ports = config.Ports;

        // Source
        mfem::Vector j_src(mesh.attributes.Max());
        j_src = 0.0;
        
        for (const auto& src : config.Sources) {
            for (int attr : src.Markers) {
                if (attr > 0 && attr <= mesh.attributes.Max()) {
                    j_src[attr - 1] = src.CurrentDensity;
                }
            }
        }
        j_coeff = std::make_unique<mfem::PWConstCoefficient>(j_src);

        // Boundary Attributes
        ess_bdr.SetSize(mesh.bdr_attributes.Max());
        ess_bdr = 0;

        std::vector<std::pair<mfem::Array<int>, double>> bcs;
        for (const auto& bc : config.BoundaryConditions) {
             mfem::Array<int> marker(mesh.bdr_attributes.Max());
             marker = 0;
             for (int attr : bc.marker) {
                 if (attr > 0 && attr <= mesh.bdr_attributes.Max()) {
                     marker[attr - 1] = 1;
                 }
             }
             bcs.push_back({marker, bc.value});
        }
        
        // Validate that BCs don't create physical conflicts
        BoundaryConditionValidator validator(mesh, *fespace);
        validator.ValidateBoundaryConditions(bcs, false);  // Strict mode - reject conflicts

        // Mark essential boundaries (only those with zero BC in this formulation)
        for (const auto& [marker, val] : bcs) {
            if (val == 0.0) {
                 for(int i=0; i<marker.Size(); i++) if(marker[i]) ess_bdr[i] = 1;
            }
        }
    }

    void Run() override {
        // Setup Complex Billinear Form
        S_AA = std::make_unique<mfem::SesquilinearForm>(fespace.get(), mfem::ComplexOperator::HERMITIAN);

        if (type == ModelType::Axisymmetric) {
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

        mfem::BilinearForm &K = S_AA->real();
        mfem::BilinearForm &M_sigma = S_AA->imag();

        int N_DOFs = fespace->GetTrueVSize();
        int N_Ports = ports.size();

        mfem::Array<int> block_offsets(3);
        block_offsets[0] = 0;
        block_offsets[1] = N_DOFs;
        block_offsets[2] = N_DOFs + N_Ports;

        mfem::Array<mfem::Vector*> port_forms(N_Ports);

        mfem::DenseMatrix G_scaled(N_Ports, N_Ports);
        G_scaled = 0.0;

        // Generate the forms
        for (int i = 0; i < N_Ports; ++i)
        {
            Port port = config.Ports[i];
            Region region = config.Regions[port.Region];
            std::vector<int> attribute_ids = region.AttributeIds;
            Material material = config.Materials[region.Material];
            double conductivity = material.Conductivity;
            // BuildPortVector is the function defined in the previous step
            port_forms[i] = BuildPortVector(fespace.get(), attribute_ids, conductivity);
            G_scaled(i, i) = -1.0 / (omega * ComputePortConductance(attribute_ids, conductivity));
        }

        PortCouplingOperator C_op(N_DOFs, N_Ports, port_forms);
        
        mfem::ScaledOperator neg_C_op(&C_op, -1.0);
        mfem::TransposeOperator neg_C_T_op(&neg_C_op);

        mfem::BlockOperator M_Re(block_offsets);
        M_Re.SetBlock(0, 0, &K.SpMat());
        M_Re.SetBlock(0, 1, &neg_C_op);
        M_Re.SetBlock(1, 0, &neg_C_T_op);
        // M_Re Block (1, 1) is 0
        
        mfem::BlockOperator M_Im(block_offsets);
        M_Im.SetBlock(0, 0, &M_sigma.SpMat());
        M_Im.SetBlock(1, 1, &G_scaled);

        mfem::ComplexOperator global_complex_system(&M_Re, &M_Im, false, false);

        // Setup the complex RHS blocks
        mfem::BlockVector B_Re(block_offsets); B_Re = 0.0;
        mfem::BlockVector B_Im(block_offsets); B_Im = 0.0;

        // Assemble the source term (J is assumed real)
        mfem::LinearForm *b_source = new mfem::LinearForm(fespace.get());
        if (type == ModelType::Axisymmetric) {
             b_source->AddDomainIntegrator(new AxisymmetricLFIntegrator(*j_coeff));
        } else {
             b_source->AddDomainIntegrator(new mfem::DomainLFIntegrator(*j_coeff));
        }
        b_source->Assemble();

        // Add source to Real part of Mesh RHS (Block 0)
        B_Re.GetBlock(0) += *b_source;
        delete b_source;

        // Set the imaginary RHS for the active port 'k'
        if (N_Ports > 0)
        {
            B_Im.GetBlock(1)(0) = -1.0 / omega;
        }

        // Construct the full RHS vector and Initial Guess
        // System size is 2 * (N_DOFs + N_Ports)
        int total_size = 2 * (N_DOFs + N_Ports);
        mfem::Vector b_combined(total_size);
        mfem::Vector x_combined(total_size);
        
        x_combined = 0.0; // Initial guess
        
        // Copy parts to combined RHS
        // convention: [Re_Mesh, Re_Port, Im_Mesh, Im_Port]
        // This assumes ComplexOperator expects [Re_Block, Im_Block] where Block is the full coupled vector
        for(int i=0; i<B_Re.Size(); i++) b_combined(i) = B_Re(i);
        for(int i=0; i<B_Im.Size(); i++) b_combined(B_Re.Size() + i) = B_Im(i);

        // Grid Function (for solution recovery later)
        A = std::make_unique<mfem::ComplexGridFunction>(fespace.get());
        *A = 0.0;
        
        // Boundaries
        // Convert BC markers to DOF list
        mfem::Array<int> ess_tdof_list;
        fespace->GetEssentialTrueDofs(ess_bdr, ess_tdof_list);

        // Solve
        mfem::OperatorHandle A_op;
        mfem::Vector B_vec, X_vec;
        
        mfem::Operator *A_op_ptr;
        
        global_complex_system.FormLinearSystem(ess_tdof_list, x_combined, b_combined, A_op_ptr, X_vec, B_vec);
        A_op.Reset(A_op_ptr, false);

#ifdef MFEM_USE_SUITESPARSE
        // Direct Complex Solver
        mfem::ComplexUMFPackSolver solver;
        solver.Control[UMFPACK_PRL] = 1;
        solver.SetOperator(*A_op.Ptr());
        solver.Mult(B_vec, X_vec);
#else
        // Iterative Complex Solver
        mfem::GMRESSolver gmres;
        gmres.SetOperator(*A_op.Ptr());
        gmres.SetPrintLevel(config.SolverPrintLevel);
        gmres.SetRelTol(config.SolverTolerance);
        gmres.SetMaxIter(config.SolverMaxIter);
        gmres.Mult(B_vec, X_vec);
#endif

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

    void Save() override {
        mfem::ParaViewDataCollection paraview("results_mqs", &mesh);
        
        paraview.RegisterField("A_Real", &A->real());
        paraview.RegisterField("A_Imag", &A->imag());

        // Derived B-Fields
        mfem::FiniteElementSpace fespace_vec(&mesh, fec.get(), mesh.Dimension());
        mfem::GridFunction B_real(&fespace_vec);
        mfem::GridFunction B_imag(&fespace_vec); 
        
        if (type == ModelType::Axisymmetric) {
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