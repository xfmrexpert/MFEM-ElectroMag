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

class MagnetoquasistaticSolver : public PhysicsSolver {
    // 1. Logic Types
    enum class SolverType { Axisymmetric, Planar };
    SolverType type = SolverType::Axisymmetric; // Default initialization

    double frequency = 60.0;
    double omega = Constants::TWO_PI * frequency;

    // 2. Smart Pointers for MFEM Objects
    // Order matters for destruction: GridFunctions depend on Spaces, Spaces depend on Collections.
    std::unique_ptr<mfem::H1_FECollection> fec;
    std::unique_ptr<mfem::FiniteElementSpace> fespace;
    
    // Complex System objects
    std::unique_ptr<mfem::SesquilinearForm> a; 
    std::unique_ptr<mfem::ComplexGridFunction> A; 
    std::unique_ptr<mfem::ComplexLinearForm> b;

    // Coefficients
    std::unique_ptr<mfem::PWConstCoefficient> nu_coeff;
    std::unique_ptr<mfem::PWConstCoefficient> omega_sigma_coeff;
    std::unique_ptr<mfem::PWConstCoefficient> j_coeff;     
    
    mfem::Array<int> ess_bdr;

public:
    // Constructor deals only with initialization, no manual nullptr assignment needed
    MagnetoquasistaticSolver(mfem::Mesh &m, const json &c) : PhysicsSolver(m, c) {}

    // Destructor is implicitly defined and handles cleanup automatically!
    // ~MagnetoquasistaticSolver() = default;

    void Setup() override {
        // 1. Config
        int order = config["simulation"].value("order", 1);
        frequency = config["simulation"].value("frequency", 60.0);
        omega = Constants::TWO_PI * frequency;

        std::string mode = config["simulation"].value("model_type", "axisymmetric");
        type = (mode == "planar") ? SolverType::Planar : SolverType::Axisymmetric;

        // 2. Spaces (make_unique)
        fec = std::make_unique<mfem::H1_FECollection>(order, mesh.Dimension());
        fespace = std::make_unique<mfem::FiniteElementSpace>(&mesh, fec.get());
        
        // 3. Materials
        InputParser parser(config);
        
        // Real Part: Reluctivity
        mfem::Vector nu_vec;
        parser.SetupReluctivity(mesh, nu_vec);
        nu_coeff = std::make_unique<mfem::PWConstCoefficient>(nu_vec);

        // Imag Part: Omega * Sigma
        mfem::Vector sigma_vec;
        parser.SetupConductivity(mesh, sigma_vec); 
        sigma_vec *= omega; 
        omega_sigma_coeff = std::make_unique<mfem::PWConstCoefficient>(sigma_vec);

        // Source
        mfem::Vector j_src;
        parser.SetupSources(mesh, j_src);
        j_coeff = std::make_unique<mfem::PWConstCoefficient>(j_src);

        // 4. Boundary Attributes
        ess_bdr.SetSize(mesh.bdr_attributes.Max());
        ess_bdr = 0;

        std::vector<std::pair<mfem::Array<int>, double>> bcs;
        parser.SetupBoundaries(mesh, bcs);

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
        // 1. Setup Complex Billinear Form
        // We use .get() to pass raw pointers to MFEM API where required
        a = std::make_unique<mfem::SesquilinearForm>(fespace.get(), mfem::ComplexOperator::HERMITIAN);

        if (type == SolverType::Axisymmetric) {
            // Real Part: Curl-Curl (1/mu) with r weight
            a->AddDomainIntegrator(new AxisymmetricCurlCurlIntegrator(*nu_coeff), nullptr);
            
            // Imag Part: Mass (sigma * omega) with r weight
            a->AddDomainIntegrator(nullptr, new AxisymmetricMassIntegrator(*omega_sigma_coeff));
        }
        else {
            // Planar 2D
            // Real Part: Diffusion (similar to Magnetostatic Planar)
            // Solves Div(nu Grad A)
            a->AddDomainIntegrator(new mfem::DiffusionIntegrator(*nu_coeff), nullptr);

            // Imag Part: Standard Mass (sigma * omega)
            a->AddDomainIntegrator(nullptr, new mfem::MassIntegrator(*omega_sigma_coeff));
        }
        
        a->Assemble();
        a->Finalize();

        // 2. Linear Form
        b = std::make_unique<mfem::ComplexLinearForm>(fespace.get());
        
        // Real source
        if (type == SolverType::Axisymmetric) {
            b->AddDomainIntegrator(new AxisymmetricLFIntegrator(*j_coeff), nullptr);
        } else {
            b->AddDomainIntegrator(new mfem::DomainLFIntegrator(*j_coeff), nullptr);
        }
        b->Assemble();

        // 3. Grid Function
        A = std::make_unique<mfem::ComplexGridFunction>(fespace.get());
        *A = 0.0;
        
        // 4. Boundaries
        mfem::ConstantCoefficient zero(0.0);
        A->ProjectBdrCoefficient(zero, zero, ess_bdr);

        // 5. Solve
        mfem::OperatorHandle A_op;
        mfem::Vector B_vec, X_vec;
        
        a->FormLinearSystem(ess_bdr, *A, *b, A_op, X_vec, B_vec);

#ifdef MFEM_USE_SUITESPARSE
        // Direct Complex Solver
        mfem::ComplexUMFPackSolver solver;
        solver.Control[UMFPACK_PRL] = 1;
        solver.SetOperator(*A_op.Ptr());
        solver.Mult(B_vec, X_vec);
#else
        // Iterative Complex Solver
        InputParser parser(config);
        mfem::GMRESSolver gmres;
        gmres.SetOperator(*A_op.Ptr());
        gmres.SetPrintLevel(parser.GetSolverPrintLevel());
        gmres.SetRelTol(parser.GetSolverTolerance());
        gmres.SetMaxIter(parser.GetSolverMaxIter());
        gmres.Mult(B_vec, X_vec);
#endif

        // 6. Extract Solution using high-level method
        a->RecoverFEMSolution(X_vec, *b, *A);
    }

    void Save() override {
        mfem::ParaViewDataCollection paraview("results_mqs", &mesh);
        
        paraview.RegisterField("A_Real", &A->real());
        paraview.RegisterField("A_Imag", &A->imag());

        // Derived B-Fields
        mfem::FiniteElementSpace fespace_vec(&mesh, fec.get(), mesh.Dimension());
        mfem::GridFunction B_real(&fespace_vec);
        mfem::GridFunction B_imag(&fespace_vec); 
        
        if (type == SolverType::Axisymmetric) {
            // Axisymmetric B = Curl(A_phi) = (-dA/dz, 1/r*d(rA)/dr)
            MagneticFieldCoefficient B_real_coeff(A->real());
            MagneticFieldCoefficient B_imag_coeff(A->imag());
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