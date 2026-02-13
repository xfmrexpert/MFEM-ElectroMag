// Copyright (c) 2026 T. C. Raymond
// SPDX-License-Identifier: MIT

#pragma once

#include <memory>
#include "mfem.hpp"
#include "physics_solver.hpp"
#include "axisymmetric_curl_curl_integrator.hpp"
#include "axisymmetric_lf_integrator.hpp"
#include "magnetic_field_coefficient.hpp"
#include "input_parser.hpp"
#include "boundary_validation.hpp"

class MagnetostaticSolver : public PhysicsSolver {
    // 1. Logic
    enum class SolverType { Axisymmetric, Planar };
    SolverType type = SolverType::Axisymmetric; 

    // 2. Resources (Order of declaration = Order of destruction)
    std::unique_ptr<mfem::H1_FECollection> fec;
    std::unique_ptr<mfem::FiniteElementSpace> fespace;
    
    std::unique_ptr<mfem::GridFunction> A; // Vector Potential (scalar in 2D)
    
    // Materials & Sources
    std::unique_ptr<mfem::PWConstCoefficient> nu_coeff;
    std::unique_ptr<mfem::PWConstCoefficient> j_coeff;
    
    // System
    std::unique_ptr<mfem::LinearForm> b;
    mfem::Array<int> ess_bdr;

public:
    MagnetostaticSolver(mfem::Mesh &m, json &c) : PhysicsSolver(m, c) {}
    
    // No manual destructor needed!

    void Setup() override {
        // 1. Config & Solver Type
        int order = config["simulation"].value("order", 1);
        int dim = mesh.Dimension();

        std::string mode = config["simulation"].value("model_type", "axisymmetric");
        type = (mode == "planar") ? SolverType::Planar : SolverType::Axisymmetric;

        // 2. FEM Spaces
        fec = std::make_unique<mfem::H1_FECollection>(order, dim);
        fespace = std::make_unique<mfem::FiniteElementSpace>(&mesh, fec.get());
        
        A = std::make_unique<mfem::GridFunction>(fespace.get());
        *A = 0.0;

        // 3. Materials (Reluctivity) & Sources (Current)
        InputParser parser(config);
        
        mfem::Vector nu_vec;
        parser.SetupReluctivity(mesh, nu_vec);
        nu_coeff = std::make_unique<mfem::PWConstCoefficient>(nu_vec);

        mfem::Vector j_vec;
        parser.SetupSources(mesh, j_vec);
        j_coeff = std::make_unique<mfem::PWConstCoefficient>(j_vec);
        
        // 4. Boundary Conditions
        ess_bdr.SetSize(mesh.bdr_attributes.Max());
        ess_bdr = 0;

        std::vector<std::pair<mfem::Array<int>, double>> bcs;
        parser.SetupBoundaries(mesh, bcs);

        // Validate that BCs don't create physical conflicts
        BoundaryConditionValidator validator(mesh, *fespace);
        validator.ValidateBoundaryConditions(bcs, false);  // Strict mode - reject conflicts

        // Apply boundary conditions
        for (const auto& [marker, val] : bcs) {
            mfem::ConstantCoefficient val_coeff(val);
            A->ProjectBdrCoefficient(val_coeff, marker);

            for(int i=0; i<marker.Size(); i++) {
                if(marker[i]) ess_bdr[i] = 1;
            }
        }
        
        // 5. Linear Form (RHS)
        b = std::make_unique<mfem::LinearForm>(fespace.get());
        
        if (type == SolverType::Axisymmetric) {
            // Your custom axisymmetric integrator (Integrates J * v * 2*pi*r)
            b->AddDomainIntegrator(new AxisymmetricLFIntegrator(*j_coeff));
        } else {
            // Standard Planar Integrator (Integrates J * v * 1.0)
            b->AddDomainIntegrator(new mfem::DomainLFIntegrator(*j_coeff));
        }
        b->Assemble();
    }

    void Run() override {
        // 6. Bilinear Form (Stiffness)
        // Stack allocation is fine for BilinearForm as it lives only during solve
        mfem::BilinearForm a(fespace.get());
        
        if (type == SolverType::Axisymmetric) {
            // Custom: Integrates (1/mu) * Curl(A) * Curl(v) * 2*pi*r
            a.AddDomainIntegrator(new AxisymmetricCurlCurlIntegrator(*nu_coeff));
        } else {
            // Planar: Div( 1/mu * Grad(A) ) 
            a.AddDomainIntegrator(new mfem::DiffusionIntegrator(*nu_coeff));
        }
        a.Assemble();

        // 8. Solve
        mfem::OperatorPtr Matrix;
        mfem::Vector B, X;
        
        mfem::Array<int> ess_tdof_list;
        fespace->GetEssentialTrueDofs(ess_bdr, ess_tdof_list);

        a.FormLinearSystem(ess_tdof_list, *A, *b, Matrix, X, B);

        // Debug check
        if (B.Norml2() < 1e-12 && X.Norml2() < 1e-12) {
            mfem::out << "WARNING: Linear system RHS is zero. "
                      << "Check that 'sources' in JSON match Mesh Attributes." << std::endl;
        }

#ifdef MFEM_USE_SUITESPARSE
        mfem::UMFPackSolver umf_solver;
        umf_solver.Control[UMFPACK_PRL] = 1; 
        umf_solver.SetOperator(*Matrix);
        umf_solver.Mult(B, X);
#else
        InputParser parser(config);
        mfem::GSSmoother M((mfem::SparseMatrix&)(*Matrix));
        mfem::PCG(*Matrix, M, B, X, parser.GetSolverPrintLevel(),
                  parser.GetSolverMaxIter(), parser.GetSolverTolerance(), 0.0);
#endif

        a.RecoverFEMSolution(X, *b, *A);
    }
    
    void Save() override {
        mfem::ParaViewDataCollection paraview("results_magnetostatic", &mesh);
        paraview.SetLevelsOfDetail(1);
        paraview.RegisterField("A", A.get());
        
        // 1. Prepare storage for B-field (common to both types)
        mfem::L2_FECollection fec_l2(fec->GetOrder() - 1, mesh.Dimension());
        mfem::FiniteElementSpace fespace_l2(&mesh, &fec_l2, mesh.Dimension()); 
        mfem::GridFunction B_gf(&fespace_l2);

        // 2. Select the correct physics coefficient
        if (type == SolverType::Axisymmetric) {
            // Custom: B_r = -dA/dz, B_z = A/r + dA/dr
            MagneticFieldCoefficient B_coeff(*A);
            B_gf.ProjectCoefficient(B_coeff);
        } 
        else {
            // Planar: B_x = dA/dy, B_y = -dA/dx
            mfem::CurlGridFunctionCoefficient B_coeff(A.get());
            B_gf.ProjectCoefficient(B_coeff);
        }

        // 3. Register and Save
        paraview.RegisterField("B", &B_gf);
        
        paraview.SetCycle(0);
        paraview.SetTime(0.0);
        paraview.Save();
    }
};