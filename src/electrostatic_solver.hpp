// Copyright (c) 2026 T. C. Raymond
// SPDX-License-Identifier: MIT

#pragma once
#include <memory>
#include "mfem.hpp"
#include "physics_solver.hpp"
#include "axisymmetric_diffusion_integrator.hpp"
#include "input_parser.hpp"
#include "boundary_validation.hpp"

class ElectrostaticSolver : public PhysicsSolver {
    ModelType type = ModelType::Axisymmetric; 

    // Primary Spaces
    std::unique_ptr<mfem::H1_FECollection> fec;
    std::unique_ptr<mfem::FiniteElementSpace> fespace;
    std::unique_ptr<mfem::GridFunction> x; // Electric Potential (V)
    
    // Physics
    std::unique_ptr<mfem::PWConstCoefficient> epsilon_coeff;
    std::unique_ptr<mfem::LinearForm> b;
    
    mfem::Array<int> ess_bdr;

public:
    ElectrostaticSolver(mfem::Mesh &m, const json &c) : PhysicsSolver(m, c) {}
    
    void Setup() override {
        // Config & solver type
        InputParser parser(config_json);
        config = parser.GetProblemConfig();
            
        int order = config.Order;
        const int dim   = mesh.Dimension();

        type = config.ModelType;

        // Spaces
        fec = std::make_unique<mfem::H1_FECollection>(order, dim);
        fespace = std::make_unique<mfem::FiniteElementSpace>(&mesh, fec.get());
        x = std::make_unique<mfem::GridFunction>(fespace.get());
        *x = 0.0;

        // Material Properties (Permittivity)
        std::vector<Material> materials = config.Materials;
        mfem::Vector epsilon_values(mesh.attributes.Max());
        epsilon_values = 0.0;
      
        for (auto& region : config.Regions) {
            for (auto attribute_id : region.AttributeIds) {
                if (attribute_id > 0 && attribute_id <= mesh.attributes.Max()) {
                    auto& material = materials[region.Material];
                    epsilon_values[attribute_id - 1] = material.RelPermittivity * Constants::EPSILON_0;
                }
            }
        }
        
        epsilon_coeff = std::make_unique<mfem::PWConstCoefficient>(epsilon_values);

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

        BoundaryConditionValidator validator(mesh, *fespace);
        validator.ValidateBoundaryConditions(bcs, /*allow_overlap=*/false);

        // Apply boundary conditions
        for (const auto& [marker, val] : bcs) {
            mfem::ConstantCoefficient val_coeff(val);
            x->ProjectBdrCoefficient(val_coeff, marker);

            for(int i=0; i<marker.Size(); i++) {
                if(marker[i]) ess_bdr[i] = 1;
            }
        }

        // 5. Linear Form (RHS)
        b = std::make_unique<mfem::LinearForm>(fespace.get());
        *b = 0.0;
    }

    void Run() override {
        // 6. Assemble Stiffness Matrix
        mfem::BilinearForm a(fespace.get());

        if (type == ModelType::Axisymmetric) {
            // Solves: Div( r * eps * Grad(V) ) = 0
            // Integrator handles 'r' weight and 'epsilon'
            a.AddDomainIntegrator(new AxisymmetricDiffusionIntegrator(*epsilon_coeff));
        }
        else {
            // Solves: Div( eps * Grad(V) ) = 0
            // Standard Cartesian Laplacian
            a.AddDomainIntegrator(new mfem::DiffusionIntegrator(*epsilon_coeff));
        }

        a.Assemble();

        // 7. Solve
        mfem::OperatorPtr A;
        mfem::Vector B, X;
        mfem::Array<int> ess_tdof_list;
        fespace->GetEssentialTrueDofs(ess_bdr, ess_tdof_list);

        a.FormLinearSystem(ess_tdof_list, *x, *b, A, X, B);

#ifdef MFEM_USE_SUITESPARSE
        mfem::UMFPackSolver umf_solver;
        umf_solver.Control[UMFPACK_PRL] = 1;
        umf_solver.SetOperator(*A);
        umf_solver.Mult(B, X);
#else
        mfem::GSSmoother M((mfem::SparseMatrix&)(*A));
        mfem::PCG(*A, M, B, X, config.SolverPrintLevel,
                  config.SolverMaxIter, config.SolverTolerance, 0.0);
#endif

        // Recover solution into GridFunction x
        a.RecoverFEMSolution(X, *b, *x);
    }
    
    void Save() override {
        mfem::ParaViewDataCollection paraview("results_electrostatic", &mesh);
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
};