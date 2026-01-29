#pragma once
#include <memory>
#include "mfem.hpp"
#include "physics_solver.hpp"
#include "axisymmetric_diffusion_integrator.hpp"
#include "input_parser.hpp"

class ElectrostaticSolver : public PhysicsSolver {
    enum class SolverType { Axisymmetric, Planar };
    SolverType type = SolverType::Axisymmetric; 

    // Primary Spaces
    std::unique_ptr<mfem::H1_FECollection> fec;
    std::unique_ptr<mfem::FiniteElementSpace> fespace;
    std::unique_ptr<mfem::GridFunction> x; // Electric Potential (V)
    
    // Physics
    std::unique_ptr<mfem::PWConstCoefficient> epsilon_coeff;
    std::unique_ptr<mfem::LinearForm> b;
    
    mfem::Array<int> ess_bdr;

public:
    ElectrostaticSolver(mfem::Mesh &m, json &c) : PhysicsSolver(m, c) {}
    
    void Setup() override {
        // 1. Config & Logic
        int order = config["simulation"].value("order", 1);
        int dim = mesh.Dimension(); // Should be 2

        std::string mode = config["simulation"].value("model_type", "axisymmetric");
        type = (mode == "planar") ? SolverType::Planar : SolverType::Axisymmetric;

        // 2. Spaces
        fec = std::make_unique<mfem::H1_FECollection>(order, dim);
        fespace = std::make_unique<mfem::FiniteElementSpace>(&mesh, fec.get());
        x = std::make_unique<mfem::GridFunction>(fespace.get());
        *x = 0.0;
        
        // 3. Material Properties (Permittivity)
        InputParser parser(config);
        mfem::Vector epsilon_values;
        parser.SetupPermittivity(mesh, epsilon_values);
        epsilon_coeff = std::make_unique<mfem::PWConstCoefficient>(epsilon_values);

        // 4. Boundary Conditions
        ess_bdr.SetSize(mesh.bdr_attributes.Max());
        ess_bdr = 0;
        
        std::vector<std::pair<mfem::Array<int>, double>> bcs;
        parser.SetupBoundaries(mesh, bcs);

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

        if (type == SolverType::Axisymmetric) {
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
        mfem::PCG(*A, M, B, X, 0, 2000, 1e-12, 0.0);
#endif

        // Recover solution into GridFunction x
        a.RecoverFEMSolution(X, *b, *x);
    }
    
    void Save() override {
        mfem::ParaViewDataCollection paraview("results_electrostatic", &mesh);
        paraview.SetLevelsOfDetail(1);
        paraview.RegisterField("V", x.get());
        
        // --- Electric Field & Visualization Spaces ---
        // Stack allocate temporaries that are managed automatically
        mfem::L2_FECollection fec_l2(fec->GetOrder() - 1, mesh.Dimension());
        
        // 1. Electric Field Vector Space
        mfem::FiniteElementSpace fespace_l2_vec(&mesh, &fec_l2, mesh.Dimension()); 
        mfem::GridFunction E(&fespace_l2_vec);
        
        mfem::GradientGridFunctionCoefficient grad_x(x.get());
        E.ProjectCoefficient(grad_x);
        
        // Negate to get E = -Grad(V)
        for (int i = 0; i < E.Size(); i++) {
            E[i] = -E[i];
        }
        paraview.RegisterField("E", &E);
        
        // 2. Permittivity Scalar Space
        mfem::FiniteElementSpace fespace_l2_scalar(&mesh, &fec_l2);
        mfem::GridFunction eps_gf(&fespace_l2_scalar);
        eps_gf.ProjectCoefficient(*epsilon_coeff);
        paraview.RegisterField("Permittivity", &eps_gf);

        paraview.SetCycle(0);
        paraview.SetTime(0.0);
        paraview.Save();
    }
};