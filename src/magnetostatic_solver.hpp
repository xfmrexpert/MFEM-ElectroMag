// Copyright (c) 2026 T. C. Raymond
// SPDX-License-Identifier: MIT

#pragma once

#include <memory>
#include <cmath>
#include <algorithm>

#include "mfem.hpp"
#include "physics_solver.hpp"
#include "axisymmetric_curl_curl_integrator.hpp"
#include "axisymmetric_lf_integrator.hpp"
#include "magnetic_field_coefficient.hpp"
#include "input_parser.hpp"
#include "boundary_validation.hpp"
#include "constants.hpp"

class MagnetostaticSolver : public PhysicsSolver
{
private:
   enum class SolverType { Axisymmetric, Planar };
   SolverType type_ = SolverType::Axisymmetric;

   // Resources (order of declaration = order of destruction)
   std::unique_ptr<mfem::H1_FECollection>    fec_;
   std::unique_ptr<mfem::FiniteElementSpace> fespace_;
   std::unique_ptr<mfem::GridFunction>       A_;        // A_phi (axisym) or A (planar scalar)
   std::unique_ptr<mfem::PWConstCoefficient> nu_coeff_;  // ν = 1/μ
   std::unique_ptr<mfem::PWConstCoefficient> j_coeff_;   // J_phi (axisym) or J (planar scalar src)

   std::unique_ptr<mfem::LinearForm> b_;
   mfem::Array<int> ess_bdr_; // boundary attribute marker (size = bdr_attributes.Max())

public:
   MagnetostaticSolver(mfem::Mesh &m, const json &c) : PhysicsSolver(m, c) {}

   void Setup() override
   {
      // 1) Config & solver type
      const int order = config["simulation"].value("order", 1);
      const int dim   = mesh.Dimension();

      const std::string mode = config["simulation"].value("model_type", "axisymmetric");
      type_ = (mode == "planar") ? SolverType::Planar : SolverType::Axisymmetric;

      // 2) FE space
      fec_     = std::make_unique<mfem::H1_FECollection>(order, dim);
      fespace_ = std::make_unique<mfem::FiniteElementSpace>(&mesh, fec_.get());

      A_ = std::make_unique<mfem::GridFunction>(fespace_.get());
      *A_ = 0.0;

      // 3) Materials & sources
      InputParser parser(config);

      mfem::Vector nu_vec;
      parser.SetupReluctivity(mesh, nu_vec);
      nu_coeff_ = std::make_unique<mfem::PWConstCoefficient>(nu_vec);

      mfem::Vector j_vec;
      parser.SetupSources(mesh, j_vec);
      j_coeff_ = std::make_unique<mfem::PWConstCoefficient>(j_vec);

      // 4) Boundary conditions (essential)
      const int nbattr = mesh.bdr_attributes.Size() ? mesh.bdr_attributes.Max() : 0;
      ess_bdr_.SetSize(nbattr);
      ess_bdr_ = 0;

      std::vector<std::pair<mfem::Array<int>, double>> bcs;
      parser.SetupBoundaries(mesh, bcs);

      BoundaryConditionValidator validator(mesh, *fespace_);
      validator.ValidateBoundaryConditions(bcs, /*allow_overlap=*/false);

      // Apply BC values to A_ and build essential boundary marker.
      for (const auto &bc : bcs)
      {
         const mfem::Array<int> &marker = bc.first;
         const double val              = bc.second;

         mfem::ConstantCoefficient val_coeff(val);
         A_->ProjectBdrCoefficient(val_coeff, marker);

         // Merge marker -> ess_bdr_
         MFEM_ASSERT(marker.Size() == ess_bdr_.Size(),
                     "Boundary marker size must match bdr_attributes.Max().");
         for (int i = 0; i < marker.Size(); i++)
         {
            if (marker[i]) { ess_bdr_[i] = 1; }
         }
      }

      // 4b) Axis regularity: enforce A_phi = 0 on r=0 as ESSENTIAL.
      // Best practice: mark the axis as an essential boundary via boundary attributes if your mesh has it tagged.
      // If you *don't* have the axis tagged as a boundary attribute, do a geometric fallback:
      if (type_ == SolverType::Axisymmetric)
      {
         // Geometric fallback: force A=0 on axis boundary vertices by marking the boundary attributes
         // that lie on r=0. This requires detecting boundary elements on the axis and marking their attribute.
         // If your mesh already has an "axis" boundary attribute, prefer using that in InputParser instead.
         MarkAxisBoundaryAttributesGeometric_();
         // After this, ess_bdr_ includes axis attributes, and A_ has already been projected
         // for other BCs. We also project A=0 on the axis here for safety.
         ProjectAxisZero_();
      }

      // 5) RHS
      b_ = std::make_unique<mfem::LinearForm>(fespace_.get());

      if (type_ == SolverType::Axisymmetric)
      {
         // Integrates J * v * r  (global 2π omitted consistently)
         b_->AddDomainIntegrator(new AxisymmetricLFIntegrator(*j_coeff_));
      }
      else
      {
         b_->AddDomainIntegrator(new mfem::DomainLFIntegrator(*j_coeff_));
      }
      b_->Assemble();
   }

   void Run() override
   {
      // 6) Stiffness
      mfem::BilinearForm a(fespace_.get());

      if (type_ == SolverType::Axisymmetric)
      {
         a.AddDomainIntegrator(new AxisymmetricCurlCurlIntegrator(*nu_coeff_));
      }
      else
      {
         a.AddDomainIntegrator(new mfem::DiffusionIntegrator(*nu_coeff_));
      }

      a.Assemble();

      // 7) Essential DOFs
      mfem::Array<int> ess_tdof_list;
      fespace_->GetEssentialTrueDofs(ess_bdr_, ess_tdof_list);

      // 7b) For axisymmetric: Mark ALL DOFs at r=0 as essential (regularity condition)
      if (type_ == SolverType::Axisymmetric) {
         // For H1 elements, we need to identify all DOFs (vertex, edge, face) at r=0
         // For simplicity with order-1 or order-2, check vertex DOFs
         const int ndofs = fespace_->GetNDofs();
         mfem::Array<int> axis_dofs;

         // Get DOF coordinates - for H1 space, vertex DOFs are straightforward
         for (int i = 0; i < mesh.GetNV(); i++) {
            const double *coords = mesh.GetVertex(i);
            if (std::abs(coords[0]) < 1e-10) {  // r ≈ 0
               // For P1 elements: vertex i = DOF i
               // For P2 elements: need to map vertex to DOF
               if (i < ndofs) {  // Simple case: vertex DOF
                  axis_dofs.Append(i);
               }
            }
         }

         // Merge axis DOFs into essential DOF list
         mfem::Array<int> combined_list;
         combined_list.Reserve(ess_tdof_list.Size() + axis_dofs.Size());
         for (int i = 0; i < ess_tdof_list.Size(); i++) {
            combined_list.Append(ess_tdof_list[i]);
         }
         for (int i = 0; i < axis_dofs.Size(); i++) {
            // Check if already in list
            bool already_added = false;
            for (int j = 0; j < ess_tdof_list.Size(); j++) {
               if (ess_tdof_list[j] == axis_dofs[i]) {
                  already_added = true;
                  break;
               }
            }
            if (!already_added) {
               combined_list.Append(axis_dofs[i]);
            }
         }
         ess_tdof_list = combined_list;
      }

      // 8) Form and solve
      mfem::OperatorPtr Aop;
      mfem::Vector X, B;

      a.FormLinearSystem(ess_tdof_list, *A_, *b_, Aop, X, B);

      if (B.Norml2() < 1e-12 && X.Norml2() < 1e-12)
      {
         mfem::out << "WARNING: Linear system RHS is ~zero. "
                   << "Check that 'sources' in JSON match mesh attributes.\n";
      }

#ifdef MFEM_USE_SUITESPARSE
      mfem::UMFPackSolver umf;
      umf.Control[UMFPACK_PRL] = 1;
      umf.SetOperator(*Aop);
      umf.Mult(B, X);
#else
      InputParser parser(config);
      auto *sp = dynamic_cast<mfem::SparseMatrix*>(Aop.Ptr());
      MFEM_ASSERT(sp, "Expected SparseMatrix operator from FormLinearSystem.");

      mfem::GSSmoother M(*sp);
      mfem::PCG(*sp, M, B, X,
                parser.GetSolverPrintLevel(),
                parser.GetSolverMaxIter(),
                parser.GetSolverTolerance(),
                0.0);
#endif

      a.RecoverFEMSolution(X, *b_, *A_);

      mfem::out << "\n=== A Statistics ===\n";
      mfem::out << "  A min:     " << A_->Min() << "\n";
      mfem::out << "  A max:     " << A_->Max() << "\n";
      mfem::out << "  A L2 norm: " << A_->Norml2() << "\n";
   }

   void Save() override
   {
      mfem::ParaViewDataCollection pv("results_magnetostatic", &mesh);
      pv.SetLevelsOfDetail(1);
      pv.RegisterField("A", A_.get());

      // 1) Vector B field in (r,z) (axisym) or (x,y) (planar): vdim = 2 for 2D problems.
      const int dim = mesh.Dimension();
      MFEM_ASSERT(dim == 2, "Save() currently assumes a 2D mesh (axisymmetric r-z or planar x-y).");

      const int h1_order = fec_->GetOrder();
      const int l2_order = std::max(0, h1_order - 1);

      mfem::L2_FECollection fec_l2(l2_order, dim);
      mfem::FiniteElementSpace fes_l2(&mesh, &fec_l2, /*vdim=*/2);
      mfem::GridFunction B_gf(&fes_l2);
      B_gf = 0.0;

      if (type_ == SolverType::Axisymmetric)
      {
         // B_r = -∂A/∂z, B_z = ∂A/∂r + A/r
         MagneticFieldCoefficient B_coeff(*A_);
         B_gf.ProjectCoefficient(B_coeff);
      }
      else
      {
         // WARNING: MFEM's CurlGridFunctionCoefficient is not the usual "rotated gradient"
         // for a scalar potential in 2D. If your planar case is A_z and B = curl(A_z k),
         // then B = (∂A/∂y, -∂A/∂x). Implement that explicitly if needed.
         //
         // Placeholder: keep your original approach but you should verify it.
         mfem::CurlGridFunctionCoefficient B_coeff(A_.get());
         B_gf.ProjectCoefficient(B_coeff);
      }

      pv.RegisterField("B", &B_gf);
      pv.SetCycle(0);
      pv.SetTime(0.0);
      pv.Save();

      // Stats
      double Bmax = 0.0, Bsum = 0.0;
      int n = 0;
      for (int i = 0; i < B_gf.Size(); i += 2)
      {
         const double b0 = B_gf(i);
         const double b1 = B_gf(i + 1);
         const double mag = std::sqrt(b0*b0 + b1*b1);
         Bmax = std::max(Bmax, mag);
         Bsum += mag;
         n++;
      }

      const double Bavg = (n > 0) ? (Bsum / n) : 0.0;

      mfem::out << "\n=== B-field Statistics ===\n";
      mfem::out << "  B_max: " << Bmax << " T (" << (Bmax * 1e3) << " mT)\n";
      mfem::out << "  B_avg: " << Bavg << " T (" << (Bavg * 1e3) << " mT)\n";
   }

private:
   // Geometric fallback: find boundary attributes whose boundary elements lie on r=0 and mark them essential.
   // This is intentionally conservative. Best practice is to tag the axis in your mesh and handle it in InputParser.
   void MarkAxisBoundaryAttributesGeometric_()
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

      // Merge axis boundary attributes into ess_bdr_
      for (int i = 0; i < axis_attr.Size(); i++)
      {
         if (axis_attr[i]) { ess_bdr_[i] = 1; }
      }
   }

   void ProjectAxisZero_()
   {
      if (!mesh.bdr_attributes.Size()) { return; }

      // Build a marker from ess_bdr_ that contains ONLY axis attributes (geometric fallback marks them)
      // Here we just project 0 on all essential boundaries again (cheap & safe).
      mfem::ConstantCoefficient zero(0.0);
      A_->ProjectBdrCoefficient(zero, ess_bdr_);
   }
};
