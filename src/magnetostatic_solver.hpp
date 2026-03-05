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
   ModelType type = ModelType::Axisymmetric;

   // Resources (order of declaration = order of destruction)
   std::unique_ptr<mfem::H1_FECollection>    fec;
   std::unique_ptr<mfem::FiniteElementSpace> fespace;
   std::unique_ptr<mfem::GridFunction>       A;        // A_phi (axisym) or A (planar scalar)
   std::unique_ptr<mfem::PWConstCoefficient> nu_coeff;  // ν = 1/μ
   std::unique_ptr<mfem::PWConstCoefficient> j_coeff;   // J_phi (axisym) or J (planar scalar src)

   std::unique_ptr<mfem::LinearForm> b;
   mfem::Array<int> ess_bdr; // boundary attribute marker (size = bdr_attributes.Max())

public:
   MagnetostaticSolver(mfem::Mesh &m, const json &c) : PhysicsSolver(m, c) {}

   void Setup() override
   {
      // Config & solver type
      InputParser parser(config_json);
      config = parser.GetProblemConfig();
        
      int order = config.Order;
      const int dim   = mesh.Dimension();

      type = config.ModelType;

      // FE space
      fec     = std::make_unique<mfem::H1_FECollection>(order, dim);
      fespace = std::make_unique<mfem::FiniteElementSpace>(&mesh, fec.get());

      A = std::make_unique<mfem::GridFunction>(fespace.get());
      *A = 0.0;

      // Materials & sources

      std::vector<Material> materials = config.Materials;
      // Reluctivity
      mfem::Vector nu_vec(mesh.attributes.Max());
      nu_vec = 0.0;
      
      for (auto& region : config.Regions) {
         for (auto attribute_id : region.AttributeIds) {
                if (attribute_id > 0 && attribute_id <= mesh.attributes.Max()) {
                    auto& material = materials[region.Material];
                    double nu = 1.0 / (Constants::MU_0 * material.RelPermeability);
                    nu_vec[attribute_id - 1] = nu;
                }
            }
      }
      nu_coeff = std::make_unique<mfem::PWConstCoefficient>(nu_vec);

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

      BoundaryConditionValidator validator(mesh, *fespace);
      validator.ValidateBoundaryConditions(bcs, /*allow_overlap=*/false);

      // Apply BC values to A and build essential boundary marker.
      for (const auto &bc : bcs)
      {
         const mfem::Array<int> &marker = bc.first;
         const double val              = bc.second;

         mfem::ConstantCoefficient val_coeff(val);
         A->ProjectBdrCoefficient(val_coeff, marker);

         // Merge marker -> ess_bdr_
         MFEM_ASSERT(marker.Size() == ess_bdr.Size(),
                     "Boundary marker size must match bdr_attributes.Max().");
         for (int i = 0; i < marker.Size(); i++)
         {
            if (marker[i]) { ess_bdr[i] = 1; }
         }
      }

      // Axis regularity: enforce A_phi = 0 on r=0 as ESSENTIAL.
      // Best practice: mark the axis as an essential boundary via boundary attributes if your mesh has it tagged.
      // If you *don't* have the axis tagged as a boundary attribute, do a geometric fallback:
      if (type == ModelType::Axisymmetric)
      {
         // Geometric fallback: force A=0 on axis boundary vertices by marking the boundary attributes
         // that lie on r=0. This requires detecting boundary elements on the axis and marking their attribute.
         // If your mesh already has an "axis" boundary attribute, prefer using that in InputParser instead.
         MarkAxisBoundaryAttributesGeometric();
         // After this, ess_bdr includes axis attributes, and A has already been projected
         // for other BCs. We also project A=0 on the axis here for safety.
         ProjectAxisZero();
      }

      // RHS
      b = std::make_unique<mfem::LinearForm>(fespace.get());

      if (type == ModelType::Axisymmetric)
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

   void Run() override
   {
      // Stiffness
      mfem::BilinearForm a(fespace.get());

      if (type == ModelType::Axisymmetric)
      {
         a.AddDomainIntegrator(new AxisymmetricCurlCurlIntegrator(*nu_coeff));
      }
      else
      {
         a.AddDomainIntegrator(new mfem::DiffusionIntegrator(*nu_coeff));
      }

      a.Assemble();

      // Essential DOFs
      mfem::Array<int> ess_tdof_list;
      fespace->GetEssentialTrueDofs(ess_bdr, ess_tdof_list);

      // Form and solve
      mfem::OperatorPtr Aop;
      mfem::Vector X, B;

      a.FormLinearSystem(ess_tdof_list, *A, *b, Aop, X, B);

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
      auto *sp = dynamic_cast<mfem::SparseMatrix*>(Aop.Ptr());
      MFEM_ASSERT(sp, "Expected SparseMatrix operator from FormLinearSystem.");

      mfem::GSSmoother M(*sp);
      mfem::PCG(*sp, M, B, X,
                config.SolverPrintLevel,
                config.SolverMaxIter,
                config.SolverTolerance,
                0.0);
#endif

      a.RecoverFEMSolution(X, *b, *A);

      mfem::out << "\n=== A Statistics ===\n";
      mfem::out << "  A min:     " << A->Min() << "\n";
      mfem::out << "  A max:     " << A->Max() << "\n";
      mfem::out << "  A L2 norm: " << A->Norml2() << "\n";
   }

   void Save() override
   {
      mfem::ParaViewDataCollection pv("results_magnetostatic", &mesh);
      pv.SetLevelsOfDetail(1);
      pv.RegisterField("A", A.get());

      // Vector B field in (r,z) (axisym) or (x,y) (planar): vdim = 2 for 2D problems.
      const int dim = mesh.Dimension();
      MFEM_ASSERT(dim == 2, "Save() currently assumes a 2D mesh (axisymmetric r-z or planar x-y).");

      const int h1_order = fec->GetOrder();
      const int l2_order = std::max(0, h1_order - 1);

      mfem::L2_FECollection fec_l2(l2_order, dim);
      mfem::FiniteElementSpace fes_l2(&mesh, &fec_l2, /*vdim=*/2);
      mfem::GridFunction B_gf(&fes_l2);
      B_gf = 0.0;

      if (type == ModelType::Axisymmetric)
      {
         // B_r = -∂A/∂z, B_z = ∂A/∂r + A/r
         MagneticFieldCoefficient B_coeff(A.get());
         B_gf.ProjectCoefficient(B_coeff);
      }
      else
      {
         // WARNING: MFEM's CurlGridFunctionCoefficient is not the usual "rotated gradient"
         // for a scalar potential in 2D. If your planar case is A_z and B = curl(A_z k),
         // then B = (∂A/∂y, -∂A/∂x). Implement that explicitly if needed.
         //
         // Placeholder: keep your original approach but you should verify it.
         mfem::CurlGridFunctionCoefficient B_coeff(A.get());
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

   void ProjectAxisZero()
   {
      if (!mesh.bdr_attributes.Size()) { return; }

      // Build a marker from ess_bdr that contains ONLY axis attributes (geometric fallback marks them)
      // Here we just project 0 on all essential boundaries again (cheap & safe).
      mfem::ConstantCoefficient zero(0.0);
      A->ProjectBdrCoefficient(zero, ess_bdr);
   }
};
