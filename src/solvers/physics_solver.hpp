// Copyright (c) 2026 T. C. Raymond
// SPDX-License-Identifier: MIT

#pragma once
#include <cmath>
#include <iomanip>
#include <limits>
#include <list>
#include <sstream>
#include "mfem.hpp"
#include "../core/problem_config.hpp"
#include "../io/gmsh_results_writer.hpp"
#include "../io/field_export.hpp"
#include "../io/matrix_writer.hpp"
#include "../io/status_reporter.hpp"
#include "amr_support.hpp"
#include "../axisym/axisymmetric_mesh.hpp"
#include "../axisym/axisymmetric_boundary_lf_integrator.hpp"
#include "../axisym/axisymmetric_curl_curl_integrator.hpp"
#include "../core/marked_boundary_condition.hpp"

/**
 * @brief Base class for physics solvers using MFEM
 *
 * @warning The mesh reference must outlive this solver instance.
 */
class PhysicsSolver {
protected:
    mfem::Mesh &mesh;
    ProblemConfig config;

    std::unique_ptr<mfem::H1_FECollection>    fec;
    std::unique_ptr<mfem::FiniteElementSpace> fespace;
    GeometryType geometry = GeometryType::Planar;
    mfem::Array<int> ess_bdr;
    mfem::Array<int> ess_tdof_list;

    // Setup-time axisymmetric mesh classification and axis boundary marker.
    // Planar problems leave it at its default.
    axisym::MeshInfo axisymmetric_mesh;
    std::vector<MarkedBoundaryCondition> closure_bcs;

    std::vector<amr::AmrIterationInfo> amr_history;

    // Non-null only while an AMR pass is in flight: the element-wise maximum of
    // the per-scenario error indicators accumulated by the current scenario loop.
    mfem::Vector* amr_errors = nullptr;

    StatusReporter& Reporter() const {
        return StatusReporter::Global();
    }

    // Publish the configuration's named entity groups into the mesh's own
    // AttributeSets containers, so a group name resolves to a marker through
    // MFEM rather than through bespoke marker-building code here.
    //
    // Each group is registered into BOTH containers because a group name alone
    // does not say whether it denotes a surface (domain) or a curve (boundary)
    // physical group, and a Gmsh model may legitimately reuse the same id for
    // one of each. Only the ids the mesh actually carries in a given container
    // are registered there, which (a) preserves the historical behaviour of
    // silently ignoring ids that do not apply to the container being queried,
    // and (b) keeps AttrToMarker's max-attribute assertion satisfied.
    void RegisterEntityGroups() {
        for (const auto& [name, group] : config.EntityGroups) {
            RegisterGroup(mesh.bdr_attribute_sets, mesh.bdr_attributes,
                          name, group.AttributeIds);
            RegisterGroup(mesh.attribute_sets, mesh.attributes,
                          name, group.AttributeIds);
        }
    }

    // Marker (1/0 over domain attributes) for a set of element attribute ids.
    // Unlike the boundary variant, a domain attribute that the mesh does not
    // carry is a configuration error: silently dropping it yields an empty
    // integration region and downstream results of exactly zero, which are much
    // harder to diagnose than a failure here.
    mfem::Array<int> DomainMarkerFromAttrs(const std::vector<int>& attrs,
                                           const std::string& context) const {
        MFEM_VERIFY(!attrs.empty(),
            "No domain attribute ids are associated with " + context + ".");
        const int max_attr = mesh.attributes.Max();
        mfem::Array<int> ids;
        ids.Reserve(static_cast<int>(attrs.size()));
        for (int a : attrs) {
            const bool is_bdr_attr = mesh.bdr_attributes.Find(a) >= 0;
            MFEM_VERIFY(a > 0 && a <= max_attr && mesh.attributes.Find(a) >= 0,
                "Domain attribute " + std::to_string(a) + " referenced by " +
                context + " does not exist in the mesh (mesh domain attributes "
                "run up to " + std::to_string(max_attr) + ")." +
                (is_bdr_attr ? " It is a boundary attribute of this mesh; the "
                               "entity group most likely names a curve physical "
                               "group instead of the surface one." : ""));
            ids.Append(a);
        }
        return mfem::AttributeSets::AttrToMarker(max_attr, ids);
    }

    // Marker (1/0 over bdr attributes) for a named entity group. Resolved through
    // the mesh's bdr_attribute_sets, populated once by RegisterEntityGroups().
    mfem::Array<int> MarkerFromGroup(const std::string& group_name) const {
        MFEM_VERIFY(config.EntityGroups.count(group_name) > 0,
            "Unknown entity group '" + group_name + "'.");
        return BoundaryMarker(group_name);
    }

    // Register only the ids that exist in @p mesh_attrs under @p name. An
    // AttributeSet holding an id the container does not know about would trip
    // AttrToMarker's max-attribute assertion at query time, and an empty set
    // cannot be queried at all (Array::Max() on an empty array), so an
    // all-foreign group is simply left unregistered; BoundaryMarker() then
    // returns an all-zero marker, matching the previous behaviour.
    static void RegisterGroup(mfem::AttributeSets& sets,
                              const mfem::Array<int>& mesh_attrs,
                              const std::string& name,
                              const std::vector<int>& attrs) {
        mfem::Array<int> present;
        present.Reserve(static_cast<int>(attrs.size()));
        for (int a : attrs) {
            if (a > 0 && mesh_attrs.Find(a) >= 0) present.Append(a);
        }
        if (present.Size() == 0) return;
        sets.SetAttributeSet(name, present);
    }

    // Boundary marker for a registered group name, or an all-zero marker when
    // the group names no boundary attribute of this mesh.
    mfem::Array<int> BoundaryMarker(const std::string& name) const {
        if (!mesh.bdr_attribute_sets.AttributeSetExists(name)) {
            mfem::Array<int> none(mesh.bdr_attributes.Max());
            none = 0;
            return none;
        }
        return mesh.bdr_attribute_sets.GetAttributeSetMarker(name);
    }

    std::vector<MarkedBoundaryCondition> BuildClosureBcs() const {
        std::vector<MarkedBoundaryCondition> bcs;
        bcs.reserve(config.BoundaryConditions.size());
        for (const auto& bc : config.BoundaryConditions) {
            MFEM_VERIFY(bc.Type != BoundaryConditionType::Robin,
                "Robin boundary conditions are reserved but not implemented. "
                "Use Dirichlet or Neumann for boundary group '" +
                bc.EntityGroupName + "'.");
            bcs.push_back({ MarkerFromGroup(bc.EntityGroupName), bc });
        }
        return bcs;
    }

    std::vector<mfem::Array<int>> DirichletClosureMarkers(
        const std::vector<MarkedBoundaryCondition>& bcs) const {
        std::vector<mfem::Array<int>> markers;
        for (const auto& bc : bcs) {
            if (bc.Condition.Type == BoundaryConditionType::Dirichlet) {
                markers.push_back(bc.Marker);
            }
        }
        return markers;
    }

    void ValidateMagneticAxisBoundaryValues() const {
        if (geometry != GeometryType::Axisymmetric ||
            !axisymmetric_mesh.TouchesAxis()) return;

        MFEM_VERIFY(fespace,
            "Magnetic axis boundary validation requires a finite element space.");

        mfem::Array<int> axis_tdofs;
        fespace->GetEssentialTrueDofs(
            axisymmetric_mesh.axis_boundary, axis_tdofs);
        mfem::Array<int> is_axis_tdof(fespace->GetTrueVSize());
        is_axis_tdof = 0;
        for (int i = 0; i < axis_tdofs.Size(); ++i) {
            is_axis_tdof[axis_tdofs[i]] = 1;
        }

        for (const auto& bc : closure_bcs) {
            if (bc.Condition.Type != BoundaryConditionType::Dirichlet ||
                bc.Condition.Value == 0.0) continue;

            mfem::Array<int> boundary_tdofs;
            fespace->GetEssentialTrueDofs(bc.Marker, boundary_tdofs);
            for (int i = 0; i < boundary_tdofs.Size(); ++i) {
                const int tdof = boundary_tdofs[i];
                MFEM_VERIFY(!is_axis_tdof[tdof],
                    "Boundary group '" + bc.Condition.EntityGroupName +
                    "' assigns a nonzero Dirichlet value at true DOF " +
                    std::to_string(tdof) + " on the magnetic symmetry axis. "
                    "Axis regularity requires A_phi = 0 at r = 0.");
            }
        }
    }

    mfem::Vector AssembleNeumannBoundaryLoad() {
        mfem::LinearForm load(fespace.get());
        std::vector<std::unique_ptr<mfem::ConstantCoefficient>> coefficients;

        for (auto& bc : closure_bcs) {
            if (bc.Condition.Type != BoundaryConditionType::Neumann ||
                bc.Condition.Value == 0.0) continue;
            coefficients.push_back(
                std::make_unique<mfem::ConstantCoefficient>(bc.Condition.Value));
            if (geometry == GeometryType::Axisymmetric) {
                load.AddBoundaryIntegrator(
                    new AxisymmetricBoundaryLFIntegrator(*coefficients.back()),
                    bc.Marker);
            } else {
                load.AddBoundaryIntegrator(
                    new mfem::BoundaryLFIntegrator(*coefficients.back()),
                    bc.Marker);
            }
        }

        load.Assemble();
        return mfem::Vector(load);
    }

    // The ordered list of (name, scenario) solves for the active analysis, so
    // every solver flows both analysis types through ONE imprint/solve loop:
    //   - Field:          the authored scenarios, as-is.
    //   - CouplingMatrix: one synthesized scenario per terminal that drives that
    //                     terminal by a unit excitation (1 V or 1 A, per physics)
    //                     while every other terminal defaults to zero. This is the
    //                     column-by-column excitation used to assemble the coupling
    //                     (capacitance / inductance) matrix.
    // Terminal order follows config.Terminals (name-sorted), which the matrix
    // assembly and WriteCouplingMatrix() also use, keeping columns aligned.
    std::vector<std::pair<std::string, Scenario>> BuildSolveScenarios() const {
        std::vector<std::pair<std::string, Scenario>> out;
        if (config.AnalysisType == AnalysisType::CouplingMatrix) {
            out.reserve(config.Terminals.size());
            for (const auto& [term_name, term] : config.Terminals) {
                Scenario sc;
                sc.Excitations.push_back({ term_name, 1.0 });
                out.push_back({ "CouplingMatrix_" + term_name, std::move(sc) });
            }
        }
        else {
            out.reserve(config.Scenarios.size());
            for (const auto& [sc_name, sc] : config.Scenarios)
                out.push_back({ sc_name, sc });
        }
        return out;
    }

    // Region/material that claims a given domain attribute, or nullptr if none
    // does. The single source of truth for attribute -> material resolution (used
    // by MaterialVector). Configuration validation rejects an attribute claimed by
    // more than one region, so the match here is unique and order-independent.
    const Material* MaterialForAttr(int attr) const {
        for (const auto& region : config.Regions) {
            const EntityGroup& group = config.EntityGroups.at(region.EntityGroupName);
            if (std::find(group.AttributeIds.begin(), group.AttributeIds.end(), attr)
                != group.AttributeIds.end())
                return &config.Materials.at(region.MaterialName);
        }
        return nullptr;
    }

    // Fill a per-domain-attribute vector from a material picker. Attributes that
    // no region claims keep `default_value`. Size = mesh.attributes.Max().
    mfem::Vector MaterialVector(double default_value,
        const std::function<double(const Material&)>& pick) const {
        const int n = mesh.attributes.Max();
        mfem::Vector v(n);
        v = default_value;
        for (int a = 1; a <= n; ++a)
            if (const Material* mat = MaterialForAttr(a))
                v[a - 1] = pick(*mat);
        return v;
    }

    // The common case: turn a material picker straight into the per-attribute
    // PWConstCoefficient the integrators consume. Use MaterialVector directly
    // only when the raw values need post-processing before wrapping.
    std::unique_ptr<mfem::PWConstCoefficient> MaterialCoefficient(
        double default_value,
        const std::function<double(const Material&)>& pick) const {
        return std::make_unique<mfem::PWConstCoefficient>(
            MaterialVector(default_value, pick));
    }

    // Essential-boundary marker: a bdr attribute is essential if ANY supplied
    // marker selects it. Unifies the (currently divergent) ess_bdr construction.
    mfem::Array<int> EssentialBdrFrom(const std::vector<mfem::Array<int>>& markers) const {
        mfem::Array<int> ess(mesh.bdr_attributes.Max());
        ess = 0;
        for (const auto& m : markers)
            for (int i = 0; i < m.Size(); ++i) if (m[i]) ess[i] = 1;
        return ess;
    }

    // Gmsh export strategy. Linear is the resample-to-linear path used today and
    // by every current consumer; HighOrder is reserved for true high-order MSH
    // output (native Lagrange elements + $InterpolationScheme) and is not wired
    // up yet.
    enum class GmshExportStrategy { Linear, HighOrder };

    // Gmsh serializer: dumb sink over a FieldExportSet. The Linear strategy
    // resamples every field onto a refined export mesh as linear views (the
    // standard tessellation approach, faithful at ref_factor ~ order). The
    // serializer computes nothing and invents no fields; magnitude/derived
    // quantities arrive as explicit FieldExport entries.
    void WriteGmshFields(const std::string& scenario_name,
                         const FieldExportSet& fields,
                         GmshExportStrategy strategy = GmshExportStrategy::Linear) const
    {
        switch (strategy) {
            case GmshExportStrategy::Linear:
                WriteGmshFieldsLinear(scenario_name, fields);
                return;
            case GmshExportStrategy::HighOrder:
                WriteGmshFieldsHighOrder(scenario_name, fields);
                return;
        }
    }

    // AMR per-iteration diagnostics from the most recent RunAdaptive(). Empty when
    // AMR is disabled. Consumed by the regression tests and useful for logging.
    const std::vector<amr::AmrIterationInfo>& GetAmrHistory() const { return amr_history; }

    // Fold the just-solved scenario's local error into the AMR indicator using an
    // element-wise maximum, so one shared mesh is refined for all scenarios.
    // Solvers call this from RunOnCurrentMesh() right after each solve; outside an
    // AMR pass it is a no-op, which is what lets a single scenario loop serve both
    // the production solve and the error estimate (no duplicate solves).
    void AccumulateScenarioError() {
        if (!amr_errors) return;

        const int ne = mesh.GetNE();
        mfem::Vector current;
        EstimateCurrentSolutionError(current);
        MFEM_VERIFY(current.Size() == ne,
            "AMR estimator returned the wrong number of element errors.");
        for (int element = 0; element < ne; ++element) {
            (*amr_errors)(element) = std::max((*amr_errors)(element), current(element));
        }
    }

    virtual void BuildOperators() = 0;
    virtual void RunOnCurrentMesh() = 0;

    // Per-element error indicator for whichever solution currently lives in the
    // solver's grid function. The scenario-wide fold lives in the base class.
    virtual void EstimateCurrentSolutionError(mfem::Vector& errors) = 0;
    virtual double ComputePeakFieldMagnitude() const = 0;

private:
    const char* PeakFieldLabel() const {
        return config.PhysicsType == PhysicsType::Electrostatics ? "|E|" : "|B|";
    }

    void RunAdaptive() {
        const AmrSettings& settings = config.Amr;
        amr_history.clear();

        const int max_iterations = std::max(1, settings.MaxIterations);
        for (int iteration = 0; iteration < max_iterations; ++iteration) {
            const long true_dofs = fespace->GetTrueVSize();

            // ONE scenario loop per AMR iteration: RunOnCurrentMesh() solves and
            // post-processes every scenario, and folds each solution's error
            // indicator into `errors` as it goes. The last iteration's results are
            // therefore already the converged-mesh results, so no extra solve pass
            // is needed after the loop.
            mfem::Vector errors(mesh.GetNE());
            errors = 0.0;
            amr_errors = &errors;
            RunOnCurrentMesh();
            amr_errors = nullptr;

            double sum_squared = 0.0;
            for (int element = 0; element < errors.Size(); ++element) {
                sum_squared += errors(element) * errors(element);
            }
            const double global_error = std::sqrt(sum_squared);
            const double peak_field = ComputePeakFieldMagnitude();
            amr_history.push_back({ true_dofs, global_error, peak_field });

            std::ostringstream diagnostic;
            diagnostic << "AMR iteration " << iteration
                << ": elements=" << mesh.GetNE()
                << ", true_dofs=" << true_dofs
                << ", global_error=" << std::scientific << std::setprecision(6)
                << global_error
                << ", peak" << PeakFieldLabel() << "=" << peak_field;
            Reporter().Diagnostic(diagnostic.str());

            if (settings.ErrorTolerance > 0.0 &&
                global_error < settings.ErrorTolerance) {
                Reporter().Diagnostic("AMR: global error below tolerance. Stop.");
                break;
            }
            if (settings.MaxDofs > 0 && true_dofs > settings.MaxDofs) {
                Reporter().Diagnostic("AMR: reached the maximum number of DOFs. Stop.");
                break;
            }
            if (iteration + 1 >= max_iterations) {
                Reporter().Diagnostic("AMR: reached the maximum number of iterations. Stop.");
                break;
            }

            mfem::Array<int> marked;
            amr::MarkElementsDorfler(errors, settings.ErrorFraction, marked);
            if (marked.Size() == 0) {
                Reporter().Diagnostic("AMR: no elements marked for refinement. Stop.");
                break;
            }

            amr::RefineConforming(mesh, marked);
            BuildOperators();
        }
    }

private:
    // Resample-to-linear Gmsh path: a refined export mesh with linear H1 (node)
    // and L2 (element-node) views. Primary scalars are sampled into H1 nodal
    // views; derived scalars/vectors are projected onto the export mesh. This
    // reproduces the historical Gmsh output (minus the auto-magnitude view,
    // which is now an explicit field when wanted).
    void WriteGmshFieldsLinear(const std::string& scenario_name,
                               const FieldExportSet& fields) const
    {
        namespace fs = std::filesystem;
        const fs::path results_dir = config.ResultsDirectory.empty()
            ? fs::path(config.MeshPath).parent_path()
            : fs::path(config.ResultsDirectory);
        const fs::path out_path = results_dir / (scenario_name + ".results.msh");

        const int ref_factor = config.ExportRefine.value_or(std::max(1, config.Order));
        mfem::Mesh export_mesh(&mesh, ref_factor, mfem::BasisType::ClosedUniform);

        const int dim = export_mesh.Dimension();
        const int sdim = export_mesh.SpaceDimension();
        MFEM_ASSERT(dim == 2, "Gmsh export assumes a 2D mesh.");

        // Linear export spaces, shared by every view of each kind. The L2 vector
        // space also defines the per-element sample points (sample_fes below).
        mfem::H1_FECollection fec_h1_lin(1, dim);
        mfem::FiniteElementSpace fes_scalar(&export_mesh, &fec_h1_lin);
        mfem::L2_FECollection fec_l2_lin(1, dim);
        mfem::FiniteElementSpace fes_vec(&export_mesh, &fec_l2_lin, sdim);

        // Resampled fields must outlive WriteGmshResults; std::list keeps element
        // addresses stable as views capture them by reference.
        std::list<mfem::GridFunction> resampled;
        std::vector<gmsh_results::View> views;

        for (const auto& f : fields.Fields()) {
            switch (f.kind) {
                case FieldExport::Kind::PrimaryScalar: {
                    resampled.emplace_back(&fes_scalar);
                    mfem::GridFunctionCoefficient sc(f.primary);
                    resampled.back().ProjectCoefficient(sc);
                    views.push_back(
                        gmsh_results::MakeScalarNodeView(f.name, resampled.back()));
                    break;
                }
                case FieldExport::Kind::DerivedScalar: {
                    resampled.emplace_back(&fes_scalar);
                    resampled.back().ProjectCoefficient(*f.scalar);
                    views.push_back(
                        gmsh_results::MakeScalarElementNodeView(f.name, resampled.back()));
                    break;
                }
                case FieldExport::Kind::DerivedVector: {
                    resampled.emplace_back(&fes_vec);
                    resampled.back() = 0.0;
                    resampled.back().ProjectCoefficient(*f.vector);
                    views.push_back(
                        gmsh_results::MakeVectorElementNodeView(f.name, resampled.back()));
                    break;
                }
            }
        }

        gmsh_results::WriteGmshResults(out_path.string(), export_mesh, fes_vec, views);
        Reporter().Diagnostic("Wrote " + out_path.string());
    }

    // True high-order MSH output (native high-order elements + per-view
    // interpolation scheme). Deferred: needs an MFEM -> Gmsh high-order node
    // ordering map. Stubbed so the strategy switch is complete and the work is
    // a single localized follow-up.
    void WriteGmshFieldsHighOrder(const std::string& scenario_name,
                                  const FieldExportSet& /*fields*/) const
    {
        Reporter().Diagnostic("WriteGmshFields(HighOrder) not implemented yet for scenario: "
                              + scenario_name + "; falling back to no output.");
    }

public:
    // ParaView serializer: dumb sink over a FieldExportSet. Primary scalars are
    // registered at native (high) order so ParaView's Lagrange cells render them
    // faithfully; derived coefficients are projected into an L2 space one order
    // below the H1 solution. Computes nothing of its own and invents no fields.
    void WriteParaviewFields(const std::string& collection_name,
                             const FieldExportSet& fields) const
    {
        const int dim = mesh.Dimension();
        MFEM_ASSERT(dim == 2, "ParaView export assumes a 2D mesh.");

        mfem::ParaViewDataCollection pv(collection_name, &mesh);

        // SetHighOrderOutput emits Lagrange cells of the solution order rather
        // than subdividing into linear cells.
        const int order = fec->GetOrder();
        pv.SetLevelsOfDetail(order);
        pv.SetHighOrderOutput(order > 1);

        // Derived fields live in L2 spaces one order below the H1 solution.
        // ParaView stores raw pointers, so the backing spaces and grid functions
        // must outlive Save(); a std::list never relocates its elements, so the
        // registered addresses stay valid as we append.
        const int l2_order = std::max(0, order - 1);
        mfem::L2_FECollection fec_l2(l2_order, dim);
        mfem::FiniteElementSpace fes_l2_scalar(&mesh, &fec_l2);
        mfem::FiniteElementSpace fes_l2_vec(&mesh, &fec_l2, dim);
        std::list<mfem::GridFunction> derived;

        for (const auto& f : fields.Fields()) {
            switch (f.kind) {
                case FieldExport::Kind::PrimaryScalar:
                    pv.RegisterField(f.name, f.primary);
                    break;
                case FieldExport::Kind::DerivedScalar:
                    derived.emplace_back(&fes_l2_scalar);
                    derived.back().ProjectCoefficient(*f.scalar);
                    pv.RegisterField(f.name, &derived.back());
                    break;
                case FieldExport::Kind::DerivedVector:
                    derived.emplace_back(&fes_l2_vec);
                    derived.back() = 0.0;
                    derived.back().ProjectCoefficient(*f.vector);
                    pv.RegisterField(f.name, &derived.back());
                    break;
            }
        }

        pv.SetCycle(0);
        pv.SetTime(0.0);
        pv.Save();

        Reporter().Diagnostic("Wrote ParaView collection " + collection_name);
    }

    // Inspect the mesh's radial extent once, at setup, for axisymmetric runs.
    // Rejects meshes reaching a materially negative radius and records whether
    // the domain closure actually touches r = 0. Annular domains (r_min > 0)
    // need no axis regularity handling at all.
    void ValidateAxisymmetricGeometry()
    {
        if (geometry != GeometryType::Axisymmetric) { return; }

        axisymmetric_mesh = axisym::ValidateMesh(mesh);

        std::ostringstream msg;
        msg << std::setprecision(6)
            << "Axisymmetric mesh radial extent: r in [" << axisymmetric_mesh.min_r
            << ", " << axisymmetric_mesh.max_r << "]; "
            << (axisymmetric_mesh.TouchesAxis()
                    ? "domain touches the symmetry axis (axis regularity enforced)."
                    : "domain is annular (no axis condition required).");
        Reporter().Diagnostic(msg.str());

        WarnOnUnderResolvedRadialQuadrature();
    }

    // The curl-curl 1/r term is integrated by a geometry-aware rule whose cost
    // is set by s = r_min/h per element (see
    // AxisymmetricCurlCurlIntegrator::RadialExtraOrder). That rule is capped, so
    // an element that is both very thin radially and very close to the axis can
    // fall outside the accuracy target. Such an element is rare and always a
    // meshing choice, but the resulting error is silent, so report it once.
    void WarnOnUnderResolvedRadialQuadrature()
    {
        int worst_element = -1;
        double worst_ratio = std::numeric_limits<double>::max();

        for (int e = 0; e < mesh.GetNE(); ++e) {
            double min_radius = 0.0;
            double radial_width = 0.0;
            AxisymmetricCurlCurlIntegrator::RadialExtent(
                *mesh.GetElementTransformation(e), min_radius, radial_width);

            // Elements meeting the axis are excluded by design: there the
            // divergent directions are removed by the A_phi = 0 constraint.
            if (!(radial_width > 0.0)) { continue; }
            if (axisymmetric_mesh.IsOnAxisGeometry(min_radius)) { continue; }

            const double ratio = min_radius / radial_width;
            if (ratio < worst_ratio) {
                worst_ratio = ratio;
                worst_element = e;
            }
        }

        if (worst_element < 0) { return; }
        if (worst_ratio >= AxisymmetricCurlCurlIntegrator::kResolvedRadiusRatio) {
            return;
        }

        std::ostringstream msg;
        msg << std::setprecision(3)
            << "Element " << worst_element << " has r_min/width = " << worst_ratio
            << ", below the ratio " << AxisymmetricCurlCurlIntegrator::kResolvedRadiusRatio
            << " at which the curl-curl 1/r quadrature reaches its accuracy "
               "target. The capped rule integrates such elements approximately; "
               "widen the innermost radial band or move it away from the axis if "
               "near-axis accuracy matters.";
        Reporter().Warning(msg.str());
    }

    double CalculateRegionArea(const std::vector<int>& attribute_ids) const {
        // AddDomainIntegrator binds the marker by non-const reference, so it
        // cannot be const here.
        mfem::Array<int> marker =
            DomainMarkerFromAttrs(attribute_ids, "a region area calculation");
        mfem::LinearForm area_form(fespace.get());
        mfem::ConstantCoefficient one(1.0);

        area_form.AddDomainIntegrator(new mfem::DomainLFIntegrator(one), marker);
        area_form.Assemble();

        double region_area = area_form.Sum();
        return region_area;
    }

public:
    PhysicsSolver(mfem::Mesh &m, const ProblemConfig &c) : mesh(m), config(c) {
        // Named groups become mesh attribute sets up front, before any Setup()
        // asks for a marker. Attribute values are refinement-invariant, so this
        // registration survives every AMR pass without being redone.
        RegisterEntityGroups();
    }

    // Virtual destructor is essential for unique_ptr polymorphism
    virtual ~PhysicsSolver() = default;

    virtual void Setup() = 0;
    void Run() {
        if (config.Amr.Enabled) {
            RunAdaptive();
        }
        else {
            amr_history.clear();
            RunOnCurrentMesh();
        }
    }
    virtual void SaveAnalysis() = 0;

    // Post-solve field recovery: each solver declares WHAT to export for the
    // just-solved scenario (primary solution fields + derived coefficients).
    // The set is format-agnostic; the writers below decide HOW to serialize it.
    virtual FieldExportSet CollectExportFields() const = 0;

    // Shared per-scenario serialization: recover the field set ONCE, then fan it
    // out to whichever formats are enabled. Solvers no longer hand-roll their own
    // ParaView/Gmsh writers.
    void SaveScenario(const std::string& scenario_name) {
        if (!config.OutputParaview && !config.OutputGmsh) return;

        FieldExportSet fields = CollectExportFields();
        if (config.OutputParaview) {
            WriteParaviewFields("results_" + std::string(ToString(config.PhysicsType))
                                + "_" + scenario_name, fields);
        }
        if (config.OutputGmsh) {
            WriteGmshFields(scenario_name, fields);
        }
    }

    // Terminal names in config (name-sorted) order. This is the same order
    // BuildSolveScenarios() drives the coupling columns in, so it labels the
    // rows/columns of any coupling matrix consistently.
    std::vector<std::string> TerminalNames() const {
        std::vector<std::string> names;
        names.reserve(config.Terminals.size());
        for (const auto& [name, term] : config.Terminals) names.push_back(name);
        return names;
    }

    // Unit label for an extracted coupling quantity.
    //
    // Axisymmetric assembly carries the full revolved measure 2*pi*r dr dz, so
    // the extracted quantity is absolute. Planar assembly integrates over the
    // (x, y) cross-section only, which is equivalent to a unit out-of-plane
    // depth: the model is translationally invariant in z and describes an
    // infinitely long structure, so the result is a per-unit-length quantity.
    // No extrusion length is configurable, so the planar label always carries
    // the "/m" suffix rather than depending on a depth setting.
    [[nodiscard]] std::string CouplingUnitLabel(const std::string& si_unit) const {
        return geometry == GeometryType::Axisymmetric
            ? "[" + si_unit + "]"
            : "[" + si_unit + "/m]";
    }

    // Shared coupling-matrix serialization: print a labeled table to the console
    // and write a CSV next to the mesh (same path convention as the field
    // writers). Each solver supplies the assembled matrix plus the human-readable
    // title (with unit) and CSV file name; rows/columns are labeled by terminal.
    void SaveCouplingMatrix(const mfem::DenseMatrix& M,
                            const std::string& title,
                            const std::string& csv_filename) const {
        matrix_io::MatrixWriter writer(title, TerminalNames());
        if (Reporter().IsMachineReadable()) {
            std::ostringstream table;
            writer.PrintConsole(M, table);
            Reporter().Status(table.str());
        }
        else {
            writer.PrintConsole(M);
        }

        namespace fs = std::filesystem;
        fs::path out_path = fs::path(config.MeshPath).parent_path() / csv_filename;
        writer.WriteCsv(M, out_path);
    }

};