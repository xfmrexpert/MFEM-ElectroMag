// Copyright (c) 2026 T. C. Raymond
// SPDX-License-Identifier: MIT

#pragma once
#include <cmath>
#include <iomanip>
#include <limits>
#include <sstream>
#include "mfem.hpp"
#include "../core/problem_config.hpp"
#include "../io/field_export.hpp"
#include "../io/matrix_writer.hpp"
#include "../io/solver_field_writer.hpp"
#include "../io/status_reporter.hpp"
#include "amr_support.hpp"
#include "../axisym/axisymmetric_mesh_validation.hpp"
#include "../axisym/axisymmetric_boundary_lf_integrator.hpp"
#include "../core/marked_boundary_condition.hpp"

/**
 * @brief Base class for physics solvers using MFEM
 *
 * @warning The mesh reference must outlive this solver instance.
 */
class PhysicsSolver {
    // ---- State --------------------------------------------------------------
protected:
	mfem::Mesh& mesh;
	ProblemConfig config;

    // Held as the BASE collection type deliberately. Every solver currently
    // builds an H1 space, but nothing in this class requires H1: the FE space
    // constructor and GetOrder() are both base-class API. Naming the concrete
    // type here would commit every present and future solver to a nodal scalar
    // discretization, which is a formulation choice that belongs to the derived
    // solver, not to the shared plumbing.
    std::unique_ptr<mfem::FiniteElementCollection> fec;
    std::unique_ptr<mfem::FiniteElementSpace> fespace;
    GeometryType geometry = GeometryType::Planar;
    mfem::Array<int> ess_bdr;
    mfem::Array<int> ess_tdof_list;

    // The prescribed boundary conditions of this solve
    // Terminals are NOT here; each solver realizes those itself. See
    // docs/boundary_and_terminal_model.md.
    BoundaryConditionSet boundary_conditions;
private:
    std::vector<amr::AmrIterationInfo> amr_history;
	// Non-null only while an AMR pass is in flight: the running root-mean-square
	// of the per-scenario error indicators accumulated by the current scenario
	// loop. Owned by RunAdaptive(); derived solvers reach it only through
	// AccumulateScenarioError().
	mfem::Vector* amr_errors = nullptr;

	// Number of scenario indicators folded into *amr_errors so far. Reset per
	// AMR iteration alongside amr_errors; drives the running-average update.
	int amr_error_samples = 0;

    // ---- Public API ---------------------------------------------------------
public:
    PhysicsSolver(mfem::Mesh& m, const ProblemConfig& c) : mesh(m), config(c) {
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
    // The set is format-agnostic; SolverFieldWriter decides HOW to serialize it.
    virtual FieldExportSet CollectExportFields() const = 0;

    // AMR per-iteration diagnostics from the most recent Run(). Empty when AMR
    // is disabled. Read-only convergence history for logging and regression tests.
    const std::vector<amr::AmrIterationInfo>& GetAmrHistory() const { return amr_history; }

    // ---- Virtual methods: each solver supplies its own physics -------------
protected:
    virtual void BuildOperators() = 0;
    virtual void RunOnCurrentMesh() = 0;

    // Per-element error indicator for whichever solution currently lives in the
    // solver's grid function. The scenario-wide fold lives in the base class.
    virtual void EstimateCurrentSolutionError(mfem::Vector& errors) = 0;
    virtual double ComputePeakFieldMagnitude() const = 0;

    // Builds ess_bdr: every boundary attribute whose DOFs this formulation pins.
    //
    // The default is the prescribed Dirichlet conditions alone. Formulations with
    // additional essential sources override this, call the base, and merge their
    // own markers in -- voltage terminals for electrostatics, axis regularity for
    // the axisymmetric magnetic solvers. Collecting the union behind one named
    // hook is what makes "what pins DOFs here?" answerable per solver, and what
    // guarantees the merges happen before BuildOperators() consumes ess_bdr.
    //
    // Called from Setup(); the result is refinement-invariant (attribute-keyed)
    // and is reused across every AMR pass.
    virtual void BuildEssentialBoundaryMarker() {
        ess_bdr = boundary_conditions.DirichletMarker(mesh.bdr_attributes.Max());
    }

    // ---- Shared helpers for derived solvers ---------------------------------

    StatusReporter& Reporter() const {
        return StatusReporter::Global();
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

    // Region/material that claims a given domain attribute, or nullptr if none
    // does. The single source of truth for attribute -> material resolution.
    // Configuration validation rejects an attribute claimed by more than one
    // region, so the match here is unique and order-independent.
    const Material* MaterialForAttr(int attr) const {
        for (const auto& region : config.Regions) {
            const EntityGroup& group = config.EntityGroups.at(region.EntityGroupName);
            if (std::find(group.AttributeIds.begin(), group.AttributeIds.end(), attr)
                != group.AttributeIds.end())
                return &config.Materials.at(region.MaterialName);
        }
        return nullptr;
    }

    // Per-domain-attribute material property for the integrators. Attributes
    // that no region claims fall back to `default_value`; @p pick reads the
    // wanted quantity off whichever material a region assigns.
    std::unique_ptr<mfem::PWConstCoefficient> MaterialCoefficient(
        double default_value,
        const std::function<double(const Material&)>& pick) const {
        mfem::Vector values = MaterialVector(default_value, pick);
        return std::make_unique<mfem::PWConstCoefficient>(values);
    }

    // Per-domain-attribute vector laid out as PWConstCoefficient expects
    // (attribute a -> element a-1): `value` on the listed attributes, zero
    // elsewhere. The companion to MaterialCoefficient() for quantities that
    // come from the configuration rather than from a material.
    mfem::Vector AttributeVector(const std::vector<int>& attrs,
                                 double value) const {
        mfem::Vector v(mesh.attributes.Max());
        v = 0.0;
        for (int attr : attrs) {
            if (attr > 0 && attr <= v.Size()) { v[attr - 1] = value; }
        }
        return v;
    }

    // The drive a scenario applies to a terminal, or 0.0 when the scenario does
    // not mention it. Excitations are prescribed as a list rather than a map, so
    // this is the single place that resolves one against a terminal name.
    //
    // The value is returned unscaled. For time-harmonic solvers it is a PEAK
    // (amplitude) phasor by convention; this function is convention-agnostic
    // and performs no rms/peak conversion. See Excitation in problem_config.hpp.
    static double ExcitationFor(const Scenario& sc,
                                const std::string& terminal_name) {
        double value = 0.0;
        for (const auto& exc : sc.Excitations) {
            if (exc.TerminalName == terminal_name) { value = exc.Value; }
        }
        return value;
    }

    // The fixed nonzero Dirichlet values. Solvers re-imprint these on every
    // scenario (the values are lifted into the RHS by FormLinearSystem, so they
    // must be re-applied after each solution reset), and this is the one place
    // that decides which conditions qualify. Terminal drives are not included:
    // they are scenario-dependent and are projected from the excitation list.
    template <typename Fn>
    void ForEachNonzeroDirichlet(Fn&& apply) {
        for (const auto& bc : boundary_conditions) {
            if (bc.IsNonzeroDirichlet()) {
                mfem::Array<int> marker(bc.Marker);
                apply(marker, bc.Condition.Value);
            }
        }
    }

    // Validate the mesh as (r,z) input, once, at setup. Delegates to
    // axisym::ValidateMesh, which ABORTS on a materially negative radius -- an
    // invalid coordinate system for any formulation, and one that would
    // otherwise produce a negative axisymmetric measure and a silently
    // indefinite operator.
    //
    // The scan is RETURNED rather than cached here: nothing in a general
    // physics solve consumes the radial extent or the axis tolerance. Both are
    // read only by the A_phi formulations, which hold onto the result
    // themselves. Returns a default-constructed value for planar runs.
    axisym::AxisGeometry ValidateAxisymmetricGeometry()
    {
        if (geometry != GeometryType::Axisymmetric) { return {}; }

        const axisym::AxisGeometry info = axisym::ValidateMesh(mesh);

        std::ostringstream msg;
        msg << std::setprecision(6)
            << "Axisymmetric mesh radial extent: r in ["
            << info.min_r << ", " << info.max_r << "].";
        Reporter().Diagnostic(msg.str());

        return info;
    }

    BoundaryConditionSet BuildBoundaryConditions() const {
        BoundaryConditionSet bcs;
        for (const auto& bc : config.BoundaryConditions) {
            MFEM_VERIFY(bc.Type != BoundaryConditionType::Robin,
                "Robin boundary conditions are reserved but not implemented. "
                "Use Dirichlet or Neumann for boundary group '" +
                bc.EntityGroupName + "'.");
            bcs.Add(MarkerFromGroup(bc.EntityGroupName), bc);
        }
        return bcs;
    }

    mfem::Vector AssembleNeumannBoundaryLoad() {
        mfem::LinearForm load(fespace.get());
        std::vector<std::unique_ptr<mfem::ConstantCoefficient>> coefficients;
        // MFEM binds the marker by non-const reference and keeps the pointer, so
        // these copies must stay alive until Assemble() has run.
        std::vector<std::unique_ptr<mfem::Array<int>>> markers;

        for (const auto& bc : boundary_conditions) {
            if (!bc.IsNeumann() || bc.Condition.Value == 0.0) continue;
            coefficients.push_back(
                std::make_unique<mfem::ConstantCoefficient>(bc.Condition.Value));
            markers.push_back(std::make_unique<mfem::Array<int>>(bc.Marker));
            if (geometry == GeometryType::Axisymmetric) {
                load.AddBoundaryIntegrator(
                    new AxisymmetricBoundaryLFIntegrator(*coefficients.back()),
                    *markers.back());
            } else {
                load.AddBoundaryIntegrator(
                    new mfem::BoundaryLFIntegrator(*coefficients.back()),
                    *markers.back());
            }
        }

        load.Assemble();
        return mfem::Vector(load);
    }

    // The ordered list of (name, scenario) solves for the active analysis, so
    // every solver flows both analysis types through ONE imprint/solve loop:
    //   - Field:          the prescribed scenarios, as-is.
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

    // Fold the just-solved scenario's local error indicator into the running
    // combination, so one shared mesh is refined for all scenarios.
    // Solvers call this from RunOnCurrentMesh() right after each solve; outside an
    // AMR pass it is a no-op, which is what lets a single scenario loop serve both
    // the production solve and the error estimate (no duplicate solves).
    //
    // The fold is a running ROOT-MEAN-SQUARE over the scenarios seen so far:
    //     e_k = sqrt( (1/N) * sum_s eta_k(s)^2 ).
    // Updated incrementally from the previous average and the sample count, so
    // no per-scenario history is retained:
    //     e_k <- sqrt( (e_k^2 * n + eta_k^2) / (n + 1) ).
    //
    // RMS rather than the element-wise maximum: the max lets a single outlier
    // scenario dictate the refinement pattern, so elements are spent resolving a
    // feature only one excitation cares about. Averaging instead accounts for
    // error that is moderately important to MANY solves. This matters most for
    // coupling-matrix runs, where one scenario is synthesized per terminal and
    // every column is equally part of the answer.
    //
    // Each scenario's indicator is normalized by that scenario's field energy
    // (see EstimateCurrentSolutionError implementations) so the samples entering
    // this average are dimensionless relative errors and therefore comparable.
    void AccumulateScenarioError() {
        if (!amr_errors) return;

        const int ne = mesh.GetNE();
        mfem::Vector current;
        EstimateCurrentSolutionError(current);
        MFEM_VERIFY(current.Size() == ne,
            "AMR estimator returned the wrong number of element errors.");

        const double n = static_cast<double>(amr_error_samples);
        for (int element = 0; element < ne; ++element) {
            const double prev = (*amr_errors)(element);
            const double eta = current(element);
            (*amr_errors)(element) =
                std::sqrt((prev * prev * n + eta * eta) / (n + 1.0));
        }
        ++amr_error_samples;
    }

    // Integral of 1 over the given domain attributes, i.e. the measure of that
    // region in the MESH's own dimension: length in 1D, area in 2D, volume in 3D.
    //
    // Note what this deliberately is NOT. It uses a plain DomainLFIntegrator, so
    // it carries no axisymmetric r-weight: on an axisymmetric mesh it returns the
    // r-z cross-section area, not the revolved volume 2*pi*Int(r dA). That is the
    // right quantity for the current caller (a conductor cross-section feeding
    // J = I/A) but it is the wrong quantity for anything that wants a true
    // physical volume, which must integrate the axisymmetric weight instead.
    double CalculateRegionMeasure(const std::vector<int>& attribute_ids) const {
        // AddDomainIntegrator binds the marker by non-const reference, so it
        // cannot be const here.
        mfem::Array<int> marker =
            DomainMarkerFromAttrs(attribute_ids, "a region measure calculation");
        mfem::LinearForm measure_form(fespace.get());
        mfem::ConstantCoefficient one(1.0);

        measure_form.AddDomainIntegrator(new mfem::DomainLFIntegrator(one), marker);
        measure_form.Assemble();

        return measure_form.Sum();
    }

    // Shared per-scenario serialization: recover the field set ONCE, then fan it
    // out to whichever formats are enabled. The writer owns the format details;
    // solvers only declare WHAT to export via CollectExportFields().
    void SaveScenario(const std::string& scenario_name) {
        if (!config.OutputParaview && !config.OutputGmsh) return;

        FieldExportSet fields = CollectExportFields();
        const SolverFieldWriter writer(mesh, config.ResultsDirectory,
                                       config.MeshPath, fec->GetOrder());
        if (config.OutputParaview) {
            writer.WriteParaview("results_" + std::string(ToString(config.PhysicsType))
                + "_" + scenario_name, fields);
        }
        if (config.OutputGmsh) {
            writer.WriteGmsh(scenario_name, fields);
        }
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

    // ---- Base-class internals -----------------------------------------------
private:

    // Publish the configuration's named entity groups into the mesh's own
    // AttributeSets containers, so a group name resolves to a marker through
    // MFEM rather than through bespoke marker-building code here. Called once
    // from the constructor; attribute values are refinement-invariant, so no
    // derived solver ever needs to repeat it.
    //
    // Each group is registered into exactly ONE container, chosen by the entity
    // dimension it declares. Registering into both would defeat the purpose of
    // that declaration: Gmsh numbers physical groups independently per
    // dimension, so a model may legitimately use the same id for a curve and a
    // surface, and a group registered into both namespaces would silently
    // resolve to whichever unrelated entity happens to share its number.
    void RegisterEntityGroups() {
        const int mesh_dim = mesh.Dimension();
        for (const auto& [name, group] : config.EntityGroups) {
            MFEM_VERIFY(group.IsDomain(mesh_dim) || group.IsBoundary(mesh_dim),
                "Entity group '" + name + "' declares dim " +
                std::to_string(group.Dim) + ", which is neither the mesh "
                "dimension (" + std::to_string(mesh_dim) + ", a domain) nor one "
                "below it (" + std::to_string(mesh_dim - 1) + ", a boundary).");
            if (group.IsBoundary(mesh_dim)) {
                RegisterGroup(mesh.bdr_attribute_sets, mesh.bdr_attributes,
                              name, group.AttributeIds);
            }
            else {
                RegisterGroup(mesh.attribute_sets, mesh.attributes,
                              name, group.AttributeIds);
            }
        }
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

    // Unchecked marker lookup behind MarkerFromGroup(): returns an all-zero
    // marker when the group names no boundary attribute of this mesh. Solvers
    // go through MarkerFromGroup(), which also rejects unknown group names.
    mfem::Array<int> BoundaryMarker(const std::string& name) const {
        if (!mesh.bdr_attribute_sets.AttributeSetExists(name)) {
            mfem::Array<int> none(mesh.bdr_attributes.Max());
            none = 0;
            return none;
        }
        return mesh.bdr_attribute_sets.GetAttributeSetMarker(name);
    }

    // Per-domain-attribute values behind MaterialCoefficient(), indexed as
    // PWConstCoefficient expects (attribute a -> element a-1). Attributes that
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

    // Terminal names in config (name-sorted) order. This is the same order
    // BuildSolveScenarios() drives the coupling columns in, so it labels the
    // rows/columns of any coupling matrix consistently.
    std::vector<std::string> TerminalNames() const {
        std::vector<std::string> names;
        names.reserve(config.Terminals.size());
        for (const auto& [name, term] : config.Terminals) names.push_back(name);
        return names;
    }

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
            amr_error_samples = 0;
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
};
