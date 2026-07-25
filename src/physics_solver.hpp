// Copyright (c) 2026 T. C. Raymond
// SPDX-License-Identifier: MIT

#pragma once
#include <list>
#include "mfem.hpp"
#include "json.hpp"
#include "problem_config.hpp"
#include "gmsh_results_writer.hpp"
#include "field_export.hpp"
#include "matrix_writer.hpp"
#include "status_reporter.hpp"

using json = nlohmann::json;

/**
 * @brief Base class for physics solvers using MFEM
 *
 * @warning The mesh and config references must outlive this solver instance.
 *          Do not destroy the mesh or config objects before the solver is done.
 */
class PhysicsSolver {
protected:
    mfem::Mesh &mesh;
    const json &config_json;
    ProblemConfig config;

    std::unique_ptr<mfem::H1_FECollection>    fec;
    std::unique_ptr<mfem::FiniteElementSpace> fespace;
    GeometryType geometry = GeometryType::Planar;
    mfem::Array<int> ess_bdr;
    mfem::Array<int> ess_tdof_list;

    StatusReporter& Reporter() const {
        return StatusReporter::Global();
    }

    mfem::Array<int> MarkerFromAttrs(const std::vector<int>& attrs) const {
        mfem::Array<int> m(mesh.bdr_attributes.Max());
        m = 0;
        for (int a : attrs) {
            if (a > 0 && a <= m.Size()) m[a - 1] = 1;
        }
        return m;
    }

    // Marker (1/0 over bdr attributes) for a named entity group.
    mfem::Array<int> MarkerFromGroup(const std::string& group_name) const {
        return MarkerFromAttrs(config.EntityGroups.at(group_name).AttributeIds);
    }

    // Closure boundary conditions as (marker, value) pairs, in config order.
    // Shared by ess_bdr construction and BoundaryConditionValidator.
    std::vector<std::pair<mfem::Array<int>, double>> BuildClosureBcs() const {
        std::vector<std::pair<mfem::Array<int>, double>> bcs;
        bcs.reserve(config.BoundaryConditions.size());
        for (const auto& bc : config.BoundaryConditions)
            bcs.push_back({ MarkerFromGroup(bc.EntityGroupName), bc.Value });
        return bcs;
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

    // First region/material that claims a given domain attribute, or nullptr if
    // none does. The single source of truth for attribute -> material resolution
    // (used by MaterialVector and TerminalConductivity).
    const Material* MaterialForAttr(int attr) const {
        for (const auto& region : config.Regions) {
            const EntityGroup& group = config.EntityGroups.at(region.EntityGroupName);
            if (std::find(group.AttributeIds.begin(), group.AttributeIds.end(), attr)
                != group.AttributeIds.end())
                return &config.Materials[region.Material];
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

        const int ref_factor = (config.ExportRefine > 0) ? config.ExportRefine
            : std::max(1, config.Order);
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
            mfem::Element* bEl = mesh.GetBdrElement(be);
            const int attr = bEl->GetAttribute();
            mfem::Array<int> v;
            bEl->GetVertices(v);

            bool on_axis = true;
            for (int i = 0; i < v.Size(); i++)
            {
                const double* vx = mesh.GetVertex(v[i]);
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

    double CalculateRegionArea(const std::vector<int>& attribute_ids) const {
        mfem::Array<int> marker(mesh.attributes.Max());
        marker = 0;
        for (int attr : attribute_ids) {
            if (attr > 0 && attr <= marker.Size()) marker[attr - 1] = 1;
        }
        mfem::LinearForm area_form(fespace.get());
        mfem::ConstantCoefficient one(1.0);

        area_form.AddDomainIntegrator(new mfem::DomainLFIntegrator(one), marker);
        area_form.Assemble();

        double region_area = area_form.Sum();
        return region_area;
    }

public:
    PhysicsSolver(mfem::Mesh &m, const json &c) : mesh(m), config_json(c) {}

    // Virtual destructor is essential for unique_ptr polymorphism
    virtual ~PhysicsSolver() = default;

    virtual void Setup() = 0;
    virtual void Run() = 0;
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