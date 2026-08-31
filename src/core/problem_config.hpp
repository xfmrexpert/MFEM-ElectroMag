// Copyright (c) 2026 T. C. Raymond
// SPDX-License-Identifier: MIT

#pragma once

#include <vector>
#include <string>
#include <map>
#include <unordered_map>
#include <optional>
#include <utility>
#include "enums.hpp"      // PhysicsType, GeometryType
#include "constants.hpp"  // Constants::DEFAULT_SOLVER_*

// What computation a run performs.
//   Field          - solve the authored Scenarios as-is and save fields.
//   CouplingMatrix - auto-generate unit scenarios (drive i = 1, rest = 0) from
//                    Terminals to assemble a C (electrostatic) or L
//                    (magnetostatic) matrix. Authored Scenarios are ignored.
enum class AnalysisType { Field, CouplingMatrix };

// The across/through quantity a terminal imposes. Room to grow: Charge, Flux.
enum class Quantity { Voltage, Current };

// How the assembled linear system is solved.
//   Iterative - preconditioned Krylov (PCG/GMRES). Lowest memory; cost scales
//               with the preconditioner's convergence rate.
//   Direct    - sparse factorization. Costs more memory, but a CouplingMatrix
//               run amortizes one factorization over every terminal's RHS and
//               its accuracy does not depend on a residual tolerance.
enum class LinearSolverType { Iterative, Direct };

enum class ConductorType { Massive, Stranded };

enum class RegionCurrentConstraint { None, Open };

enum class BoundaryConditionType { Dirichlet, Neumann, Robin };

// A name bound to a set of mesh attribute ids of one entity dimension.
//
// Dim is the topological dimension of the entities the ids refer to: 1 for
// curves, 2 for surfaces, 3 for volumes. This is the mesher's own notion, not
// a role: Gmsh numbers physical groups independently per dimension, so the
// same id may name both a curve and a surface in one model, and MFEM preserves
// that split as bdr_attributes vs. attributes. An id alone is therefore
// ambiguous and the group must state which dimension it means.
//
// Whether a group is a boundary or a domain is DERIVED by comparing Dim to the
// mesh dimension; it is not authored. See IsBoundary()/IsDomain() below.
struct EntityGroup {
	int Dim = 0;                    // 1 = curve, 2 = surface, 3 = volume
	std::vector<int> AttributeIds;  // ids within the namespace selected by Dim

	// Entities of the mesh's own dimension carry the PDE and a material.
	bool IsDomain(int mesh_dim) const { return Dim == mesh_dim; }

	// Codimension-one entities bound the solved region.
	bool IsBoundary(int mesh_dim) const { return Dim == mesh_dim - 1; }
};

struct Region {
	std::string EntityGroupName;   // mesh domain (element) group name (validated)
	std::string MaterialName;      // must match a ProblemConfig::Materials key (validated)
	RegionCurrentConstraint CurrentConstraint = RegionCurrentConstraint::None;
};

// A driven/measured excitation site. Single primitive for both physics
//   DriveQuantity == Voltage -> AttributeIds are BOUNDARY attrs (essential BC)
//   DriveQuantity == Current -> AttributeIds are DOMAIN   attrs (RHS source)
struct Terminal {
    Quantity DriveQuantity = Quantity::Voltage;
    ConductorType Conductor = ConductorType::Massive;
	std::string EntityGroupName;   // mesh boundary (essential BC) or domain (RHS source) group name (validated)
};

// One scenario's setting of one terminal.
struct Excitation {
	std::string TerminalName;   // must match a Terminal::Name (validated)
	double Value = 0.0;         // volts (Voltage terminal) | amps (Current terminal)
};

struct Material {
	double Conductivity    = 0.0;   // sigma  [S/m]
	double RelPermittivity = 1.0;   // epsilon_r (vacuum default)
	double RelPermeability = 1.0;   // mu_r     (vacuum default)
};

// An authored boundary condition. For Neumann data, Value is the prescribed
// outward natural flux n.material_flux; zero is the implicit natural condition.
// Robin metadata is retained for forward compatibility but is not yet assembled.
// Terminals and magnetic axis regularity are modeled separately.
struct BoundaryCondition {
	BoundaryConditionType Type;
	std::string EntityGroupName;     // mesh boundary group name (validated)
	double Value = 0.0;
	double RobinCoeff = 0.0;         // Reserved: natural_flux + RobinCoeff*u = Value

	BoundaryCondition(BoundaryConditionType t, const std::string& g,
					  double v, double rc = 0.0)
		: Type(t), EntityGroupName(g), Value(v), RobinCoeff(rc) {}
};

// One solve: a parameter point plus per-terminal excitations.
// A terminal omitted from Excitations defaults to zero of its quantity
// (grounded for Voltage, open for Current).
struct Scenario {
	double Frequency = 0.0; // Hz; required and positive for MQS, ignored by ES/MS
	std::vector<Excitation> Excitations;
};

// Adaptive mesh refinement (AMR) controls. Parsed from the optional
// "simulation.amr" block. When Enabled is false (the default, and the case when
// the block is absent) the solver performs the legacy single solve and emits
// byte-identical output to the pre-AMR release.
struct AmrSettings {
	bool   Enabled       = false;    // Master switch.
	int    MaxIterations = 5;        // Max refine -> re-solve iterations.
	long   MaxDofs       = 2000000;  // Stop once global DOFs exceed this (<=0 disables).
	double ErrorFraction = 0.7;      // Bulk (Dorfler) marking fraction in (0, 1].
	double ErrorTolerance = 0.0;     // Absolute stop threshold on global error (<=0 ignores).
	bool   Conforming    = true;     // Require conforming output (always true for now).
};

struct ProblemConfig {
	int Order = 1;
	::PhysicsType  PhysicsType  = ::PhysicsType::Electrostatics;
	::GeometryType GeometryType = ::GeometryType::Planar;
	::AnalysisType AnalysisType = ::AnalysisType::Field;

	double SolverTolerance  = Constants::DEFAULT_SOLVER_TOLERANCE;
	int    SolverMaxIter    = Constants::DEFAULT_SOLVER_MAX_ITER;
	int    SolverPrintLevel = Constants::DEFAULT_SOLVER_PRINT_LEVEL;
	::LinearSolverType LinearSolver = ::LinearSolverType::Direct;

	std::string MeshPath;
	bool OutputParaview = false;
	bool OutputGmsh = false;
	std::string ResultsDirectory;  // Optional Gmsh results directory (empty = mesh directory)
	std::optional<int> ExportRefine;  // Export mesh refinement factor (unset = default to Order)
	AmrSettings Amr;               // Adaptive mesh refinement controls (disabled by default)
	std::unordered_map<std::string, EntityGroup> EntityGroups;
	std::vector<Region> Regions;
	std::map<std::string, Material> Materials;
	std::map<std::string, Terminal> Terminals;
	std::vector<BoundaryCondition> BoundaryConditions;
	std::vector<std::pair<std::string, Scenario>> Scenarios;
};