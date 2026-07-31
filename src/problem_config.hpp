// Copyright (c) 2026 T. C. Raymond
// SPDX-License-Identifier: MIT

#pragma once

#include <vector>
#include <string>
#include <map>
#include <unordered_map>
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

enum class EntityDim { Boundary, Domain };

enum class ConductorType { Massive, Stranded };

enum class RegionCurrentConstraint { None, Open };

struct EntityGroup {
	EntityDim Dim;              // boundary or domain
	std::vector<int> AttributeIds;   // mesh attribute ids (boundary or domain, depending on context)
};

struct Region {
	std::string EntityGroupName;   // mesh domain (element) group name (validated)
	int Material = -1;               // index into ProblemConfig::Materials
	RegionCurrentConstraint CurrentConstraint = RegionCurrentConstraint::None;
};

// A driven/measured excitation site. Single primitive for both physics
//   Excitation == Voltage -> AttributeIds are BOUNDARY attrs (essential BC)
//   Excitation == Current -> AttributeIds are DOMAIN   attrs (RHS source)
struct Terminal {
    Quantity Excitation = Quantity::Voltage;
    ConductorType Conductor = ConductorType::Massive;
	std::string EntityGroupName;   // mesh boundary (essential BC) or domain (RHS source) group name (validated)
};

// One scenario's setting of one terminal.
struct Excitation {
	std::string TerminalName;   // must match a Terminal::Name (validated)
	double Value = 0.0;         // volts (Voltage terminal) | amps (Current terminal)
	bool Floating = false;      // explicit; Value ignored. Voltage terminals only.
};

struct Material {
	double Conductivity    = 0.0;   // sigma  [S/m]
	double RelPermittivity = 1.0;   // epsilon_r (vacuum default)
	double RelPermeability = 1.0;   // mu_r     (vacuum default)
};

// A "closure" condition that makes the BVP well-posed (far-field, symmetry,
// axis regularity). Terminals are NOT modeled here.
struct BoundaryCondition {
	std::string Type;                // "Dirichlet", "Neumann", "Robin"
	std::string EntityGroupName;     // mesh boundary group name (validated)
	double Value = 0.0;
	double RobinCoeff = 0.0;         // Robin: alpha*u + beta*du/dn = value

	BoundaryCondition(const std::string& t, const std::string& g,
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

	std::string MeshPath;
	bool OutputParaview = false;
	bool OutputGmsh = false;
	std::string ResultsDirectory;  // Optional Gmsh results directory (empty = mesh directory)
	int ExportRefine = -1;         // Refinement factor for export mesh (<0 = default to Order)
	AmrSettings Amr;               // Adaptive mesh refinement controls (disabled by default)
	std::unordered_map<std::string, EntityGroup> EntityGroups;
	std::vector<Region> Regions;
	std::vector<Material> Materials;
	std::map<std::string, Terminal> Terminals;
	std::vector<BoundaryCondition> BoundaryConditions;
	std::vector<std::pair<std::string, Scenario>> Scenarios;
};