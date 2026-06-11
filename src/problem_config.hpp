// Copyright (c) 2026 T. C. Raymond
// SPDX-License-Identifier: MIT

#pragma once

#include <vector>
#include <string>
#include "enums.hpp"      // ModelType
#include "constants.hpp"  // Constants::DEFAULT_SOLVER_*

// What deliverable a run produces.
//   Field          - solve the authored Scenarios as-is and save fields.
//   CouplingMatrix - auto-generate unit scenarios (drive i = 1, rest = 0) from
//                    Terminals to assemble a C (electrostatic) or L
//                    (magnetostatic) matrix. Authored Scenarios are ignored.
enum class StudyType { Field, CouplingMatrix };

// The across/through quantity a terminal imposes. Room to grow: Charge, Flux.
enum class Quantity { Voltage, Current };

struct Region {
	std::vector<int> AttributeIds;   // mesh domain (element) attribute ids
	int Material = -1;               // index into ProblemConfig::Materials
};

// A driven/measured excitation site. Single primitive for both physics
//   Excitation == Voltage -> AttributeIds are BOUNDARY attrs (essential BC)
//   Excitation == Current -> AttributeIds are DOMAIN   attrs (RHS source)
struct Terminal {
    std::string Name;
    Quantity Excitation = Quantity::Voltage;
	std::vector<int> AttributeIds;
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
	std::vector<int> AttributeIds;   // mesh boundary attribute ids
	double Value = 0.0;
	double RobinCoeff = 0.0;         // Robin: alpha*u + beta*du/dn = value

	BoundaryCondition(const std::string& t, const std::vector<int>& m,
					  double v, double rc = 0.0)
		: Type(t), AttributeIds(m), Value(v), RobinCoeff(rc) {}
};

// One solve: a parameter point (Frequency, ...) plus per-terminal excitations.
// A terminal omitted from Excitations defaults to zero of its quantity
// (grounded for Voltage, open for Current).
struct Scenario {
	std::string Name;
	std::vector<Excitation> Excitations;
};

struct ProblemConfig {
	int Order = 1;
	::ModelType ModelType = ::ModelType::Planar;
	::StudyType StudyType = ::StudyType::Field;
	double Frequency = 60.0;     // MQS only; constant across the study (ignored by ES/MS)

	double SolverTolerance  = Constants::DEFAULT_SOLVER_TOLERANCE;
	int    SolverMaxIter    = Constants::DEFAULT_SOLVER_MAX_ITER;
	int    SolverPrintLevel = Constants::DEFAULT_SOLVER_PRINT_LEVEL;

	std::string MeshPath;
	bool OutputParaview = false;
	bool OutputGmsh = false;
	std::string ResultsFile;       // Optional Gmsh MSH 2.2 results path (empty = derive from config)
	int ExportRefine = -1;         // Refinement factor for export mesh (<0 = default to Order)
	std::vector<Region> Regions;
	std::vector<Material> Materials;
	std::vector<Terminal> Terminals;
	std::vector<BoundaryCondition> BoundaryConditions;
	std::vector<Scenario> Scenarios;
};