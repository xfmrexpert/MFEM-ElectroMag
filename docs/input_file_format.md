# Input File Format

The solver reads one JSON configuration file. This document defines every
section, every term, and the rules a valid file must satisfy.

Companion document: [boundary_and_terminal_model.md](boundary_and_terminal_model.md)
covers the boundary/terminal/coupling vocabulary in depth and records why those
terms were chosen. This document is the schema reference; that one is the
rationale.

---

## 1. The model in one page

A configuration answers five questions, in this order. The order matters,
because each layer only refers to names defined by the layer above it.

| # | Question | Section | Refers to |
|---|----------|---------|-----------|
| 1 | How is the run configured? | `simulation` | - |
| 2 | Which parts of the mesh have names? | `entity_groups` | mesh attribute ids |
| 3 | What are those parts made of? | `materials`, `regions` | entity groups |
| 4 | How is the model constrained and driven? | `boundary_conditions`, `terminals` | entity groups |
| 5 | What is actually solved? | `scenarios` | terminals |

The mesh is the substrate. Everything else is a layer of naming and physics
placed on top of it, and nothing in the config refers to raw mesh numbers except
`entity_groups`.

```
mesh attribute ids
		|
		v
  entity_groups ------------------+--------------------+
		|                         |                    |
		v                         v                    v
	regions --> materials  boundary_conditions     terminals
														|
														v
												  scenarios
												 (excitations)
```

---

## 2. Glossary

These are the load-bearing terms. Each is defined once here and used
consistently everywhere else.

| Term | Definition |
|------|------------|
| **Attribute id** | An integer tag Gmsh/MFEM attaches to mesh elements. Domain attributes tag area/volume elements; boundary attributes tag curve/surface elements. The two numbering spaces are independent -- attribute `1` as a domain is unrelated to attribute `1` as a boundary. |
| **Entity group** | A *name* bound to a set of attribute ids of one entity dimension. The sole bridge between mesh numbering and physics. Every other section refers to groups by name, never to raw ids. |
| **Domain** | A full-dimensional part of the model (an area in 2D). Carries a material and the PDE. |
| **Boundary** | A codimension-one part (a curve in 2D) bounding the solved region. |
| **Material** | A named set of physical properties: `sigma`, `epsilon_r`, `mu_r`. Purely a property bundle; it has no location. |
| **Region** | The assignment of a material to a domain entity group. This is what gives a material a location. |
| **Boundary condition** | A condition imposed on a boundary entity group: Dirichlet (prescribed value) or Neumann (prescribed flux). Authored once; fixed for the entire run. |
| **Terminal** | A named connection through which the model is driven or measured. Its value is *not* fixed here -- scenarios supply it. |
| **Excitation** | One scenario's setting of one terminal. Volts for a voltage terminal, amps for a current terminal. |
| **Scenario** | One solve: a set of excitations plus, for MQS, a frequency. |
| **Coupling matrix** | The terminal-by-terminal result extracted by driving each terminal in turn (capacitance, inductance, or R/L). |

### Distinctions that matter

**Material vs. region.** A material is *what something is made of*; a region is
*where that material lives*. Two regions may share one material. This split is
why materials can be reused across a model without duplication.

**Boundary condition vs. terminal.** Both can pin DOFs, so the distinction is
not "essential vs. not". It is **who owns the value**:

- a boundary condition's value is authored in the file and fixed for the run;
- a terminal's value is supplied per scenario.

A grounded far-field boundary is a boundary condition. A driven electrode is a
terminal. See the companion document for the full argument.

**Terminal vs. excitation.** The terminal is the *site* and persists across the
whole run; the excitation is a *value at that site* in one scenario.

---

## 3. Sections

### 3.1 `simulation` (required)

Run-level settings. The only required section.

| Key | Type | Default | Meaning |
|-----|------|---------|---------|
| `physics_type` | string | **required** | `electrostatics`, `magnetostatics`, `magnetoquasistatics` |
| `mesh` | string | **required** | Mesh path; relative paths resolve against the config file's directory |
| `geometry_type` | string | `planar` | `planar` or `axisymmetric` |
| `analysis_type` | string | `field` | `field` or `coupling_matrix` |
| `order` | integer | `1` | FE polynomial order, 1-10 |
| `linear_solver` | string | `direct` | `direct` or `iterative` |
| `solver_tolerance` | number | see constants | Relative residual tolerance, in (0, 1) |
| `solver_max_iter` | integer | see constants | Iterative solver cap, >= 1 |
| `solver_print_level` | integer | see constants | Solver verbosity |
| `output_paraview` | bool | `false` | Write ParaView `.pvd`/`.vtu` |
| `output_gmsh` | bool | `false` | Write Gmsh `.msh` results |
| `results_path` | string | mesh directory | Output directory |
| `export_refine` | integer | follows `order` | Export subdivision factor, >= 1 |
| `amr` | object | absent = disabled | See below |

`direct` is the default solver because a coupling-matrix run amortizes one
factorization over every terminal's right-hand side, and its accuracy does not
depend on a residual tolerance.

#### `simulation.amr`

| Key | Type | Default | Meaning |
|-----|------|---------|---------|
| `enabled` | bool | `false` | Master switch |
| `max_iterations` | integer | `5` | Refine/re-solve passes |
| `max_dofs` | integer | `2000000` | Stop past this DOF count (<= 0 disables) |
| `error_fraction` | number | `0.7` | Bulk (Dorfler) marking fraction, in (0, 1] |
| `error_tolerance` | number | `0.0` | Absolute stop threshold (<= 0 ignores) |
| `conforming` | bool | `true` | Require conforming output |

When AMR is absent or disabled the solver performs a single solve.

### 3.2 `entity_groups`

Array of named attribute-id sets. The only section that mentions raw mesh
numbers.

| Key | Type | Required | Meaning |
|-----|------|----------|---------|
| `name` | string | yes | Unique group name |
| `dim` | integer | yes | Entity dimension: `1` = curve, `2` = surface, `3` = volume |
| `attribute_ids` | array of int | yes | Positive mesh attribute ids |

```json
"entity_groups": [
  { "name": "Air",         "dim": 2, "attribute_ids": [1] },
  { "name": "Dielectric",  "dim": 2, "attribute_ids": [2] },
  { "name": "TopPlate",    "dim": 1, "attribute_ids": [1] },
  { "name": "FarField",    "dim": 1, "attribute_ids": [3] }
]
```

`dim` is the topological dimension of the entities the ids name. It is stated,
not inferred, because it selects *which numbering space the ids live in*. Gmsh
numbers physical groups independently per dimension, so one model may use id `1`
for both a curve and a surface; MFEM preserves that split as `bdr_attributes`
versus `attributes`. In the example above `Air` and `TopPlate` both use id `1`
and refer to entirely different entities. Without `dim` the ids are ambiguous.

Boundary versus domain is **derived** from `dim`, not authored:

| Relation to mesh dimension | Role |
|---|---|
| `dim == mesh dimension` | domain -- carries a material and the PDE |
| `dim == mesh dimension - 1` | boundary -- bounds the solved region |

So in a 2D mesh `dim: 2` is a domain and `dim: 1` is a boundary; in 3D those
become `3` and `2`. Any other value names entities the solver cannot use and is
rejected, naming both acceptable values. The derived role is then checked
against use: a boundary condition must reference a boundary group, a region must
reference a domain group.

### 3.3 `materials`

Array of named property bundles. A material has no location.

| Key | Type | Required | Meaning |
|-----|------|----------|---------|
| `name` | string | yes | Unique material name |
| `properties.sigma` | number | no (0.0) | Conductivity [S/m] |
| `properties.epsilon_r` | number | no (1.0) | Relative permittivity |
| `properties.mu_r` | number | no (1.0) | Relative permeability |

Defaults are vacuum. Only the properties the physics needs are read:
electrostatics uses `epsilon_r`, magnetostatics `mu_r`, MQS `mu_r` and `sigma`.

### 3.4 `regions`

Array binding materials to domain entity groups.

| Key | Type | Required | Meaning |
|-----|------|----------|---------|
| `entity_group` | string | yes | A `domain` group |
| `material` | string | yes | A defined material name |
| `current_constraint` | string | no (`none`) | `none` or `open` |

**Every domain attribute in the mesh must be claimed by exactly one region.**
Unclaimed attributes and doubly-claimed attributes are both errors. This is
strict on purpose: an unclaimed domain silently receives default (vacuum)
properties and produces a plausible but wrong answer.

`current_constraint: open` marks a conductor carrying no net current -- an
induced-current-only body.

### 3.5 `boundary_conditions`

Array of boundary conditions.

| Key | Type | Required | Meaning |
|-----|------|----------|---------|
| `entity_group` | string | yes | A `boundary` group |
| `type` | string | yes | `dirichlet`, `neumann`, or `robin` |
| `value` | number | yes | Prescribed value or flux |
| `name` | string | no | Documentation only |
| `robin_coefficient` | number | Robin only | Reserved |

- **dirichlet** prescribes the solution: potential `V` (electrostatics) or
  `A_phi` (magnetics).
- **neumann** prescribes the outward natural flux. Value `0` is the implicit
  natural condition, so a boundary with no entry behaves as homogeneous Neumann.
- **robin** is parsed and reserved but **not implemented**; the solver rejects it.
  `robin_coefficient` is required for Robin entries and rejected on all others.

Entries omitted entirely are homogeneous Neumann. Axis regularity on `r = 0` in
axisymmetric magnetic runs is imposed automatically and must **not** be authored.

### 3.6 `terminals`

Array of named drive/measurement sites.

| Key | Type | Required | Meaning |
|-----|------|----------|---------|
| `name` | string | yes | Unique terminal name |
| `quantity` | string | yes | `voltage` or `current` |
| `entity_group` | string | yes | Group; role depends on quantity |
| `conductor_type` | string | no (`massive`) | `massive` or `stranded` |

`quantity` decides which role of group is required, and how the terminal is
realized:

| `quantity` | Group role | Realization |
|------------|------------|-------------|
| `voltage` | boundary | Essential constraint on the boundary |
| `current` | domain | Source term over the conductor domain |

`quantity` has no default. A wrong guess here would silently change the physics
while still solving, so it must be stated.

`conductor_type` applies to current terminals in MQS: `stranded` imposes uniform
current density (litz/fine-wire, eddy currents suppressed); `massive` solves for
the true current distribution including skin and proximity effects.

### 3.7 `scenarios`

Array of solves. Ignored when `analysis_type` is `coupling_matrix`, which
synthesizes its own unit-drive scenarios.

| Key | Type | Required | Meaning |
|-----|------|----------|---------|
| `name` | string | yes | Scenario name; labels output |
| `excitations` | array | no | Per-terminal values |
| `frequency` | number/array/object | MQS only | See below |

Each excitation:

| Key | Type | Meaning |
|-----|------|---------|
| `terminal` | string | A defined terminal name |
| `value` | number | Volts or amps, per that terminal's type |

**A terminal omitted from `excitations` defaults to zero** of its quantity:
grounded for voltage, open for current. Omission is meaningful, not an error.

#### Frequency (MQS only)

Three forms, all expanding to one solve per point:

```json
"frequency": 60.0                                              // single
"frequency": [10.0, 100.0, 1000.0]                             // explicit list
"frequency": { "scale": "log", "start": 10.0,
			   "stop": 1000.0, "points": 5 }                   // sweep
```

`scale` is `log` or `linear`. `points` includes both endpoints; `points: 1`
solves only at `start`. Swept scenarios are named
`<name>_f<n>_<frequency>Hz`. Frequency is required and positive for MQS and
ignored by the static solvers.

---

## 4. Analysis types

### `field`

Solves the authored scenarios and writes the resulting fields.

### `coupling_matrix`

Ignores authored scenarios. Drives each terminal to unity in turn, all others at
zero, and assembles the terminal-by-terminal matrix:

| Physics | Matrix | Unit |
|---------|--------|------|
| Electrostatics | Capacitance | F |
| Magnetostatics | Inductance | H |
| MQS | Resistance + inductance, per frequency point | Ohm, H |

Rows and columns follow the terminal definition order. For MQS, each scenario
supplies a frequency point and a separate labeled R/L pair is written per point.

---

## 5. Validation rules

Checked before solving, and reported together rather than one at a time.

**Structural**
- Root is an object; `simulation` is present.
- Every section has the expected JSON type.
- Superseded names are rejected with the replacement named:
  `physics` -> `physics_type`, `geometry` -> `geometry_type`,
  `results_file` -> `results_path`, `boundaries` -> `boundary_conditions`,
  entity group `kind` -> `dim`, terminal `excitation_type`/`excitation` ->
  `quantity`. A capitalized boundary condition type is reported with its
  lowercase spelling rather than a generic rejection.

**Names and references**
- Entity group names are unique.
- Every `entity_group` reference resolves to a defined group.
- Every group's `dim` is the mesh dimension (domain) or one below it
  (boundary); anything else is rejected, naming both acceptable values.
- Every referenced group has the role its use requires.
- Every `material` reference resolves to a defined material.
- Every `terminal` reference in an excitation resolves to a defined terminal.

**Coverage and conflicts**
- Every mesh domain attribute is claimed by exactly one region.
- Attribute ids are positive integers.
- No DOF is pinned to two different values by two boundary conditions.
- No boundary entity carries two boundary conditions.
- No DOF belongs to two terminals, or to a terminal and a Dirichlet boundary
  condition at once. Because terminal values vary per scenario, any shared DOF
  is guaranteed to conflict in some scenario.

**Ranges**
- `order` in [1, 10]; `solver_tolerance` in (0, 1); `solver_max_iter` >= 1;
  `export_refine` >= 1; `amr.error_fraction` in (0, 1];
  `amr.max_iterations` and `amr.error_tolerance` non-negative.

**Physics-specific**
- MQS requires a positive frequency on every scenario.
- Robin boundary conditions are rejected as unimplemented.

---

## 6. Worked example

Axisymmetric parallel-plate capacitor with a dielectric, in air, truncated by a
grounded far-field boundary.

```json
{
  "simulation": {
	"physics_type": "electrostatics",
	"geometry_type": "axisymmetric",
	"mesh": "capacitor.mesh",
	"order": 2,
	"analysis_type": "coupling_matrix",
	"output_paraview": true
  },

  "entity_groups": [
	{ "name": "Air",         "dim": 2, "attribute_ids": [1] },
	{ "name": "Dielectric",  "dim": 2, "attribute_ids": [2] },
	{ "name": "TopPlate",    "dim": 1, "attribute_ids": [1] },
	{ "name": "BottomPlate", "dim": 1, "attribute_ids": [2] },
	{ "name": "FarField",    "dim": 1, "attribute_ids": [3] }
  ],

  "materials": [
	{ "name": "Air",        "properties": { "epsilon_r": 1.0 } },
	{ "name": "Dielectric", "properties": { "epsilon_r": 2.5 } }
  ],

  "regions": [
	{ "entity_group": "Air",        "material": "Air" },
	{ "entity_group": "Dielectric", "material": "Dielectric" }
  ],

  "boundary_conditions": [
	{ "entity_group": "FarField", "type": "dirichlet", "value": 0.0 }
  ],

  "terminals": [
	{ "name": "Top",    "quantity": "voltage", "entity_group": "TopPlate" },
	{ "name": "Bottom", "quantity": "voltage", "entity_group": "BottomPlate" }
  ]
}
```

Note what is *absent*: no `scenarios`, because `coupling_matrix` synthesizes
them. Switching to `"analysis_type": "field"` would require adding:

```json
"scenarios": [
  {
	"name": "Charged",
	"excitations": [ { "terminal": "Top", "value": 1.0 } ]
  }
]
```

`Bottom` is omitted and therefore grounded.

---

## 7. Terminology review

Terms were assessed against one test: **does the name say what the thing is?**

### Kept

| Term | Why |
|------|-----|
| **entity group** | Role-neutral, so one concept covers both boundary and domain naming. `physical_group` (Gmsh) would tie the schema to one mesher; `attribute_set` (MFEM) leaks the FE library into the user's file. |
| **region** | Standard for a material-assigned subdomain, and correctly distinct from the material itself. |
| **material** | Universal. |
| **terminal** | Matches COMSOL. Deliberately not `port`, which means a modal/S-parameter boundary and would mislead in a static/quasi-static code. |
| **excitation** | Matches Ansys, where prescribed drives live under Excitations rather than boundary conditions. Correctly implies a scenario-supplied value. |
| **scenario** | Neutral and accurate: one parameterized solve. `case` and `step` both carry unwanted baggage (`step` implies time-stepping). |
| **coupling matrix** | Names the content -- terminal-to-terminal coupling -- rather than merely the index set. |

### Corrected

Four names disagreed with the model they described. Each was renamed rather than
kept, since the internal model already used the better name in every case and the
JSON surface was the only thing out of step. All four superseded spellings are
now rejected with a message naming the replacement.

| Was | Now | Why |
|-----|-----|-----|
| `boundaries` | `boundary_conditions` | The section defines *conditions*, not geometry -- and the geometry it references is already named by entity groups. The internal type was already `BoundaryCondition`. |
| `excitation_type` | `quantity` | The value is `voltage` or `current` -- a quantity, not a type of excitation. Matches the internal `Quantity` enum. |
| `"Dirichlet"` | `"dirichlet"` | Every other enum value in the schema is lowercase. A capitalized value now reports its lowercase spelling rather than a bare rejection. |

One internal name moved with them: `Terminal::ExcitationType` ->
`Terminal::DriveQuantity` (`Quantity` alone would collide with the enum type).

`entity_groups[].dim` was briefly renamed to `kind: "boundary" | "domain"` and
then reverted. The role is derived from the dimension by comparing it to the
mesh dimension, so authoring it separately stored a fact the mesh already
determines and allowed the two to disagree. More importantly, `dim` names the
mesher's own concept -- Gmsh and MFEM both key their attribute namespaces on
entity dimension -- and so keeps its meaning across tools and in 3D, where a
boundary is a surface rather than a curve. `EntityKind` was likewise reverted to
a plain `int Dim` with derived `IsDomain()`/`IsBoundary()` helpers.

### Deliberate, not an inconsistency

**`terminals[].name` vs. `excitations[].terminal`.** Definition sites use `name`;
reference sites use the referenced type as the key. This holds throughout the
schema and is noted only to confirm it is intentional.

### Rejected

| Proposed | Rejected because |
|----------|------------------|
| `TerminalMatrix` for coupling matrix | Names the index set, not the content |
| Separate exterior/interior boundary kinds | Nothing branches on it; entity group names already convey it to the reader |
| `port` for terminal | Reserved for modal/S-parameter use |
| Folding terminals into boundary conditions | Only voltage terminals are boundary conditions; current terminals are domain sources |
