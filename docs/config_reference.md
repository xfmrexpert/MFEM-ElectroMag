# Configuration File Reference

Exhaustive listing of every key the solver reads, its type, whether it is
required, and its default. This document is the normative reference; for the
model behind these keys and the reasoning that shaped them, see
[input_file_format.md](input_file_format.md).

Defaults shown are the literal values applied in code. A key listed as
**required** is enforced by `ConfigValidator`; a key with a default may be
omitted entirely.

---

## Contents

- [`simulation`](#simulation)
- [`simulation.amr`](#simulationamr)
- [`entity_groups`](#entity_groups)
- [`materials`](#materials)
- [`regions`](#regions)
- [`boundary_conditions`](#boundary_conditions)
- [`terminals`](#terminals)
- [`scenarios`](#scenarios)
- [Validation rules](#validation-rules)
- [Superseded names](#superseded-names)

---

## `simulation`

Object. **Required** -- the only required section.

| Key | Type | Default | Range / values |
|-----|------|---------|----------------|
| `physics_type` | string | **required** | `electrostatics`, `magnetostatics`, `magnetoquasistatics` |
| `mesh` | string | **required** | Path; relative resolves against the config file's directory |
| `geometry_type` | string | `planar` | `planar`, `axisymmetric` |
| `analysis_type` | string | `field` | `field`, `coupling_matrix` |
| `order` | integer | `1` | 1-10 |
| `linear_solver` | string | `direct` | `direct`, `iterative` |
| `solver_tolerance` | number | `1e-12` | (0, 1) |
| `solver_max_iter` | integer | `1000` | >= 1 |
| `solver_print_level` | integer | `1` | -- |
| `output_paraview` | bool | `false` | -- |
| `output_gmsh` | bool | `false` | -- |
| `results_path` | string | mesh directory | Relative resolves against the config file's directory |
| `amr` | object | absent = disabled | See below |

`direct` is the default linear solver because a coupling-matrix run amortizes
one factorization over every terminal's right-hand side, and its accuracy does
not depend on a residual tolerance.

`frequency` is **not** valid here. It belongs on each scenario; see
[`scenarios`](#scenarios).

## `simulation.amr`

Object. Optional; absent or `enabled: false` means a single solve.

| Key | Type | Default | Range / values |
|-----|------|---------|----------------|
| `enabled` | bool | `false` | Master switch |
| `max_iterations` | integer | `5` | >= 0; refine/re-solve passes |
| `max_dofs` | integer | `2000000` | <= 0 disables the cap |
| `error_fraction` | number | `0.7` | (0, 1]; bulk (Dorfler) marking fraction |
| `error_tolerance` | number | `0.0` | >= 0; <= 0 ignores the threshold |
| `conforming` | bool | `true` | Require conforming output |

Unknown keys in this block are ignored, so a producer and the solver can evolve
independently.

---

## `entity_groups`

Array of objects. The only section that names raw mesh attribute ids.

| Key | Type | Required | Meaning |
|-----|------|----------|---------|
| `name` | string | yes | Unique group name |
| `dim` | integer | yes | `1` = curve, `2` = surface, `3` = volume |
| `attribute_ids` | array of integer | yes | Positive mesh attribute ids |

Boundary versus domain is **derived** from `dim`, never prescribed:

| Relation to mesh dimension | Role |
|---|---|
| `dim == mesh dimension` | domain -- carries a material and the PDE |
| `dim == mesh dimension - 1` | boundary -- bounds the solved region |

Any other value is rejected. `dim` is stated rather than inferred because it
selects which numbering space the ids live in: Gmsh numbers physical groups
independently per dimension, so id `1` may name both a curve and a surface.

---

## `materials`

Array of objects. A material is a property bundle with no location.

| Key | Type | Required | Default | Meaning |
|-----|------|----------|---------|---------|
| `name` | string | yes | -- | Unique material name |
| `properties.sigma` | number | no | `0.0` | Conductivity [S/m] |
| `properties.epsilon_r` | number | no | `1.0` | Relative permittivity |
| `properties.mu_r` | number | no | `1.0` | Relative permeability |

Defaults are vacuum. Only the properties the physics needs are read:
electrostatics uses `epsilon_r`, magnetostatics `mu_r`, MQS `mu_r` and `sigma`.

---

## `regions`

Array of objects binding materials to domain entity groups.

| Key | Type | Required | Default | Meaning |
|-----|------|----------|---------|---------|
| `entity_group` | string | yes | -- | Must be a domain group |
| `material` | string | yes | -- | A defined material name |
| `current_constraint` | string | no | `none` | `none`, `open` |

Every domain attribute in the mesh must be claimed by exactly one region.
`current_constraint: open` marks a conductor carrying no net current -- an
induced-current-only body.

---

## `boundary_conditions`

Array of objects.

| Key | Type | Required | Meaning |
|-----|------|----------|---------|
| `entity_group` | string | yes | Must be a boundary group |
| `type` | string | yes | `dirichlet`, `neumann`, `robin` |
| `value` | number | yes | Prescribed value or flux |
| `name` | string | no | Documentation only |
| `robin_coefficient` | number | Robin only | Reserved; rejected on non-Robin entries |

- **dirichlet** prescribes the solution: potential `V` (electrostatics) or
  `A_phi` (magnetics).
- **neumann** prescribes the outward natural flux. Value `0` is the implicit
  natural condition.
- **robin** is parsed and reserved but **not implemented**; the solver rejects it.

Boundaries with no entry are homogeneous Neumann. Axis regularity on `r = 0` in
axisymmetric magnetic runs is imposed automatically and must **not** be
prescribed here.

---

## `terminals`

Array of objects naming drive/measurement sites.

| Key | Type | Required | Default | Meaning |
|-----|------|----------|---------|---------|
| `name` | string | yes | -- | Unique terminal name |
| `quantity` | string | yes | none | `voltage`, `current` |
| `entity_group` | string | yes | -- | Role depends on `quantity` |
| `conductor_type` | string | no | `massive` | `massive`, `stranded` |

| `quantity` | Required group role | Realization |
|------------|---------------------|-------------|
| `voltage` | boundary | Essential constraint on the boundary |
| `current` | domain | Source term over the conductor domain |

`quantity` has no default: a wrong guess would silently change the physics while
still solving.

`conductor_type` applies to current terminals in MQS. `stranded` imposes uniform
current density (litz/fine-wire, eddy currents suppressed); `massive` solves for
the true current distribution including skin and proximity effects.

---

## `scenarios`

Array of objects, one per solve. Ignored when `analysis_type` is
`coupling_matrix`, which synthesizes its own unit-drive scenarios.

| Key | Type | Required | Meaning |
|-----|------|----------|---------|
| `name` | string | yes | Scenario name; labels output |
| `excitations` | array | no | Per-terminal values |
| `frequency` | number / array / object | MQS only | See below |

Each entry of `excitations`:

| Key | Type | Required | Meaning |
|-----|------|----------|---------|
| `terminal` | string | yes | A defined terminal name |
| `value` | number | yes | Volts or amps, per that terminal's `quantity` |

A terminal omitted from `excitations` defaults to zero of its quantity:
grounded for voltage, open for current. Omission is meaningful, not an error.

> **Excitation values are PEAK (amplitude) phasors in time-harmonic runs.**
> There is no rms/peak selector and no conversion. See
> [faq.md](faq.md#are-excitations-peak-or-rms).

### `frequency` (MQS only)

Three accepted forms, all expanding to one solve per point:

```json
"frequency": 60.0                                              // single
"frequency": [10.0, 100.0, 1000.0]                             // explicit list
"frequency": { "scale": "log", "start": 10.0,
			   "stop": 1000.0, "points": 5 }                   // sweep
```

Sweep object keys:

| Key | Type | Required | Values |
|-----|------|----------|--------|
| `scale` | string | yes | `log`, `linear` |
| `start` | number | yes | > 0 |
| `stop` | number | yes | > 0 |
| `points` | integer | yes | >= 1; includes both endpoints |

With `points: 1` only `start` is solved. Swept scenarios are named
`<name>_f<n>_<frequency>Hz`. Frequency is required and positive for MQS and
ignored by the static solvers.

---

## Validation rules

All checks run before solving, and every failure is reported together rather
than one at a time.

**Structural**
- Root is an object; `simulation` is present.
- Every section has the expected JSON type.
- Superseded names are rejected with the replacement named (see below).

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
  `amr.error_fraction` in (0, 1];
  `amr.max_iterations` and `amr.error_tolerance` non-negative.

**Physics-specific**
- MQS requires a positive frequency on every scenario.
- `simulation.frequency` is rejected for MQS; frequency belongs on scenarios.
- Robin boundary conditions are rejected as unimplemented.

---

## Superseded names

Rejected with a message naming the replacement:

| Superseded | Use instead |
|------------|-------------|
| `simulation.physics` | `simulation.physics_type` |
| `simulation.type` | `simulation.physics_type` |
| `simulation.geometry` | `simulation.geometry_type` |
| `simulation.results_file` | `simulation.results_path` |
| `boundaries` | `boundary_conditions` |
| entity group `kind` | `dim` |
| terminal `excitation_type` | `quantity` |
| terminal `excitation` | `quantity` |

A capitalized boundary condition type (for example `"Dirichlet"`) is reported
with its lowercase spelling rather than a generic rejection.
