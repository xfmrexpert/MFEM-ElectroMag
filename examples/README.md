# MFEM-ElectroMag Examples

This directory contains example problems demonstrating electromagnetic field simulations
using MFEM-ElectroMag. All examples use the axisymmetric (r-z) formulation with second
order (`"order": 2`) elements.

## Examples

| Directory | Physics | Mesh | Notes |
|-----------|---------|------|-------|
| `simple_capacitor/` | Electrostatics | `capacitor_from_sl.mesh` | Parallel plate capacitor, fringing fields, capacitance |
| `solenoid/` | Magnetostatics | `solenoid.mesh` | Axisymmetric solenoid, inductance |
| `current_loop/` | Magnetostatics | `loop.mesh` | Single current loop with an analytical hand calculation notebook |
| `eddy_current/` | Magnetoquasistatics | `eddy_current.mesh` | Conducting cylinder in an AC field, skin effect, losses |

Each directory contains a `config.json` and its `.geo` source geometry. The generated
`.mesh` files are committed, so the examples run without installing Gmsh.
`simple_capacitor/`, `solenoid/`, and `eddy_current/` also include a `README.md` with the
detailed problem description; `current_loop/` includes `hand_calc.ipynb` instead.

## Quick Start

### 1. Build the solver

```bash
cmake -S . -B build
cmake --build build --config Release
```

On Windows the repository also ships CMake presets (`x64-debug`, `x64-release`) that
Visual Studio picks up automatically:

```pwsh
cmake --preset x64-release
cmake --build --preset x64-release
```

MFEM is fetched and built by CMake, so the first configure takes a while.

### 2. Run an example

The config path is the only required argument. Relative paths inside the config
(mesh file, `results_path`) are resolved against the config file's directory, so the
examples can be run from anywhere:

```bash
./build/mfem-electromag examples/simple_capacitor/config.json
./build/mfem-electromag examples/solenoid/config.json
./build/mfem-electromag examples/current_loop/config.json
./build/mfem-electromag examples/eddy_current/config.json
```

Useful options:

| Option | Purpose |
|--------|---------|
| `--results-path <dir>` | Override `simulation.results_path` (relative paths resolve against the CWD) |
| `--verbosity <0\|1\|2>` | `0` status/timing, `1` solver output, `2` diagnostics |
| `--machine-readable` | Emit JSON Lines progress on stdout for tooling |
| `--version` | Print version/build information |
| `-h`, `--help` | Show usage |

### 3. Output

Output is off by default and enabled per format in the config's `simulation` block:

```json
{
  "simulation": {
	"output_paraview": true,
	"output_gmsh": true,
	"results_path": "results"
  }
}
```

- `output_paraview` writes a ParaView collection per scenario into
  `results_<physics>_<scenario>/`.
- `output_gmsh` writes `<scenario>.results.msh` (Gmsh MSH 2.2 ASCII). The mesh is
  emitted as native high-order Lagrange elements at the solution order, and the file
  includes an `$InterpolationScheme` block so a consumer can reconstruct fields at
  arbitrary points using the same basis MFEM used.
- If `results_path` is not set, output is written next to the mesh file.

Of the shipped examples only `simple_capacitor` enables output (`output_paraview`);
add the flags above to the other configs if you want files written.

Visualize with:

```bash
paraview examples/simple_capacitor/"results_electrostatics_Top Plate"/data.pvd
gmsh <scenario>.results.msh
```

## Regenerating Meshes

The committed `.mesh` files are sufficient to run the examples. Regenerate only if you
change a `.geo` file.

`generate_meshes.sh` (bash; requires Gmsh on `PATH`) regenerates the eddy current,
capacitor, and solenoid meshes:

```bash
cd examples
./generate_meshes.sh
```

It runs `gmsh -2 -format msh2` on each `.geo`, then uses MFEM's `convert-mesh` utility if
available, otherwise simply renames the `.msh` to `.mesh` (MFEM reads Gmsh MSH2 directly).

Note: the script writes `simple_capacitor/capacitor.mesh`, but that example's config uses
`capacitor_from_sl.mesh`. `current_loop/loop.mesh` is not covered by the script.

To regenerate a single mesh manually:

```bash
cd examples/solenoid
gmsh -2 -format msh2 solenoid.geo -o solenoid.mesh
```

Use MSH2 format; MFEM does not read all MSH4 element types.

### Mesh refinement

Each `.geo` sets mesh size fields for appropriate refinement:

- **Eddy current**: fine mesh near the conductor surface to resolve skin depth
- **Capacitor**: refined near plate edges to capture fringing fields
- **Solenoid**: uniform in the coil region, graded toward the far field

Adjust the `lc_*` parameters in the `.geo` files to change mesh density: increase them for
a fast coarse test, decrease them for higher accuracy.

## Modifying Examples

### Geometry

Edit the `.geo` parameters, then regenerate the mesh:

```c
// In eddy_current.geo
cyl_radius = 0.05;      // Change cylinder size
coil_r_inner = 0.06;    // Change coil position
```

### Materials / physics

Edit `config.json`:

```json
{
  "materials": [
	{
	  "name": "Conductor",
	  "properties": { "sigma": 5.8e7 }
	}
  ]
}
```

Material attribute names must match the physical groups in the mesh. The solver validates
the config up front and reports unknown or mismatched fields by JSON path.

## Troubleshooting

**`gmsh: command not found`** — install Gmsh (`sudo apt-get install gmsh`,
`brew install gmsh`, or download from https://gmsh.info/).

**"Unknown element type"** — regenerate with `-format msh2`.

**No output files** — set `output_paraview` and/or `output_gmsh` to `true` in the config;
both default to `false`.

**Results look wrong** — check mesh quality in Gmsh, verify material attributes match the
mesh physical groups, confirm the boundary conditions, and refine the mesh if the solution
looks under-resolved.

## Physics Summary

| Example | Type | Key physics | Typical runtime |
|---------|------|-------------|-----------------|
| simple_capacitor | Electrostatics | Fringing fields, capacitance | < 1 min |
| solenoid | Magnetostatics | Magnetic field, inductance | < 1 min |
| current_loop | Magnetostatics | Loop field vs. analytical result | < 1 min |
| eddy_current | Magnetoquasistatics | Skin effect, eddy losses | 2-5 min |

## Further Information

See the individual example README files for detailed problem descriptions, analytical
solutions for validation, expected results, and parameter studies.
