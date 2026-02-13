# MFEM-ElectroMag Examples

This directory contains example problems demonstrating electromagnetic field simulations using MFEM-ElectroMag.

## Examples

1. **simple_capacitor/** - Parallel plate capacitor (electrostatics)
2. **solenoid/** - Axisymmetric solenoid (magnetostatics)
3. **eddy_current/** - Conducting cylinder in AC field (magnetoquasistatics)

Each example includes:
- `config.json` - Simulation configuration file
- `README.md` - Detailed problem description and physics
- `*.geo` - Gmsh geometry script for mesh generation

## Quick Start

### 1. Generate Meshes

First, generate the meshes for all examples:

```bash
cd examples
./generate_meshes.sh
```

This will create:
- `eddy_current/eddy_current.mesh`
- `simple_capacitor/capacitor.mesh`
- `solenoid/solenoid.mesh`

**Requirements:**
- Gmsh must be installed: `sudo apt-get install gmsh` (Ubuntu/Debian) or `brew install gmsh` (macOS)
- Alternatively, download from https://gmsh.info/

### 2. Build the Solver

```bash
cd ..
mkdir -p build && cd build
cmake ..
make
```

### 3. Run Examples

```bash
# From build directory
./mfem-electromag ../examples/simple_capacitor/config.json
./mfem-electromag ../examples/solenoid/config.json
./mfem-electromag ../examples/eddy_current/config.json
```

Results will be written to `results_*/` directories.

### 4. Visualize Results

```bash
paraview results_electrostatic/results_electrostatic.pvd
paraview results_magnetostatic/results_magnetostatic.pvd
paraview results_mqs/results_mqs.pvd
```

## Mesh Generation Details

### Manual Mesh Generation

To generate a single mesh:

```bash
cd examples/solenoid
gmsh -2 -format msh2 solenoid.geo -o solenoid.msh
```

The `-format msh2` flag outputs ASCII MSH2 format, which MFEM can read directly.

### Mesh Refinement

Each `.geo` file includes mesh size fields for appropriate refinement:

- **Eddy current**: Fine mesh (0.5 mm) near conductor surface to resolve skin depth
- **Capacitor**: Refined mesh near plate edges to capture fringing fields
- **Solenoid**: Uniform mesh in coil region, graded toward far field

Modify the `lc_*` parameters in the `.geo` files to adjust mesh density.

### Viewing Meshes in Gmsh

To visualize the mesh before running simulations:

```bash
gmsh eddy_current/eddy_current.geo
# In Gmsh: Tools → Statistics to see element count
# Press 'e' to show element edges
# Press '2' to show surface elements
```

## Modifying Examples

### Change Geometry

Edit the `.geo` file parameters:
```c
// In eddy_current.geo
cyl_radius = 0.05;      // Change cylinder size
coil_r_inner = 0.06;    // Change coil position
```

Then regenerate the mesh:
```bash
gmsh -2 -format msh2 eddy_current.geo -o eddy_current.msh
```

### Change Materials/Physics

Edit the `config.json` file:
```json
{
  "materials": [
    {
      "name": "Conductor",
      "properties": {
        "sigma": 5.8e7  // Change to copper conductivity
      }
    }
  ]
}
```

### Change Mesh Resolution

For higher accuracy or faster testing:

**Quick test mesh** (coarse):
```bash
# In .geo file, increase all lc_* values by factor of 2-3
lc_conductor_surface = 0.002;  // Was 0.0005
```

**High accuracy** (fine):
```bash
# Decrease all lc_* values by factor of 2
lc_conductor_surface = 0.00025;  // Was 0.0005
```

## Troubleshooting

### "gmsh: command not found"

Install Gmsh:
- **Ubuntu/Debian**: `sudo apt-get install gmsh`
- **macOS**: `brew install gmsh`
- **Windows**: Download from https://gmsh.info/

### Mesh generation is slow

- Reduce mesh refinement in `.geo` files (increase `lc_*` values)
- The eddy current mesh is intentionally fine to resolve skin depth

### "Unknown element type" error

Make sure Gmsh outputs MSH2 format: `gmsh -2 -format msh2 ...`

MFEM may not support newer MSH4 format elements.

### Simulation results look wrong

1. Check mesh quality: Open in Gmsh and look for badly shaped elements
2. Verify material attributes match between mesh and config.json
3. Check boundary conditions are applied correctly
4. Increase mesh resolution if solution is under-resolved

## Example Physics Summary

| Example | Type | Key Physics | Typical Runtime |
|---------|------|-------------|-----------------|
| simple_capacitor | Electrostatics | Fringing fields, capacitance | < 1 min |
| solenoid | Magnetostatics | Magnetic field, inductance | < 1 min |
| eddy_current | Magnetoquasistatics | Skin effect, eddy losses | 2-5 min |

All examples use axisymmetric formulation for efficiency.

## Further Information

See individual example README files for:
- Detailed problem descriptions
- Analytical solutions for validation
- Expected results and visualization tips
- Parameter studies and variations
