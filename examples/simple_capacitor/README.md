# Example: Parallel Plate Capacitor

This example demonstrates electrostatic field calculation for a simple parallel plate capacitor in axisymmetric geometry.

## Problem Description

**Geometry:**
- Two circular plates of radius R = 0.1 m
- Separation distance d = 0.01 m (1 cm)
- Axisymmetric about central axis

**Boundary Conditions:**
- Top plate: V = 1000 V (Dirichlet)
- Bottom plate: V = 0 V (Dirichlet, ground)
- Far field: V = 0 V

**Material:**
- Dielectric between plates: εᵣ = 2.5 (e.g., paper, FR4)
- Surrounding air: εᵣ = 1.0

## Analytical Solution

For an ideal parallel plate capacitor:

**Capacitance:**
```
C = ε₀ εᵣ A / d = ε₀ εᵣ π R² / d
```

where:
- ε₀ = 8.854 × 10⁻¹² F/m
- εᵣ = 2.5
- R = 0.1 m
- d = 0.01 m

```
C = 8.854e-12 * 2.5 * π * (0.1)² / 0.01
  = 6.95 × 10⁻¹¹ F
  ≈ 69.5 pF
```

**Electric Field (between plates):**
```
E = V / d = 1000 V / 0.01 m = 100,000 V/m
```

**Energy Stored:**
```
U = ½ C V² = ½ × 69.5e-12 × 1000² = 34.75 μJ
```

## Running the Example

```bash
# From project root
mkdir -p build && cd build
cmake ..
make

# Run simulation
./mfem-electromag ../examples/simple_capacitor/config.json

# Results will be in results_electrostatic/
```

## Expected Results

The simulation should produce:
1. **Electric potential V:** Linear variation from 1000V to 0V between plates
2. **Electric field E:** Nearly uniform field of ~100 kV/m between plates
3. **Fringing fields:** Edge effects near plate boundaries

## Visualization

```bash
# Open in ParaView
paraview results_electrostatic/results_electrostatic.pvd
```

**Suggested visualizations:**
- Contour plot of electric potential V
- Vector field of electric field E
- Color map of field magnitude |E|

## Mesh Generation

The mesh should capture:
- Plate geometry (thin conductors)
- Dielectric region between plates
- Air region around plates
- Far-field boundary

You can create a mesh using:
- **Gmsh:** See `capacitor.geo` for geometry definition
- **MFEM mesh converter:** Convert from other formats

## Validation

Compare FEM results to analytical:
```
|E_FEM - E_analytical| / E_analytical < 1%  (between plates, away from edges)
```

Edge effects will cause local deviations near r = R.

## Variations

Try modifying:
1. **Dielectric constant:** Change εᵣ (1.0 to 10.0)
2. **Voltage:** Vary applied voltage (100V to 10kV)
3. **Geometry:** Change R or d to see capacitance scaling
4. **Mesh refinement:** Increase element order or mesh density
