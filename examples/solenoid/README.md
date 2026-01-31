# Example: Axisymmetric Solenoid

This example calculates the magnetostatic field of a current-carrying solenoid coil.

## Problem Description

**Geometry:**
- Solenoid: inner radius r_in = 0.05 m, outer radius r_out = 0.08 m
- Length: L = 0.2 m
- Number of turns: N = 1000
- Current per turn: I = 1 A
- Axisymmetric about central axis

**Boundary Conditions:**
- Far field: A_φ = 0 (Dirichlet)
- Axis (r = 0): Natural boundary (∂A/∂r = 0)

**Materials:**
- Coil region: Current density J_φ = N × I / A_coil
- Air core: μᵣ = 1.0
- Surrounding air: μᵣ = 1.0

**Optional:** Add iron core (μᵣ = 1000) to demonstrate field concentration

## Analytical Solution

For a long solenoid, the magnetic field inside is approximately uniform:

**Magnetic field inside solenoid:**
```
B = μ₀ N I / L = 4π × 10⁻⁷ × 1000 × 1 / 0.2
  = 6.28 × 10⁻³ T
  = 6.28 mT
```

**Current density in coil:**
```
J = N I / A_coil
```

where `A_coil = (r_out - r_in) × L = 0.03 × 0.2 = 0.006 m²`

```
J = 1000 × 1 / 0.006 = 1.67 × 10⁵ A/m²
```

## Running the Example

```bash
# From project root
mkdir -p build && cd build
cmake ..
make

# Run simulation
./mfem-electromag ../examples/solenoid/config.json

# Results will be in results_magnetostatic/
```

## Expected Results

The simulation should produce:
1. **Vector potential A_φ:** Increases inside coil, drops to zero at far field
2. **Magnetic flux density B:** ~6 mT axial field inside solenoid
3. **Field lines:** Closed loops through coil and return path in air
4. **Fringing:** Field spreads near ends of solenoid

## Visualization

```bash
# Open in ParaView
paraview results_magnetostatic/results_magnetostatic.pvd
```

**Suggested visualizations:**
- Vector field of B (magnetic flux density)
- Streamlines showing field lines
- Contour plot of |B| magnitude
- Line plot along axis to show field uniformity

## Mesh Generation

The mesh should include:
- Coil region (source of current)
- Air core inside coil
- Air region outside coil
- Far-field boundary (large enough to avoid boundary effects)

Recommended mesh:
- Element size in coil: ~5 mm
- Element size in air: graded from 5 mm to 50 mm
- Far-field radius: at least 5 × coil radius

## Validation

Compare FEM results to analytical:

**Inside solenoid (center):**
```
|B_z,FEM - B_analytical| / B_analytical < 2%
```

**End effects:** Near ends, expect ~50% of central field value.

## With Iron Core

Add iron core (μᵣ = 1000) inside solenoid:

```json
{
  "name": "IronCore",
  "attributes": [3],
  "properties": {
    "mu_r": 1000.0
  }
}
```

**Expected changes:**
- Field inside core increases by factor of ~μᵣ
- Field concentration along axis
- Reduced fringing outside coil

## Variations

Try:
1. **Current:** Vary I (0.1 A to 10 A) - field scales linearly
2. **Turns:** Change N - field scales with N
3. **Core material:** Compare air (μᵣ=1) vs iron (μᵣ=1000)
4. **Geometry:** Vary length-to-diameter ratio
5. **Mesh refinement:** Increase order or mesh density

## Advanced: Inductance Calculation

Calculate inductance from stored magnetic energy:

```
L = 2 U / I²
```

where `U = ∫ (B²/2μ) dV` is the magnetic energy.

Extract from simulation:
```
U = ∫ (B·H/2) dV
L ≈ N Φ / I
```

For this geometry, expect L ≈ 10-100 mH depending on core.
