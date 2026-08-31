# Example: Eddy Current in Conducting Cylinder

This example demonstrates time-harmonic eddy current analysis for a conducting cylinder in an alternating magnetic field.

## Problem Description

**Geometry:**
- Conducting cylinder: radius R = 0.05 m, length L = 0.1 m
- Surrounding coil: creates time-varying magnetic field
- Axisymmetric configuration

**Physics:**
- Frequency sweep: f = 10 Hz to 1 kHz (five logarithmically spaced points)
- The 60 Hz calculations below provide a power-frequency reference
- Angular frequency: ω = 2πf = 377 rad/s
- Conductor: Aluminum (σ = 3.5 × 10⁷ S/m, μᵣ = 1.0)
- Coil: stranded winding driven with a total current of 240 A, which over the
  0.0024 m² coil cross-section gives J_source = 1 × 10⁵ A/m²

**Boundary Conditions:**
- Far field: A_φ = 0 (Dirichlet)

**Skin Depth:**
```
δ = √(2 / (ω μ σ))
  = √(2 / (377 × 4π×10⁻⁷ × 3.5×10⁷))
  = 1.47 mm
```

The field penetrates only ~1.5 mm into the aluminum at 60 Hz.

## Physical Phenomena

### 1. Skin Effect
- Induced currents flow primarily near surface
- Exponential decay with depth: `e^(-x/δ)`
- Current density maximum at surface

### 2. Eddy Current Losses
- Power dissipation: `P = ∫ (σ/2) |E|² dV`
- Joule heating in conductor
- Phase lag between applied field and induced current

### 3. Shielding Effect
- Eddy currents generate opposing magnetic field
- Reduced field penetration into conductor
- Effective permeability appears < μ_r at high frequency

## Analytical Approximation

For a thin conducting sheet:

**Surface resistance:**
```
R_s = √(ω μ / (2 σ))
    = √(377 × 4π×10⁻⁷ / (2 × 3.5×10⁷))
    = 5.2 × 10⁻⁵ Ω/square
```

**Power loss per unit area:**
```
P/A = ½ R_s |H_surface|²
```

## Running the Example

```bash
# From project root
mkdir -p build && cd build
cmake ..
make

# Run simulation
./mfem-electromag ../examples/eddy_current/config.json

# Results will be in results_mqs/
```

## Expected Results

The simulation produces complex-valued fields:

1. **A_real, A_imag:** In-phase and quadrature components of vector potential
2. **B_real, B_imag:** Real and imaginary parts of magnetic flux density
3. **B_magnitude:** RMS magnitude = √(B_real² + B_imag²)
4. **Power loss:** Concentrated near conductor surface

## Visualization

```bash
# Open in ParaView
paraview results_mqs/results_mqs.pvd
```

**Key visualizations:**
- B_Magnitude: Shows field concentration and shielding
- A_Real vs A_Imag: Phase relationship
- Line plot through conductor: Demonstrates skin effect
- Time animation: Use Temporal Interpolator filter

### Time-Domain Animation

To animate the time-harmonic solution:

1. In ParaView, add "Calculator" filter:
   ```
   B_time = B_Real*cos(2*pi*60*t) - B_Imag*sin(2*pi*60*t)
   ```

2. Add "Programmable Source" to generate time steps

3. Animate to see oscillating field

## Mesh Considerations

**Critical:** Mesh must resolve skin depth!

- **In conductor:** Element size ≤ δ/3 ≈ 0.5 mm near surface
- **Boundary layer:** Use graded mesh from surface inward
- **Air region:** Coarser mesh acceptable (5-10 mm)

**Typical mesh:**
- Surface layer: 5-10 elements within 3δ ≈ 4.5 mm
- Geometric growth ratio: 1.2-1.5
- Total elements: 10,000-50,000 depending on order

## Validation

### 1. Skin Depth Check
Plot |B| vs depth into conductor:

```
|B(x)| / |B(0)| ≈ e^(-x/δ)
```

At depth x = δ, field should drop to ~37% (1/e) of surface value.

### 2. Power Loss
Compare computed losses to analytical for simple geometry.

### 3. Phase Relationship
Inside conductor:
- Current lags applied field by ~45°
- Phase increases with depth

## Frequency Sweep

Each MQS scenario requires either a positive scalar frequency or an inclusive
linear/logarithmic range. The example uses:

```json
"frequency": {
  "scale": "log",
  "start": 10.0,
  "stop": 1000.0,
  "points": 5
}
```

Use `"scale": "linear"` for uniform spacing. Both endpoints are included;
`"points": 1` solves only `start`. The excitation list is copied to every
expanded frequency point, and output scenario names include the point and
frequency.

For a coupling-matrix sweep, set `"analysis_type": "coupling_matrix"` and
provide one or more frequency scenarios. Terminal excitations in those
scenarios are ignored because the solver synthesizes each unit-current column;
each frequency receives its own labeled resistance and inductance CSV files.

**Expected trends:**
- Higher f → smaller δ (stronger skin effect)
- Higher f → greater power loss
- Higher f → better shielding

**At 1 kHz:**
```
δ = 1.47 mm / √(1000/60) = 0.36 mm
```

Much more confined to surface!

## Material Variations

Compare different conductors:

| Material | σ (S/m) | δ @ 60 Hz | Application |
|----------|---------|-----------|-------------|
| Aluminum | 3.5×10⁷ | 1.5 mm | Lightweight |
| Copper   | 5.8×10⁷ | 1.1 mm | High conductivity |
| Steel    | 1.0×10⁶ | 7.0 mm | Structural, magnetic |
| Carbon   | 1.0×10⁴ | 230 mm | Composite materials |

## Applications

This example models:
- **Induction heating:** Controlled power deposition in conductors
- **Eddy current testing (NDT):** Defect detection via field distortion
- **Transformer cores:** Laminations reduce eddy current losses
- **Electromagnetic braking:** Drag force from induced currents
- **Shielding:** Protection from EMI

## Advanced Analysis

### 1. Loss Calculation

Integrate power loss density:

```
P_total = ∫_conductor (σ ω² / 2) |A|² dV
```

Extract from simulation output.

### 2. Force Calculation

Lorentz force on conductor:

```
F = ∫ J × B dV
```

Causes repulsion or attraction depending on phase.

### 3. Impedance

For coil + conductor system:

```
Z = R + jωL
```

where R includes eddy current losses, L is effective inductance.
