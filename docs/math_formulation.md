# Mathematical Formulation

This document describes the mathematical formulation for the three physics types supported by MFEM-ElectroMag: Electrostatics, Magnetostatics, and Magnetoquasistatics, all in the axisymmetric coordinate system.

## Coordinate System

We use cylindrical coordinates `(r, z, φ)` where:
- `r` is the radial distance from the axis of symmetry
- `z` is the axial coordinate
- `φ` is the azimuthal angle

For axisymmetric problems, all field quantities are independent of `φ` (`∂/∂φ = 0`).

## 1. Electrostatics

### Strong Form

The electrostatic problem solves for the electric potential `V(r,z)`:

```
-∇ · (ε ∇V) = ρ/ε₀
```

where:
- `V` is the electric potential [V]
- `ε = ε_r ε₀` is the permittivity [F/m]
- `ρ` is the space charge density [C/m³]
- `ε₀ = 8.854 × 10⁻¹² F/m` is the permittivity of free space

### Weak Form

Find `V ∈ H¹` such that for all test functions `v`:

```
∫_Ω ε ∇V · ∇v dΩ = ∫_Ω ρ v dΩ + ∫_∂Ω g v dS
```

### Axisymmetric Form

In cylindrical coordinates with axisymmetry:

```
∫₀^{z_max} ∫₀^{r_max} ε (∂V/∂r · ∂v/∂r + ∂V/∂z · ∂v/∂z) · 2πr dr dz
```

The factor `2πr` comes from the volume element in cylindrical coordinates when integrating over the full azimuthal angle.

### Derived Quantities

After solving for `V`:

**Electric field:**
```
E⃗ = -∇V = -(∂V/∂r r̂ + ∂V/∂z ẑ)
```

**Energy density:**
```
u = ½ ε |E⃗|²
```

### Boundary Conditions

- **Dirichlet:** `V = V₀` on `∂Ω_D` (e.g., electrode surfaces)
- **Neumann:** `n̂ · (ε ∇V) = g` on `∂Ω_N`. The configured `value`
  is this outward natural flux and is added to the weak-form boundary RHS.
  A zero value is the implicit natural condition and requires no assembled term.
- **Robin:** Reserved in the input schema but not yet implemented by the solvers.

For axisymmetric problems, a nonzero Neumann load is integrated with the
meridional boundary measure `2πr ds` (the implementation omits the common
global `2π` factor consistently from every term).

## 2. Magnetostatics

### Strong Form

The magnetostatic problem solves for the magnetic vector potential `A(r,z)`. In 2D axisymmetry, only the azimuthal component `A_φ` is non-zero:

```
∇ × (ν ∇ × A⃗) = J⃗
```

where:
- `A⃗ = A_φ(r,z) φ̂` is the magnetic vector potential [Wb/m or T·m]
- `ν = 1/μ` is the reluctivity [H⁻¹/m]
- `μ = μ_r μ₀` is the permeability [H/m]
- `μ₀ = 4π × 10⁻⁷ H/m` is the permeability of free space
- `J⃗ = J_φ(r,z) φ̂` is the current density [A/m²]

### Weak Form

Find `A_φ ∈ H¹` such that for all test functions `v`:

```
∫_Ω ν (∇ × A⃗) · (∇ × v) dΩ = ∫_Ω J⃗ · v dΩ
```

### Axisymmetric Form

The curl of the azimuthal vector potential is:

```
∇ × A⃗ = ∇ × (A_φ φ̂) = -∂A_φ/∂z r̂ + (1/r ∂(rA_φ)/∂r) ẑ
                       = -∂A_φ/∂z r̂ + (A_φ/r + ∂A_φ/∂r) ẑ
```

The weak form becomes:

```
∫₀^{z_max} ∫₀^{r_max} ν [(∂A_φ/∂r)(∂v/∂r) + (∂A_φ/∂z)(∂v/∂z) + A_φ·v/r²] · 2πr dr dz
     = ∫₀^{z_max} ∫₀^{r_max} J_φ · v · 2πr dr dz
```

### Singularity at r=0

Near the axis (`r → 0`), the term `A_φ/r` is potentially singular. However, for physical solutions, `A_φ → 0` as `r → 0`, and the limit:

```
lim_{r→0} A_φ/r = ∂A_φ/∂r|_{r=0}
```

exists and is finite. Our integrators handle this by substituting the derivative for the ratio when `r < tolerance`.

### Derived Quantities

**Magnetic flux density:**
```
B⃗ = ∇ × A⃗ = [-∂A_φ/∂z, (A_φ/r + ∂A_φ/∂r), 0] in (r, z, φ) components
```

**Magnetic field:**
```
H⃗ = ν B⃗
```

**Energy density:**
```
u = ½ B⃗ · H⃗ = ½ ν |B⃗|²
```

### Boundary Conditions

- **Dirichlet:** `A_φ = A₀` on `∂Ω_D` (e.g., far-field boundaries, symmetry planes)
- **Neumann:** The configured `value` is the outward natural flux `g`; zero is
  implicit and nonzero data contributes to the boundary RHS. In the planar
  scalar formulation, `g = ν n̂ · ∇A_z`. In the axisymmetric curl-curl
  formulation,
  `g = ν[n_r(∂A_φ/∂r + A_φ/r) + n_z ∂A_φ/∂z]`. The `A_φ/r` contribution is
  part of the radial natural flux and cannot be replaced by `n̂ · (ν∇A_φ)`.
- **Robin:** Reserved in the input schema but not yet implemented.

For an axisymmetric magnetic problem that reaches `r = 0`, regularity requires
`A_φ = 0` on the axis. The solver detects that boundary and applies this
essential constraint automatically; it is not a natural Neumann condition.

## 3. Magnetoquasistatics (Eddy Currents)

### Strong Form

The time-harmonic magnetoquasistatic problem includes eddy current effects:

```
∇ × (ν ∇ × A⃗) + jω σ A⃗ = J⃗_source
```

where:
- `ω = 2πf` is the angular frequency [rad/s]
- `σ` is the electrical conductivity [S/m]
- `J⃗_source` is the applied current density [A/m²]
- `j = √(-1)` is the imaginary unit

For time-harmonic fields (`e^{jωt}`), `A⃗` becomes complex: `A⃗ = A⃗_real + j A⃗_imag`.

### Weak Form

Find complex `A_φ ∈ H¹` such that for all test functions `v`:

```
∫_Ω ν (∇ × A⃗) · (∇ × v*) dΩ + jω ∫_Ω σ A⃗ · v* dΩ = ∫_Ω J⃗_source · v* dΩ
```

where `v*` is the complex conjugate of the test function.

### Axisymmetric Form

This becomes a complex-valued system:

```
[K + jωM] {A} = {F}
```

where:
- `K` is the stiffness matrix (curl-curl term)
- `M` is the mass matrix (conductivity term)
- `{F}` is the source term

In integral form:

```
∫_Ω ν (∇ × A⃗_real) · (∇ × v) · 2πr dr dz
    + jω ∫_Ω σ A⃗_real · v · 2πr dr dz (real part)

∫_Ω ν (∇ × A⃗_imag) · (∇ × v) · 2πr dr dz
    + jω ∫_Ω σ A⃗_imag · v · 2πr dr dz (imaginary part)
```

### Physical Interpretation

The real and imaginary parts of `A_φ` represent:
- **Real part:** Component in phase with the driving current
- **Imaginary part:** Component 90° out of phase (eddy current losses)

### Derived Quantities

**Complex magnetic flux density:**
```
B⃗ = B⃗_real + j B⃗_imag = ∇ × A⃗
```

**RMS magnitude:**
```
|B⃗| = √(|B⃗_real|² + |B⃗_imag|²)
```

**Time-domain fields:**
```
B⃗(t) = Re{B⃗ e^{jωt}} = B⃗_real cos(ωt) - B⃗_imag sin(ωt)
```

**Power loss density (Joule heating):**
```
P = ½ σ ω² |A⃗|²
```

### Boundary Conditions

Same essential/natural split as magnetostatics:
- **Dirichlet:** The configured real value constrains `A_φ`; the imaginary
  boundary value is currently zero.
- **Neumann:** The configured real outward natural flux is assembled into the
  real field RHS. A zero value remains implicit.
- **Robin:** Reserved in the input schema but not yet implemented.

## Finite Element Discretization

All three formulations use H¹-conforming (Lagrange) finite elements:

- **Basis functions:** `φ_i(r, z)` with `C⁰` continuity
- **Degrees of freedom:** Nodal values
- **Integration:** Gauss quadrature with elevated order to handle the `1/r` and `2πr` terms

### Special Considerations for Axisymmetry

1. **Axis handling:** Near `r = 0`, use L'Hôpital's rule to evaluate `A/r` limits
2. **Integration weight:** Include `2πr` factor in all volume integrals
3. **Boundary conditions:** Electrostatic axis symmetry is natural, while the
   magnetic `A_φ` formulation automatically enforces `A_φ = 0` on `r = 0`
