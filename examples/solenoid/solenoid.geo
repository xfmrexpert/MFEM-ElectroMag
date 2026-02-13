// Gmsh script for axisymmetric solenoid
// Axisymmetric geometry: r-z plane (x=r, y=z in Gmsh)
//
// Geometry:
//   - Solenoid coil: inner radius 0.05 m, outer radius 0.08 m, length 0.2 m
//   - Air core and surrounding region
//
// Material attributes:
//   1 = Air (core and surrounding)
//   2 = Coil (current source: J = 166,667 A/m^2)
//
// Boundary attributes:
//   1 = Far field (Dirichlet: A_phi = 0)

SetFactory("OpenCASCADE");

// Geometry parameters (all in meters)
coil_r_inner = 0.05;     // Inner radius
coil_r_outer = 0.08;     // Outer radius
coil_length = 0.2;       // Coil length
far_field_r = 0.40;      // Far field radius (5x coil outer radius)
far_field_z = 0.50;      // Far field half-height (2.5x coil length)

// Mesh parameters
lc_coil = 0.005;         // 5 mm in coil
lc_core = 0.005;         // 5 mm in air core
lc_air = 0.01;           // 10 mm in air near coil
lc_far = 0.05;           // 50 mm at far field

// Coil centered vertically at z=0
coil_z_start = -coil_length / 2;
coil_z_end = coil_length / 2;

// Air core (inside coil)
Rectangle(1) = {0, coil_z_start, 0, coil_r_inner, coil_length};

// Coil region
Rectangle(2) = {coil_r_inner, coil_z_start, 0, coil_r_outer - coil_r_inner, coil_length};

// Air outside coil (radial)
Rectangle(3) = {coil_r_outer, coil_z_start, 0, far_field_r - coil_r_outer, coil_length};

// Air above coil
// Inner region (r < coil_inner)
Rectangle(4) = {0, coil_z_end, 0, coil_r_inner, far_field_z - coil_z_end};

// Middle region (coil_inner < r < coil_outer)
Rectangle(5) = {coil_r_inner, coil_z_end, 0, coil_r_outer - coil_r_inner, far_field_z - coil_z_end};

// Outer region (r > coil_outer)
Rectangle(6) = {coil_r_outer, coil_z_end, 0, far_field_r - coil_r_outer, far_field_z - coil_z_end};

// Air below coil
// Inner region (r < coil_inner)
Rectangle(7) = {0, -far_field_z, 0, coil_r_inner, far_field_z + coil_z_start};

// Middle region (coil_inner < r < coil_outer)
Rectangle(8) = {coil_r_inner, -far_field_z, 0, coil_r_outer - coil_r_inner, far_field_z + coil_z_start};

// Outer region (r > coil_outer)
Rectangle(9) = {coil_r_outer, -far_field_z, 0, far_field_r - coil_r_outer, far_field_z + coil_z_start};

// Define physical surfaces (material attributes)
Physical Surface("Coil", 2) = {2};  // Coil with current density
Physical Surface("Air", 1) = {1, 3, 4, 5, 6, 7, 8, 9};  // All air regions (including core)

// Define physical curves for boundaries
// Far field boundary (outer edges)
Physical Curve("FarField", 1) = {
    26,  // Right edge bottom (r=far_field_r)
    30,  // Right edge middle bottom
    11,  // Right edge beside coil
    22,  // Right edge middle top
    18,  // Right edge top
    19,  // Top edge (z=far_field_z)
    25   // Bottom edge (z=-far_field_z)
};

// Mesh refinement
// Fine mesh in coil
Field[1] = Distance;
Field[1].SurfacesList = {2};  // Distance from coil

Field[2] = Threshold;
Field[2].InField = 1;
Field[2].SizeMin = lc_coil;
Field[2].SizeMax = lc_air;
Field[2].DistMin = 0;
Field[2].DistMax = 0.05;  // Grade over 5 cm

// Fine mesh in core
Field[3] = Distance;
Field[3].SurfacesList = {1};  // Distance from core

Field[4] = Threshold;
Field[4].InField = 3;
Field[4].SizeMin = lc_core;
Field[4].SizeMax = lc_air;
Field[4].DistMin = 0;
Field[4].DistMax = 0.05;

// Coarse mesh at far field
Field[5] = Distance;
Field[5].SurfacesList = {7, 9, 4, 6};  // Distance from far regions

Field[6] = Threshold;
Field[6].InField = 5;
Field[6].SizeMin = lc_air;
Field[6].SizeMax = lc_far;
Field[6].DistMin = 0.1;
Field[6].DistMax = 0.2;

// Combine fields
Field[7] = Min;
Field[7].FieldsList = {2, 4, 6};
Background Field = 7;

// Mesh settings
Mesh.Algorithm = 6;        // Frontal-Delaunay
Mesh.RecombineAll = 1;     // Recombine triangles into quads
Mesh.RecombineAll = 1;
Mesh.RecombinationAlgorithm = 2;
Mesh.ElementOrder = 1;     // Linear elements
Mesh.SecondOrderIncomplete = 0;

// Generate 2D mesh
Mesh 2;

// Optimize mesh quality
Mesh.Optimize = 1;
Mesh.OptimizeNetgen = 1;
