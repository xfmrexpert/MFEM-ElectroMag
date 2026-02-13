// Gmsh script for eddy current problem - conducting cylinder in AC field
// Axisymmetric geometry: r-z plane (x=r, y=z in Gmsh)
//
// Geometry:
//   - Conducting cylinder: radius 0.05 m, length 0.1 m
//   - Coil region: inner radius 0.06 m, outer radius 0.08 m, centered on cylinder
//   - Air region extending to far field
//
// Material attributes:
//   1 = Air
//   2 = Aluminum conductor
//   3 = Coil (current source)
//
// Boundary attributes:
//   1 = Far field (Dirichlet: A_phi = 0)

SetFactory("OpenCASCADE");

// Geometry parameters (all in meters)
cyl_radius = 0.05;      // Cylinder radius
cyl_length = 0.1;       // Cylinder length
coil_r_inner = 0.06;    // Coil inner radius
coil_r_outer = 0.08;    // Coil outer radius
coil_length = 0.12;     // Coil length (slightly longer than cylinder)
far_field_r = 0.25;     // Far field radius
far_field_z = 0.30;     // Far field half-height

// Mesh parameters
lc_conductor_surface = 0.0005;  // 0.5 mm at conductor surface (< skin_depth/3)
lc_conductor_bulk = 0.005;      // 5 mm in conductor interior
lc_coil = 0.005;                // 5 mm in coil
lc_air = 0.010;                 // 10 mm in air near conductor
lc_far = 0.030;                 // 30 mm at far field

// Cylinder center at origin, extending from z=0 to z=cyl_length
Rectangle(1) = {0, 0, 0, cyl_radius, cyl_length};

// Coil region (annular region around cylinder)
coil_z_start = (cyl_length - coil_length) / 2;  // Center coil on cylinder
coil_z_end = coil_z_start + coil_length;
Rectangle(2) = {coil_r_inner, coil_z_start, 0, coil_r_outer - coil_r_inner, coil_length};

// Inner air region (between cylinder and coil)
Rectangle(3) = {cyl_radius, 0, 0, coil_r_inner - cyl_radius, cyl_length};

// Air around coil (extends to full cylinder length)
Rectangle(4) = {coil_r_outer, 0, 0, far_field_r - coil_r_outer, cyl_length};

// Air above cylinder
Rectangle(5) = {0, cyl_length, 0, far_field_r, far_field_z - cyl_length};

// Air below cylinder
Rectangle(6) = {0, -far_field_z, 0, far_field_r, far_field_z};

// Air beside coil (top part)
Rectangle(7) = {coil_r_inner, cyl_length, 0, coil_r_outer - coil_r_inner, coil_z_end - cyl_length};

// Air beside coil (bottom part)
Rectangle(8) = {coil_r_inner, coil_z_start, 0, coil_r_outer - coil_r_inner, -coil_z_start};

// Define physical surfaces (material attributes)
Physical Surface("Conductor", 2) = {1};  // Aluminum cylinder
Physical Surface("Coil", 3) = {2};        // Coil (current source)
Physical Surface("Air", 1) = {3, 4, 5, 6, 7, 8};  // All air regions

// Define physical curves for boundaries
// Far field boundary (outer edges)
Physical Curve("FarField", 1) = {
    9,   // Right edge of rectangle 6 (r=far_field_r, z=-far_field_z to 0)
    10,  // Right edge of rectangle 4 (r=far_field_r, z=0 to cyl_length)
    14,  // Right edge of rectangle 5 (r=far_field_r, z=cyl_length to far_field_z)
    15,  // Top edge of rectangle 5 (z=far_field_z, r=0 to far_field_r)
    5    // Bottom edge of rectangle 6 (z=-far_field_z, r=0 to far_field_r)
};

// Mesh refinement
// Fine mesh at conductor surface (skin depth resolution)
MeshSize {1, 2} = lc_conductor_surface;  // Vertices of conductor surface

// Graded mesh in conductor
Field[1] = Distance;
Field[1].SurfacesList = {1};  // Distance from conductor

Field[2] = Threshold;
Field[2].InField = 1;
Field[2].SizeMin = lc_conductor_surface;
Field[2].SizeMax = lc_conductor_bulk;
Field[2].DistMin = 0;
Field[2].DistMax = cyl_radius * 0.5;  // Grade over half the radius

// Medium mesh in coil region
Field[3] = Distance;
Field[3].SurfacesList = {2};  // Distance from coil

Field[4] = Threshold;
Field[4].InField = 3;
Field[4].SizeMin = lc_coil;
Field[4].SizeMax = lc_air;
Field[4].DistMin = 0;
Field[4].DistMax = 0.02;

// Coarse mesh at far field
Field[5] = Distance;
Field[5].SurfacesList = {5, 6};  // Distance from far field regions

Field[6] = Threshold;
Field[6].InField = 5;
Field[6].SizeMin = lc_air;
Field[6].SizeMax = lc_far;
Field[6].DistMin = 0.05;
Field[6].DistMax = 0.15;

// Combine fields
Field[7] = Min;
Field[7].FieldsList = {2, 4, 6};
Background Field = 7;

// Mesh settings
Mesh.Algorithm = 6;        // Frontal-Delaunay for quads
Mesh.RecombineAll = 1;     // Recombine triangles into quads
Mesh.RecombinationAlgorithm = 2;  // Simple full-quad recombination
Mesh.ElementOrder = 1;     // Linear elements (MFEM will handle higher order)
Mesh.SecondOrderIncomplete = 0;

// Generate 2D mesh
Mesh 2;

// Optional: Optimize mesh quality
Mesh.Optimize = 1;
Mesh.OptimizeNetgen = 1;
