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

// Create geometries
// Conductor (aluminum cylinder)
Rectangle(1) = {0, 0, 0, cyl_radius, cyl_length};

// Coil region (annular region around cylinder)
coil_z_start = (cyl_length - coil_length) / 2;  // Center coil on cylinder
coil_z_end = coil_z_start + coil_length;
Rectangle(2) = {coil_r_inner, coil_z_start, 0, coil_r_outer - coil_r_inner, coil_length};

// Air domain (entire computational domain)
Rectangle(3) = {0, -far_field_z, 0, far_field_r, 2*far_field_z};

// Use BooleanDifference to cut out conductor and coil from air
// This ensures proper connectivity at interfaces
BooleanDifference(4) = { Surface{3}; Delete; }{ Surface{1}; };
BooleanDifference(5) = { Surface{4}; Delete; }{ Surface{2}; };

// Define physical surfaces (material attributes)
// Note: Surface IDs change after Boolean operations
Physical Surface("Conductor", 2) = {1};  // Aluminum cylinder
Physical Surface("Coil", 3) = {2};        // Coil (current source)
Physical Surface("Air", 1) = {5};         // Air (after cutting out conductor and coil)

// Define physical curves for boundaries
// Note: Curve IDs will need to be identified after Boolean operations
// For now, mark the outer boundary of the final air domain
// The outer rectangle (3) has curves that form the far field boundary
// After Boolean operations, these curves remain on surface 5
// Typically: bottom, right, top edges (left edge at r=0 is axis)
Physical Curve("FarField", 1) = {9, 10, 11};  // Will need to verify these IDs

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
