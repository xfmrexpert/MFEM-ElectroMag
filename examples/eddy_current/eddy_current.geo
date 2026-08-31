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
// Skin depth in aluminum (sigma = 3.5e7 S/m, mu_r = 1) at the highest swept
// frequency of 1 kHz is delta = sqrt(1/(pi*f*mu*sigma)) ~= 2.7 mm. Resolving
// with ~3 elements per skin depth means ~0.9 mm at the conductor surface;
// anything finer only inflates the system without improving accuracy.
lc_conductor_surface = 0.0009;  // 0.9 mm at conductor surface (~delta/3)
lc_conductor_bulk = 0.006;      // 6 mm in conductor interior
lc_coil = 0.006;                // 6 mm in coil
lc_air = 0.015;                 // 15 mm in air near conductor
lc_far = 0.040;                 // 40 mm at far field

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

eps = 1e-6;

// Define physical surfaces (material attributes).
// Identify each region by its bounding box so the assignment does not depend
// on ids that shift after the boolean operations.
conductor_surfaces() = Surface In BoundingBox{
    -eps, -eps, -eps,
    cyl_radius + eps, cyl_length + eps, eps};
coil_surfaces() = Surface In BoundingBox{
    coil_r_inner - eps, coil_z_start - eps, -eps,
    coil_r_outer + eps, coil_z_end + eps, eps};
Physical Surface("Conductor", 2) = {conductor_surfaces()};
Physical Surface("Coil", 3) = {coil_surfaces()};

all_surfaces() = Surface{:};
air_surfaces() = {};
For i In {0 : #all_surfaces() - 1}
    is_other = 0;
    For j In {0 : #conductor_surfaces() - 1}
        If (all_surfaces(i) == conductor_surfaces(j))
            is_other = 1;
        EndIf
    EndFor
    For j In {0 : #coil_surfaces() - 1}
        If (all_surfaces(i) == coil_surfaces(j))
            is_other = 1;
        EndIf
    EndFor
    If (is_other == 0)
        air_surfaces() += all_surfaces(i);
    EndIf
EndFor
Physical Surface("Air", 1) = {air_surfaces()};

// Define physical curves for boundaries, selected geometrically rather than
// by hard-coded ids that are invalidated by the boolean operations.
axis_curves() = Curve In BoundingBox{
    -eps, -far_field_z - eps, -eps,
    eps, far_field_z + eps, eps};
right_curves() = Curve In BoundingBox{
    far_field_r - eps, -far_field_z - eps, -eps,
    far_field_r + eps, far_field_z + eps, eps};
bottom_curves() = Curve In BoundingBox{
    -eps, -far_field_z - eps, -eps,
    far_field_r + eps, -far_field_z + eps, eps};
top_curves() = Curve In BoundingBox{
    -eps, far_field_z - eps, -eps,
    far_field_r + eps, far_field_z + eps, eps};

// The symmetry axis needs its own attribute so A_phi = 0 is enforced there
// without contaminating the far-field boundary.
Physical Curve("FarField", 1) = {right_curves(), bottom_curves(), top_curves()};
Physical Curve("Axis", 2) = {axis_curves()};

// Mesh refinement
// Fine mesh at conductor surface (skin depth resolution)
MeshSize {1, 2} = lc_conductor_surface;  // Vertices of conductor surface

// Graded mesh in conductor
Field[1] = Distance;
Field[1].SurfacesList = {conductor_surfaces()};  // Distance from conductor

Field[2] = Threshold;
Field[2].InField = 1;
Field[2].SizeMin = lc_conductor_surface;
Field[2].SizeMax = lc_conductor_bulk;
Field[2].DistMin = 0;
Field[2].DistMax = cyl_radius * 0.5;  // Grade over half the radius

// Medium mesh in coil region
Field[3] = Distance;
Field[3].SurfacesList = {coil_surfaces()};  // Distance from coil

Field[4] = Threshold;
Field[4].InField = 3;
Field[4].SizeMin = lc_coil;
Field[4].SizeMax = lc_air;
Field[4].DistMin = 0;
Field[4].DistMax = 0.02;

// Coarse mesh at far field
Field[5] = Distance;
Field[5].CurvesList = {right_curves(), bottom_curves(), top_curves()};

Field[6] = Threshold;
Field[6].InField = 5;
// Distance is measured from the far-field boundary, so the coarse size applies
// at distance 0 and grades back down toward the interior.
Field[6].SizeMin = lc_far;
Field[6].SizeMax = lc_air;
Field[6].DistMin = 0.0;
Field[6].DistMax = 0.15;

// Combine fields
Field[7] = Min;
Field[7].FieldsList = {2, 4, 6};
Background Field = 7;

// Mesh settings
Mesh.Algorithm = 6;        // Frontal-Delaunay
// Triangles rather than recombined quads: the strongly graded skin-depth mesh
// makes recombination unreliable here, and MFEM handles triangles natively.
Mesh.RecombineAll = 0;
Mesh.ElementOrder = 1;     // Linear elements (MFEM will handle higher order)
Mesh.SecondOrderIncomplete = 0;

// Generate 2D mesh
Mesh 2;

// Optional: Optimize mesh quality
Mesh.Optimize = 1;
Mesh.OptimizeNetgen = 1;
