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

// Make the independently created rectangles share edges and nodes so the
// mesh is conforming across every material interface.
BooleanFragments{ Surface{1:9}; Delete; }{}

// Define physical surfaces (material attributes)
// The coil is the only non-air region; identify it by its bounding box so the
// attribute assignment survives any renumbering done by the boolean operation.
eps = 1e-6;
coil_surfaces() = Surface In BoundingBox{
    coil_r_inner - eps, coil_z_start - eps, -eps,
    coil_r_outer + eps, coil_z_end + eps, eps};
Physical Surface("Coil", 2) = {coil_surfaces()};

all_surfaces() = Surface{:};
air_surfaces() = {};
For i In {0 : #all_surfaces() - 1}
    is_coil = 0;
    For j In {0 : #coil_surfaces() - 1}
        If (all_surfaces(i) == coil_surfaces(j))
            is_coil = 1;
        EndIf
    EndFor
    If (is_coil == 0)
        air_surfaces() += all_surfaces(i);
    EndIf
EndFor
Physical Surface("Air", 1) = {air_surfaces()};

// Define physical curves for boundaries.
// Select geometrically rather than by hard-coded ids, which are not stable
// across Gmsh versions or boolean operations.
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
// Fine mesh in coil
Field[1] = Distance;
Field[1].SurfacesList = {coil_surfaces()};  // Distance from coil

Field[2] = Threshold;
Field[2].InField = 1;
Field[2].SizeMin = lc_coil;
Field[2].SizeMax = lc_far;
Field[2].DistMin = 0;
Field[2].DistMax = 0.05;  // Grade over 5 cm

// Fine mesh in core
core_surfaces() = Surface In BoundingBox{
    -eps, coil_z_start - eps, -eps,
    coil_r_inner + eps, coil_z_end + eps, eps};

Field[3] = Distance;
Field[3].SurfacesList = {core_surfaces()};  // Distance from core

Field[4] = Threshold;
Field[4].InField = 3;
Field[4].SizeMin = lc_core;
Field[4].SizeMax = lc_far;
Field[4].DistMin = 0;
Field[4].DistMax = 0.05;

// Coarse mesh at far field
Field[5] = Distance;
Field[5].CurvesList = {right_curves(), bottom_curves(), top_curves()};  // Distance from far field

Field[6] = Threshold;
Field[6].InField = 5;
// Distance is now measured from the far-field boundary, so the coarse size
// applies at distance 0 and grades back down toward the interior.
Field[6].SizeMin = lc_far;
Field[6].SizeMax = lc_air;
Field[6].DistMin = 0.0;
Field[6].DistMax = 0.2;

// Combine fields
Field[7] = Min;
Field[7].FieldsList = {2, 4, 6};
Background Field = 7;

// Mesh settings
Mesh.Algorithm = 6;        // Frontal-Delaunay
// Triangles are used rather than recombined quads: the graded sizes here make
// recombination fail on some subregions, and MFEM handles triangles natively.
Mesh.RecombineAll = 0;
Mesh.ElementOrder = 1;     // Linear elements
Mesh.SecondOrderIncomplete = 0;

// Generate 2D mesh
Mesh 2;

// Optimize mesh quality
Mesh.Optimize = 1;
Mesh.OptimizeNetgen = 1;
