// Single circular current loop - axisymmetric test case
// Loop at radius a = 0.1 m, very thin cross-section

SetFactory("OpenCASCADE");

// Parameters
a = 0.1;              // Loop radius [m]
delta = 0.002;        // Loop thickness [m] (very thin)
domain_r = 0.5;       // Outer boundary radius
domain_z = 0.5;       // Domain half-height

// Mesh size
lc_loop = 0.0005;     // Very fine mesh in loop (4× refinement)
lc_domain = 0.02;     // Finer mesh in domain

// Current loop (thin annulus)
Rectangle(1) = {a - delta/2, -delta/2, 0, delta, delta};

// Air domain
Rectangle(2) = {0, -domain_z, 0, domain_r, 2*domain_z};

BooleanDifference(3) = { Surface{2}; Delete; }{ Surface{1}; };

eps = 1e-6;

// Physical surfaces. Identify the loop by its bounding box so the assignment
// does not depend on ids that shift after the boolean operation.
loop_surfaces() = Surface In BoundingBox{
    a - delta/2 - eps, -delta/2 - eps, -eps,
    a + delta/2 + eps, delta/2 + eps, eps};
Physical Surface("Loop", 2) = {loop_surfaces()};

all_surfaces() = Surface{:};
air_surfaces() = {};
For i In {0 : #all_surfaces() - 1}
    is_loop = 0;
    For j In {0 : #loop_surfaces() - 1}
        If (all_surfaces(i) == loop_surfaces(j))
            is_loop = 1;
        EndIf
    EndFor
    If (is_loop == 0)
        air_surfaces() += all_surfaces(i);
    EndIf
EndFor
Physical Surface("Air", 1) = {air_surfaces()};

// Physical boundaries, selected geometrically rather than by hard-coded ids.
axis_curves() = Curve In BoundingBox{
    -eps, -domain_z - eps, -eps,
    eps, domain_z + eps, eps};
right_curves() = Curve In BoundingBox{
    domain_r - eps, -domain_z - eps, -eps,
    domain_r + eps, domain_z + eps, eps};
bottom_curves() = Curve In BoundingBox{
    -eps, -domain_z - eps, -eps,
    domain_r + eps, -domain_z + eps, eps};
top_curves() = Curve In BoundingBox{
    -eps, domain_z - eps, -eps,
    domain_r + eps, domain_z + eps, eps};

// The symmetry axis gets its own attribute so A_phi = 0 there is distinct
// from the far-field Dirichlet condition.
Physical Curve("Outer", 1) = {right_curves(), bottom_curves(), top_curves()};
Physical Curve("Axis", 2) = {axis_curves()};

// Mesh size control
// Fine mesh near loop
Field[1] = Distance;
Field[1].SurfacesList = {loop_surfaces()};  // Distance from loop

Field[2] = Threshold;
Field[2].InField = 1;
Field[2].SizeMin = lc_loop;
Field[2].SizeMax = lc_domain;
Field[2].DistMin = 0.0;
Field[2].DistMax = 0.05;  // Grade over 5 cm

Background Field = 2;

// Mesh settings
Mesh.Algorithm = 6;
Mesh.ElementOrder = 1;

Mesh 2;
