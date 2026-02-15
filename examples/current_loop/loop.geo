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

// Physical surfaces (note: surface IDs will change after fragmentation)
// Will need to identify correct IDs after fragmentation
// For now, mark all surfaces and assign attributes based on geometry
Physical Surface("Loop", 2) = {1};
Physical Surface("Air", 1) = {3};

// Physical boundaries (outer edges)
Physical Curve("Outer", 1) = {15, 19, 12, 7};  // Outer boundary (A=0)

// Mesh size control
// Fine mesh near loop
Field[1] = Distance;
Field[1].SurfacesList = {1};  // Distance from loop

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
