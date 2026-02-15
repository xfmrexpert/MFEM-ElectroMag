// Single circular current loop - axisymmetric test case
// Loop at radius a = 0.1 m, very thin cross-section

SetFactory("OpenCASCADE");

// Parameters
a = 0.1;              // Loop radius [m]
delta = 0.002;        // Loop thickness [m] (very thin)
domain_r = 0.5;       // Outer boundary radius
domain_z = 0.5;       // Domain half-height

// Mesh size
lc_loop = 0.002;      // Fine mesh in loop
lc_domain = 0.05;     // Coarse mesh in domain

// Current loop (thin annulus)
Rectangle(1) = {a - delta/2, -delta/2, 0, delta, delta};

// Air domain
Rectangle(2) = {0, -domain_z, 0, a - delta/2, 2*domain_z};        // Inner region
Rectangle(3) = {a + delta/2, -domain_z, 0, domain_r - a - delta/2, 2*domain_z}; // Outer region
Rectangle(4) = {a - delta/2, -domain_z, 0, delta, delta/2 + domain_z};  // Below loop
Rectangle(5) = {a - delta/2, delta/2, 0, delta, delta/2 + domain_z};    // Above loop

// Fragment all rectangles to ensure proper connectivity at interfaces
BooleanFragments{ Surface{1, 2, 3, 4, 5}; Delete; }{}

// Physical surfaces (note: surface IDs will change after fragmentation)
// Will need to identify correct IDs after fragmentation
// For now, mark all surfaces and assign attributes based on geometry
Physical Surface("Loop", 2) = {1};  // Will need to update this
Physical Surface("Air", 1) = {2, 3, 4, 5};  // Will need to update this

// Physical boundaries (outer edges)
Physical Curve("Outer", 1) = {15, 19, 12, 7};  // Outer boundary (A=0)

// Mesh settings
Mesh.Algorithm = 6;
Mesh.RecombineAll = 1;
Mesh.ElementOrder = 1;

Mesh 2;
