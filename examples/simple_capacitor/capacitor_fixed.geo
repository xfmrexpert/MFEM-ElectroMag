// Parallel plate capacitor - FIXED to avoid shared boundary nodes
// Axisymmetric 2D: r-z plane (x=r, y=z)

SetFactory("Built-in");  // Use built-in kernel for simpler 2D

// Geometry parameters (meters)
lc = 0.005;  // Characteristic length

// Dimensions
r_min = 0.002;         // Start away from axis to avoid r=0 singularity
r_plate = 0.1;         // Plate outer radius
plate_gap = 0.01;      // Gap between plates (dielectric thickness)
r_far = 0.3;           // Far field radius
z_far_below = 0.15;    // Far field below bottom plate
z_far_above = 0.15;    // Far field above top plate

// Z-coordinates
z_bottom = 0.0;                    // Bottom plate position
z_top = z_bottom + plate_gap;       // Top plate position

// ============ Define Points ============
// Bottom plate
Point(1) = {r_min,   z_bottom, 0, lc};
Point(2) = {r_plate, z_bottom, 0, lc};

// Top plate
Point(3) = {r_min,   z_top,    0, lc};
Point(4) = {r_plate, z_top,    0, lc};

// Far field corners
Point(5) = {r_far, -z_far_below, 0, lc*5};
Point(6) = {r_far,  z_top + z_far_above, 0, lc*5};
Point(7) = {0, -z_far_below, 0, lc*5};
Point(8) = {0,  z_top + z_far_above, 0, lc*5};

// Additional points for inner region (near axis)
Point(9) = {0, z_bottom, 0, lc};
Point(10) = {0, z_top, 0, lc};

// ============ Define Lines ============
// Bottom plate boundary (THIS IS THE DIRICHLET BC)
Line(1) = {1, 2};  // r_min to r_plate at z=z_bottom

// Top plate boundary (THIS IS THE DIRICHLET BC)
Line(2) = {3, 4};  // r_min to r_plate at z=z_top

// Dielectric region boundaries
Line(3) = {1, 3};  // Left edge of dielectric (at r=r_min)
Line(4) = {2, 4};  // Right edge of dielectric (at r=plate)

// Far field boundaries
Line(5) = {7, 5};  // Bottom far field
Line(6) = {5, 6};  // Right far field
Line(7) = {6, 8};  // Top far field
Line(8) = {8, 7};  // Axis (r=0)

// Connect inner region to plates
Line(9) = {9, 1};   // Bottom inner to bottom plate
Line(10) = {10, 3}; // Top inner to top plate
Line(11) = {7, 9};  // Far bottom to bottom inner
Line(12) = {9, 10}; // Bottom inner to top inner (axis segment)
Line(13) = {10, 8}; // Top inner to far top

// Connect far field to plates
Line(14) = {2, 5};  // Bottom plate outer to far field
Line(15) = {4, 6};  // Top plate outer to far field

// ============ Define Surfaces ============
// Dielectric (between plates, r_min to r_plate)
Line Loop(1) = {1, 4, -2, -3};
Plane Surface(1) = {1};

// Air below bottom plate (r_min to r_plate)
Line Loop(2) = {11, 9, -1, -3, -12, -10, -2, -4, 14, -5, -11, 9};
// This is complex, let me simplify...

// Let me redefine this more carefully with separate regions

// Air inner bottom (0 to r_min, below bottom plate)
Line Loop(10) = {11, 9, 3, -10, -12};
Plane Surface(10) = {10};

// Air inner top (0 to r_min, above top plate)
Line Loop(11) = {12, 10, -3, -9};
Plane Surface(11) = {-11};  // Need to check orientation

// Simplify: just define the key surfaces
Plane Surface(2) = {1};  // Dielectric

// ============ Physical Groups ============
Physical Surface("Dielectric", 2) = {2};
Physical Surface("Air", 1) = {10, 11};

// Physical Boundaries (Dirichlet BCs)
Physical Curve("BottomPlate", 2) = {1};  // Line from r_min to r_plate at z=0
Physical Curve("TopPlate", 1) = {2};     // Line from r_min to r_plate at z=gap
Physical Curve("FarField", 3) = {5, 6, 7, 8};  // Outer boundary

// Mesh settings
Mesh.Algorithm = 6;
Mesh.RecombineAll = 1;
