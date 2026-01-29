// ---------------------------------------------------------
// test_loop.geo
// Axisymmetric model of a square wire loop in a vacuum box.
// ---------------------------------------------------------

// 1. Parameters (Units: meters)
R_dom  = 2.0;   // Radius of the air domain
Z_dom  = 1.0;   // Half-height of the air domain (Total height = 2.0)

R_coil = 0.5;   // Radius of the loop center
W_coil = 0.1;   // Width of the wire square cross-section
H_coil = 0.1;   // Height of the wire square cross-section

// 2. Mesh Sizing
lc_far  = 0.2;   // Coarse mesh at far-field
lc_coil = 0.02;  // Fine mesh inside the wire
lc_axis = 0.1;   // Moderate mesh along the axis

// 3. Points

// -- Air Box Corners --
Point(1) = {0,      -Z_dom, 0, lc_axis};  // Bottom-Left (Axis)
Point(2) = {R_dom,  -Z_dom, 0, lc_far};   // Bottom-Right
Point(3) = {R_dom,   Z_dom, 0, lc_far};   // Top-Right
Point(4) = {0,       Z_dom, 0, lc_axis};  // Top-Left (Axis)

// -- Coil Corners --
// Calculated relative to R_coil
x_min = R_coil - W_coil/2.0;
x_max = R_coil + W_coil/2.0;
y_min = -H_coil/2.0;
y_max =  H_coil/2.0;

Point(5) = {x_min, y_min, 0, lc_coil};
Point(6) = {x_max, y_min, 0, lc_coil};
Point(7) = {x_max, y_max, 0, lc_coil};
Point(8) = {x_min, y_max, 0, lc_coil};

// 4. Lines

// -- Air Box --
Line(1) = {1, 2}; // Bottom
Line(2) = {2, 3}; // Right (FarField)
Line(3) = {3, 4}; // Top
Line(4) = {4, 1}; // Left (Axis of Symmetry)

// -- Coil --
Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 5};

// 5. Surfaces

// Loop for Air (Outer Box)
Curve Loop(1) = {1, 2, 3, 4};

// Loop for Coil
Curve Loop(2) = {5, 6, 7, 8};

// Surface 1: Vacuum (Air Box minus Coil hole)
Plane Surface(1) = {1, 2};

// Surface 2: Coil (The hole filled back in)
Plane Surface(2) = {2};

// 6. Physical Groups
// These IDs match the "attributes" in your JSON input file.

// -- Volumes (2D Surfaces) --
Physical Surface(1) = {1}; // Name: "Vacuum"
Physical Surface(2) = {2}; // Name: "Coil_Source"

// -- Boundaries (1D Curves) --
// Assigning IDs 1, 2, 3, 4 allows the JSON to group them easily.
Physical Curve(1) = {1}; // Bottom
Physical Curve(2) = {2}; // Right (FarField)
Physical Curve(3) = {3}; // Top
Physical Curve(4) = {4}; // Axis (Critical for MQS!)

// 7. Compatibility Settings
// MFEM works best with the legacy Gmsh 2.2 format
Mesh.MshFileVersion = 2.2;
