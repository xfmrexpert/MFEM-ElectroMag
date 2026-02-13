// Simple parallel plate capacitor - no shared boundary nodes
// Axisymmetric geometry: r-z plane

SetFactory("OpenCASCADE");

// Parameters
r_inner = 0.005;  // Inner radius (avoid r=0)
r_plate = 0.1;    // Plate radius
gap = 0.01;        // Gap between plates
r_outer = 0.3;     // Far field
z_extent = 0.15;   // Far field extent

// Z coordinates
z_bottom_plate = 0.0;
z_top_plate = gap;

// Create regions using boxes
// Bottom plate region
Box(1) = {r_inner, -z_extent, 0, r_plate - r_inner, z_extent, 1};

// Gap (dielectric) region
Box(2) = {r_inner, z_bottom_plate, 0, r_plate - r_inner, gap, 1};

// Top plate region
Box(3) = {r_inner, z_top_plate, 0, r_plate - r_inner, z_extent, 1};

// Outer air region
Box(4) = {r_plate, -z_extent, 0, r_outer - r_plate, gap + 2*z_extent, 1};

// Inner axis region
Box(5) = {0, -z_extent, 0, r_inner, gap + 2*z_extent, 1};

// Extrude to create volumes (needed for OpenCASCADE)
// But we only want 2D surfaces for the 2D mesh

// Physical groups
Physical Surface("Air", 1) = {1, 3, 4, 5};
Physical Surface("Dielectric", 2) = {2};

// Boundaries - define explicitly by creating curves
// Bottom plate: line from (r_inner, 0) to (r_plate, 0)
Point(100) = {r_inner, z_bottom_plate, 0};
Point(101) = {r_plate, z_bottom_plate, 0};
Line(100) = {100, 101};
Physical Curve("BottomPlate", 2) = {100};

// Top plate: line from (r_inner, gap) to (r_plate, gap)
Point(102) = {r_inner, z_top_plate, 0};
Point(103) = {r_plate, z_top_plate, 0};
Line(101) = {102, 103};
Physical Curve("TopPlate", 1) = {101};

// Far field: outer boundary
Point(104) = {r_outer, -z_extent, 0};
Point(105) = {r_outer, z_top_plate + z_extent, 0};
Line(102) = {104, 105};
Point(106) = {0, -z_extent, 0};
Point(107) = {r_outer, -z_extent, 0};
Line(103) = {106, 107};
Point(108) = {0, z_top_plate + z_extent, 0};
Point(109) = {r_outer, z_top_plate + z_extent, 0};
Line(104) = {108, 109};
Physical Curve("FarField", 3) = {102, 103, 104};

// Mesh settings
Mesh.CharacteristicLengthMin = 0.002;
Mesh.CharacteristicLengthMax = 0.03;
Mesh.Algorithm = 6;
Mesh.RecombineAll = 1;
