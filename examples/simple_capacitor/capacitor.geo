// Gmsh script for parallel plate capacitor
// Axisymmetric geometry: r-z plane (x=r, y=z in Gmsh)
//
// Geometry:
//   - Two circular plates: radius 0.1 m, separation 0.01 m
//   - Dielectric between plates
//   - Air surrounding
//
// Material attributes:
//   1 = Air
//   2 = Dielectric (epsilon_r = 2.5)
//
// Boundary attributes:
//   1 = Top plate (V = 1000 V)
//   2 = Bottom plate (V = 0 V)
//   3 = Far field (V = 0 V)

SetFactory("OpenCASCADE");

// Geometry parameters (all in meters)
r_min = 0.001;           // Plate bore radius; the field domain includes the axis
plate_radius = 0.1;      // Plate outer radius
plate_thickness = 0.001; // Thin plate (1 mm)
separation = 0.01;       // Separation between plates (1 cm)
far_field_r = 0.3;       // Far field radius (3x plate radius)
far_field_z_top = 0.15;  // Far field distance above top plate
far_field_z_bot = 0.15;  // Far field distance below bottom plate

// Mesh parameters
lc_plate = 0.002;        // 2 mm at plates
lc_dielectric = 0.002;   // 2 mm in dielectric
lc_air = 0.01;           // 10 mm in air near plates
lc_far = 0.03;           // 30 mm at far field

// Z-coordinates
z_bot_plate = 0;
z_bot_plate_top = z_bot_plate + plate_thickness;
z_top_plate_bot = z_bot_plate_top + separation;
z_top_plate_top = z_top_plate_bot + plate_thickness;

// Bottom plate (conductor) - starts at r_min to avoid axis
Rectangle(1) = {r_min, z_bot_plate, 0, plate_radius - r_min, plate_thickness};

// Dielectric between plates
Rectangle(2) = {r_min, z_bot_plate_top, 0, plate_radius - r_min, separation};

// Top plate (conductor) - starts at r_min to avoid axis
Rectangle(3) = {r_min, z_top_plate_bot, 0, plate_radius - r_min, plate_thickness};

// Air near axis (r < r_min)
// Below bottom plate
Rectangle(11) = {0, -far_field_z_bot, 0, r_min, far_field_z_bot};
// Bottom plate height
Rectangle(12) = {0, z_bot_plate, 0, r_min, plate_thickness};
// Dielectric height
Rectangle(13) = {0, z_bot_plate_top, 0, r_min, separation};
// Top plate height
Rectangle(14) = {0, z_top_plate_bot, 0, r_min, plate_thickness};
// Above top plate
Rectangle(15) = {0, z_top_plate_top, 0, r_min, far_field_z_top - z_top_plate_top};

// Air inside plate radius (r_min < r < plate_radius)
// Below bottom plate
Rectangle(4) = {r_min, -far_field_z_bot, 0, plate_radius - r_min, far_field_z_bot};

// Above top plate
Rectangle(5) = {r_min, z_top_plate_top, 0, plate_radius - r_min, far_field_z_top - z_top_plate_top};

// Air outside capacitor radius (r > plate_radius)
// Bottom region
Rectangle(6) = {plate_radius, -far_field_z_bot, 0, far_field_r - plate_radius, far_field_z_bot};

// Region beside bottom plate
Rectangle(7) = {plate_radius, z_bot_plate, 0, far_field_r - plate_radius, plate_thickness};

// Region beside dielectric
Rectangle(8) = {plate_radius, z_bot_plate_top, 0, far_field_r - plate_radius, separation};

// Region beside top plate
Rectangle(9) = {plate_radius, z_top_plate_bot, 0, far_field_r - plate_radius, plate_thickness};

// Top region
Rectangle(10) = {plate_radius, z_top_plate_top, 0, far_field_r - plate_radius, far_field_z_top - z_top_plate_top};

BooleanFragments{ Surface{1:15}; Delete; }{}
bottom_plate_curves[] = Boundary{ Surface{1}; };
top_plate_curves[] = Boundary{ Surface{3}; };
bottom_plate_curves[] = Abs(bottom_plate_curves[]);
top_plate_curves[] = Abs(top_plate_curves[]);
Delete{ Surface{1, 3}; }

// Define physical surfaces (material attributes)
Physical Surface("Dielectric", 2) = {2, 13};  // Dielectric between plates (including axis region)
Physical Surface("Air", 1) = {4, 5, 6, 7, 8, 9, 10, 11, 12, 14, 15};  // All air regions

// Select curves by bounding box (with small tolerance)
tol = 0.0001;

far_field_curves[] = Curve In BoundingBox {far_field_r - tol, -far_field_z_bot - tol, -tol,
                                             far_field_r + tol, far_field_z_top + tol, tol};
far_field_top[] = Curve In BoundingBox {-tol, far_field_z_top - tol, -tol,
                                      far_field_r + tol, far_field_z_top + tol, tol};
far_field_bottom[] = Curve In BoundingBox {-tol, -far_field_z_bot - tol, -tol,
                                         far_field_r + tol, -far_field_z_bot + tol, tol};

// Define physical curves
Physical Curve("TopPlate", 1) = {top_plate_curves[]};
Physical Curve("BottomPlate", 2) = {bottom_plate_curves[]};
Physical Curve("FarField", 3) = {far_field_curves[], far_field_top[], far_field_bottom[]};

// Mesh refinement
// Fine mesh at plates and in dielectric
Field[1] = Distance;
Field[1].SurfacesList = {2, 13};

Field[2] = Threshold;
Field[2].InField = 1;
Field[2].SizeMin = lc_dielectric;
Field[2].SizeMax = lc_air;
Field[2].DistMin = 0;
Field[2].DistMax = 0.05;  // Grade over 5 cm

// Coarse mesh at far field
Field[3] = Distance;
Field[3].SurfacesList = {4, 5};  // Distance from far regions

Field[4] = Threshold;
Field[4].InField = 3;
Field[4].SizeMin = lc_air;
Field[4].SizeMax = lc_far;
Field[4].DistMin = 0.05;
Field[4].DistMax = 0.1;

// Refinement near plate edges (fringing field)
Field[5] = Distance;
Field[5].CurvesList = {top_plate_curves[], bottom_plate_curves[]};

Field[6] = Threshold;
Field[6].InField = 5;
Field[6].SizeMin = lc_plate;
Field[6].SizeMax = lc_air;
Field[6].DistMin = 0;
Field[6].DistMax = 0.03;

// Combine fields
Field[7] = Min;
Field[7].FieldsList = {2, 4, 6};
Background Field = 7;

// Mesh settings
Mesh.Algorithm = 6;        // Frontal-Delaunay
Mesh.RecombineAll = 1;     // Recombine triangles into quads
Mesh.RecombinationAlgorithm = 2;
Mesh.ElementOrder = 1;
Mesh.SecondOrderIncomplete = 0;

// Generate 2D mesh
Mesh 2;

// Optimize mesh quality
Mesh.Optimize = 1;
Mesh.OptimizeNetgen = 1;
