Mesh.Algorithm = 8; // Delaunay for quads

// Define points on surface
ms = 0.002;
ms_l = 0.002;
ms_r = 0.02;

// Dimensions in x y and z
x_d = 1.0;
y_d = 0.1;
z_d = 0.1;
// Length of crack
c_d = 0.5;

// Back plane
Point(1) = {0, 0, 0, ms};
Point(2) = {c_d, 0, 0, ms};
Point(3) = {x_d, 0, 0, ms_r};
Point(4) = {x_d, y_d, 0, ms_r};
Point(5) = {c_d, y_d, 0, ms};
Point(6) = {0, y_d, 0, ms};
Point(7) = {0, y_d/2, 0, ms};
// Inner point
Point(8) = {c_d, y_d/2, 0, ms};
// Front plane
Point(9) = {0, 0, z_d, ms};
Point(10) = {c_d, 0, z_d, ms};
Point(11) = {x_d, 0, z_d, ms_r};
Point(12) = {x_d, y_d, z_d, ms_r};
Point(13) = {c_d, y_d, z_d, ms};
Point(14) = {0, y_d, z_d, ms};
Point(15) = {0, y_d/2, z_d, ms};
// Inner point
Point(16) = {c_d, y_d/2, z_d, ms};

// Back lines
Line(1) = {1,2}; Line(2) = {2,3};
Line(3) = {3,4}; Line(4) = {4,5};
Line(5) = {5,6}; Line(6) = {6,7};
Line(7) = {7,1}; Line(8) = {2,8};
Line(9) = {8,5}; Line(10) = {8,7};
// Front Lines
Line(11) = {9,10}; Line(12) = {10,11};
Line(13) = {11,12}; Line(14) = {12,13};
Line(15) = {13,14}; Line(16) = {14,15};
Line(17) = {15,9}; Line(18) = {10,16};
Line(19) = {16,13}; Line(20) = {16,15};
// Connection lines
Line(21) = {1,9}; Line(22) = {2,10};
Line(23) = {3,11}; Line(24) = {4,12};
Line(25) = {5,13}; Line(26) = {6,14};
Line(27) = {7,15}; Line(28) = {8,16};

// Curve Loops front and back
Curve Loop(1) = {1,8,10,7};
Curve Loop(2) = {9,5,6,-10};
Curve Loop(3) = {2,3,4,-9,-8};
Curve Loop(4) = {11,18,20,17};
Curve Loop(5) = {19,15,16,-20};
Curve Loop(6) = {12,13,14,-19,-18};

// Curve Loops connecting front to back
Curve Loop(7) = {21,11,-22,-1};  // Bottom left
Curve Loop(8) = {22,12,-23,-2};  // Bottom right
Curve Loop(9) = {-3,23,13,-24}; // Right
Curve Loop(10) = {4,25,-14,-24}; // Top right
Curve Loop(11) = {5,26,-15,-25}; // Top left
Curve Loop(12) = {6,27,-16,-26}; // Left Top
Curve Loop(13) = {7,21,-17,-27}; // Left Bottom
// Curve Loop(10) = {10,27,-20,-28}; // Crack surface

Plane Surface(1) = {1};
Plane Surface(2) = {2};
Plane Surface(3) = {3};
Plane Surface(4) = {4};
Plane Surface(5) = {5};
Plane Surface(6) = {6};
Plane Surface(7) = {7};
Plane Surface(8) = {8};
Plane Surface(9) = {9};
Plane Surface(10) = {10};
Plane Surface(11) = {11};
Plane Surface(12) = {12};
Plane Surface(13) = {13};

// Combine planes
Surface Loop(100) = {1,2,3,4,5,6,7,8,9,10,11,12,13};

// Convert surface to volume
Volume(1000) = {100};

/*
// Transfinite curve
Transfinite Surface{10} = {1};
Transfinite Surface{11} = {2};
Transfinite Surface{12} = {3};
Transfinite Surface{13} = {4};
Transfinite Surface{14} = {5};
Transfinite Surface{15} = {6};

Recombine Surface{10};
Recombine Surface{11};
Recombine Surface{12};
Recombine Surface{13};
Recombine Surface{14};
Recombine Surface{15};


// Crack lines in surface
//Point{17} In Surface{1};
//Point{18} In Surface{2};

*/
