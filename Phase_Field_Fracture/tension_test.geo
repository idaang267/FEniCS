// Rectangle with slit in center

Mesh.Algorithm = 8; // Delaunay for quads
Mesh.RecombinationAlgorithm = 2; // or 3

// Rectangle length and height
L = 2.0;  // Length scale along x axis
H = 0.4;  // Length scale along y axis

// Mesh resolution
res = 0.009;

// DEFINE POINTS
// Points for top line
Point(1) = {0   , 0.5*H, 0, res};
Point(2) = {0.9 , 0.5*H, 0, res};
Point(3) = {1.1,  0.5*H, 0, res};
Point(4) = {L   , 0.5*H, 0, res};
// Points for upper off-center line
Point(5) = {0,   9E-3, 0, res};
Point(6) = {0.9, 9E-3, 0, res};
Point(7) = {1.1, 9E-3, 0, res};
Point(8) = {L,   9E-3, 0, res};
// Points for lower off-center line
Point(9)  = {0,   -9E-3, 0, res};
Point(10) = {0.9, -9E-3, 0, res};
Point(11) = {1.1, -9E-3, 0, res};
Point(12) = {L,   -9E-3, 0, res};
// Points along bottom line
Point(13) = {0,   -0.5*H, 0, res};
Point(14) = {0.9, -0.5*H, 0, res};
Point(15) = {1.1, -0.5*H, 0, res};
Point(16) = {L,   -0.5*H, 0, res};

// DEFINE LINES
// Top lines
Line(25) = {1,2};
Line(26) = {2,3};
Line(27) = {3,4};
// Upper off-center lines
Line(28) = {5,6};
Line(29) = {6,7};
Line(30) = {7,8};
// Lower off-center lines
Line(31) = {9,10};
Line(32) = {10,11};
Line(33) = {11,12};
// Bottom lines
Line(34) = {13,14};
Line(35) = {14,15};
Line(36) = {15,16};
// Left Lines
Line(37) = {1,5};
Line(38) = {5,9};
Line(39) = {9,13};
// Left off-center lines
Line(40) = {2,6};
Line(41) = {6,10};
Line(42) = {10,14};
// Right off-center lines
Line(43) = {3,7};
Line(44) = {7,11};
Line(45) = {11,15};
// Right line_search
Line(46) = {4,8};
Line(47) = {8,12};
Line(48) = {12,16};

// Connect lines
// Upper sections above the crack
Line Loop(50) = {25,40,-28,-37};
Line Loop(52) = {26,43,-29,-40};
Line Loop(54) = {27,46,-30,-43};
// Sections around the crack
Line Loop(56) = {28,41,-31,-38};
Line Loop(58) = {30,47,-33,-44};
// Lower sections below the crack
Line Loop(60) = {31,42,-34,-39};
Line Loop(62) = {32,45,-35,-42};
Line Loop(64) = {33,48,-36,-45};

// Surfaces
Surface(51) = {50};
Surface(53) = {52};
Surface(55) = {54};
Surface(57) = {56};
Surface(59) = {58};
Surface(61) = {60};
Surface(63) = {62};
Surface(65) = {64};

// Loop
Surface Loop(100) = {50,52,54,56,58,64,62,60};
