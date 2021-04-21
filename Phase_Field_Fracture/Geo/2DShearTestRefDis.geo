// Embedded Crack for discrete crack simulations
// Origin is at (0,0,0)

Mesh.Algorithm = 8; // Delaunay for quads

// Dimensions in x y
x_d = 6.0; x_m = 1.0;
y_d = 1.0; y_m = 0.25;
// Setting mesh resolutions
ms_r = 0.01; // Mesh refinement
ms = 0.05;

// Bottom Outer points
Point(1) = {-x_d/2, -0.000001, 0, ms_r};
Point(2) = {-x_d/2, -y_m/2, 0, ms_r};
Point(3) = {-x_d/2, -y_d/2, 0, ms};
Point(4) = {0, -y_d/2, 0, ms};
Point(5) = {x_d/2, -y_d/2, 0, ms};
Point(6) = {x_d/2, -y_m/2, 0, ms};
// Top Outer points
Point(7) = {x_d/2, y_m/2, 0, ms};
Point(8) = {x_d/2, y_d/2, 0, ms};
Point(9) = {0, y_d/2, 0, ms};
Point(10) = {-x_d/2, y_d/2, 0, ms};
Point(11) = {-x_d/2, y_m/2, 0, ms_r};
Point(12) = {-x_d/2, +0.000001, 0, ms_r};
// Inner Points
Point(13) = {0, 0, 0, ms_r};
Point(14) = {0, -y_m/2, 0, ms_r};
Point(15) = {x_m, -y_m/2, 0, ms_r};
Point(16) = {x_m, y_m/2, 0, ms_r};
Point(17) = {0, y_m/2, 0, ms_r};

// Outer Lines
Line(1) ={1,2}; Line(2) ={2,3};
Line(3) ={3,4}; Line(4) ={4,5};
Line(5) ={5,6}; Line(6) ={6,7};
Line(7) ={7,8}; Line(8) ={8,9};
Line(9) ={9,10}; Line(10) = {10,11};
Line(11) = {11,12};
// Lines around crack
Line(12) = {12,13}; Line(13) = {13,1};
// Inner Circumference Lines
Line(14) = {2,14}; Line(15) = {14,15};
Line(16) = {15,6}; Line(17) = {7,16};
Line(18) = {16,17}; Line(19) = {17,11};
// Innermost lines
Line(20) = {13,14};
Line(21) = {15,16};
Line(22) = {17,13};

//Make line loops and surfaces
Line Loop(1) = {1,14,-20,13};
Line Loop(2) = {11,12,-22,19};
Line Loop(3) = {20,15,21,18,22};
Line Loop(4) = {-21,16,6,17};
// Outer line loops
Line Loop(5) = {2,3,4,5,-16,-15,-14};
Line Loop(6) = {10,-19,-18,-17,7,8,9};

// Convert line loop to surface
Plane Surface(50) = {1};
Plane Surface(51) = {2};
Plane Surface(52) = {3};
Plane Surface(53) = {4};
Plane Surface(54) = {5};
Plane Surface(55) = {6};
