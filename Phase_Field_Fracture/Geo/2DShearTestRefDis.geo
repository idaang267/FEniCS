// Embedded Crack for discrete crack simulations
// Origin is at (0,0,0)

Mesh.Algorithm = 8; // Delaunay for quads

// Dimensions in x y
x_d = 6.0; x_m = 1.0;
y_d = 1.0; y_m = 0.15;
// Setting mesh resolutions
ms_r = 0.005; // Mesh refinement
ms = 0.15;

// Bottom Outer points
Point(1) = {-x_d/2, -0.000001, 0, ms_r};
Point(2) = {-x_d/2, -y_m/2, 0, ms_r};
Point(3) = {-x_d/2, -y_d/2, 0, ms};
Point(4) = {0, -y_d/2, 0, ms};
Point(5) = {x_d/2, -y_d/2, 0, ms};
Point(6) = {x_d/2, -y_m/2, 0, ms_r};
// Mid point for BC
Point(7) = {x_d/2, 0, 0, ms_r};
// Top Outer points
Point(8) = {x_d/2, y_m/2, 0, ms_r};
Point(9) = {x_d/2, y_d/2, 0, ms};
Point(10) = {0, y_d/2, 0, ms};
Point(11) = {-x_d/2, y_d/2, 0, ms};
Point(12) = {-x_d/2, y_m/2, 0, ms_r};
Point(13) = {-x_d/2, +0.000001, 0, ms_r};
// Inner Points
Point(14) = {0, 0, 0, ms_r};
Point(15) = {0, -y_m/2, 0, ms_r};
Point(16) = {x_m, -y_m/2, 0, ms_r};
Point(17) = {x_m, 0, 0, ms_r};
Point(18) = {x_m, y_m/2, 0, ms_r};
Point(19) = {0, y_m/2, 0, ms_r};

// Outer Lines
Line(1) ={1,2}; Line(2) ={2,3};
Line(3) ={3,4}; Line(4) ={4,5};
Line(5) ={5,6}; Line(6) ={6,7};
Line(7) ={7,8}; Line(8) ={8,9};
Line(9) ={9,10}; Line(10) = {10,11};
Line(11) = {11,12}; Line(12) = {12,13};
// Lines around crack
Line(13) = {13,14}; Line(14) = {14,1};
// Inner Circumference Lines
Line(15) = {2,15}; Line(16) = {15,16};
Line(17) = {16,6}; Line(18) = {8,18};
Line(19) = {18,19}; Line(20) = {19,12};
// Innermost lines
Line(21) = {14,15}; Line(22) = {16,17};
Line(23) = {17,18}; Line(24) = {19,14};
Line(25) = {14,17}; Line(26) = {17,7};

//Make line loops and surfaces
Line Loop(1) = {12,13,-24,20};
Line Loop(2) = {1,15,-21,14};
Line Loop(3) = {24,25,23,19};
Line Loop(4) = {21,16,22,-25};
Line Loop(5) = {-23,26,7,18};
Line Loop(6) = {-22,17,6,-26};
// Outer line loops
Line Loop(7) = {2,3,4,5,-17,-16,-15};
Line Loop(8) = {11,-20,-19,-18,8,9,10};

// Convert line loop to surface
Plane Surface(50) = {1};
Plane Surface(51) = {2};
Plane Surface(52) = {3};
Plane Surface(53) = {4};
Plane Surface(54) = {5};
Plane Surface(55) = {6};
Plane Surface(56) = {7};
Plane Surface(57) = {8};
