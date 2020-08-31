// 2D version of Shear Test
// Determine N = x_d/ms_r for entry into .py code

//Mesh.Algorithm = 8; // Delaunay for quads

// Define points on surface
// 0.05 and 0.1 Div
// 0.04 and 0.08 Div
// 0.01 and 0.05 Con
// 0.008 and 0.05 Con
ms_r = 0.01; // More refined
ms = 0.05;

// Dimensions in x y
x_d = 15;
y_d = 1.5;
pf = 0.2;
// Length of crack
c_d = 5;

// Points for plane defined counter-clockwise
Point(1) = {0, 0, 0, ms};
Point(2) = {x_d, 0, 0, ms};
Point(3) = {x_d, y_d/2-pf, 0, ms_r};
Point(4) = {x_d, y_d/2+pf, 0, ms_r};
Point(5) = {x_d, y_d, 0, ms};
Point(6) = {0, y_d, 0, ms};
Point(7) = {0, y_d/2+pf, 0, ms_r};
Point(8) = {0, y_d/2-pf, 0, ms_r};
// Points for crack
Point(9) = {0, y_d/2, 0, ms_r};
Point(10) = {c_d, y_d/2, 0, ms_r};

// Circumference lines
Line(1) = {1,2}; Line(2) = {2,3};
Line(3) = {3,4}; Line(4) = {4,5};
Line(5) = {5,6}; Line(6) = {6,7};
Line(7) = {7,8}; Line(8) = {8,1};
// Mid lines
Line(9) = {8,3};
Line(10) = {4,7};
// Line making up crack edge
Line(11) = {9,10};

// Join lines together to form a surface
Line Loop(1) = {1,2,-9,8};
Line Loop(2) = {9,3,10,7};
Line Loop(3) = {4,5,6,-10};
// Convert line loop to surface
Plane Surface(50) = {1};
Plane Surface(51) = {2};
Plane Surface(52) = {3};
// Place crack in surface
Line{11} In Surface {51};
// Combine planes in loop
Surface Loop(100) = {50,51,52};
// Convert surface to volume
//Volume(1000) = {50};
