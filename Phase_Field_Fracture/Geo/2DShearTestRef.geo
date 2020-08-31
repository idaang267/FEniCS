// 2D version of Shear Test

// Mesh.Algorithm = 8; // Delaunay for quads

// Define points on surface
// 0.02 0.03 and 0.09 Converges
// 0.01 0.02 and 0.08 Converges
// 0.009 0.01 and 0.03 Converges
// 0.008 0.009 and 0.02 Converges
// 0.007 0.009 and 0.02 Converges
// 0.006 0.01 and 0.02 Converges
// 0.006 0.008 and 0.02 Div
// 0.001 0.01 and 0.02 ?
ms = 0.008;
ms_l = 0.01;
ms_r = 0.02;

// Dimensions in x y
x_d = 15;
y_d = 1.5;
// Length of crack
c_d = 5;

// Points for plane defined counter-clockwise
Point(1) = {0, 0, 0, ms_l};
Point(2) = {c_d, 0, 0, ms};
Point(3) = {x_d, 0, 0, ms_r};
Point(4) = {x_d, y_d, 0, ms_r};
Point(5) = {c_d, y_d, 0, ms};
Point(6) = {0, y_d, 0, ms_l};
Point(7) = {0, y_d/2, 0, ms_l};
Point(8) = {c_d, y_d/2, 0, ms};

// Back lines
Line(1) = {1,2}; Line(2) = {2,3};
Line(3) = {3,4}; Line(4) = {4,5};
Line(5) = {5,6}; Line(6) = {6,7};
Line(7) = {7,1};
// Line making up crack edge
Line(8) = {7,8};

// Join lines together to form a surface
Line Loop(1) = {1,2,3,4,5,6,7};
// Convert line loop to surface
Plane Surface(50) = {1};
// Place crack in surface
Line{8} In Surface {50};
// Convert surface to volume
//Volume(1000) = {50};
