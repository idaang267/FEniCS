// 2D version of Shear Test

Mesh.Algorithm = 8; // Delaunay for quads

// Define points on surface
ms = 0.05;

// Dimensions in x y
x_d = 15;
y_d = 1.5;

// Points for plane defined counter-clockwise
Point(1) = {0, 0, 0, ms};
Point(2) = {x_d, 0, 0, ms};
Point(3) = {x_d, y_d, 0, ms};
Point(4) = {0, y_d, 0, ms};

// Perimeter lines
Line(1) = {1,2}; Line(2) = {2,3};
Line(3) = {3,4}; Line(4) = {4,1};

// Join lines together to form a surface
Line Loop(1) = {1,2,3,4};
// Convert line loop to surface
Plane Surface(50) = {1};
// Place crack in surface
// Line{8} In Surface {50};
// Convert surface to volume
// Volume(1000) = {50};
