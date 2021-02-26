// 2D version of unit square shear test

Mesh.Algorithm = 8; // Delaunay for quads

// Dimensions in x y
x_d = 0.2;
y_d = 0.2;
// Dimension related to the center point of refinement
pf = 0.02;
// Length of crack
c_d = 0.5;

// Determine N = x_d/ms_r for entry into .py code
ms_r = 0.001; // More refined
ms = 0.001;

// Points for crack - origin is center of domain
Point(1) = {-x_d/2, 0, 0, ms_r};
Point(2) = {0, 0, 0, ms_r};
// Points for plane defined counter-clockwise
Point(3) = {-x_d/2, -pf, 0, ms_r};
Point(4) = {-x_d/2, -y_d/2, 0, ms};
Point(5) = {x_d/2, -y_d/2, 0, ms};
Point(6) = {x_d/2, -pf, 0, ms_r};
Point(7) = {x_d/2, pf, 0, ms_r};
Point(8) = {x_d/2, y_d/2, 0, ms};
Point(9) = {-x_d/2, y_d/2, 0, ms};
Point(10) = {-x_d/2, pf, 0, ms_r};
// Define Points for point in surface
Point(11) = {0, y_d/2, 0, ms};
Point(12) = {0, -y_d/2, 0, ms};

// Circumference lines
Line(1) = {1,3}; Line(2) = {3,4};
Line(3) = {4,5}; Line(4) = {5,6};
Line(5) = {6,7}; Line(6) = {7,8};
Line(7) = {8,9}; Line(8) = {9,10};
Line(9) = {10,1};
// Mid lines
Line(10) = {3,6};
Line(11) = {7,10};
// Crack Line
Line(12) = {1,2};

// Join lines together to form a surface
Line Loop(1) = {1,10,5,11,9};
Line Loop(2) = {2,3,4,-10};
Line Loop(3) = {6,7,8,-11};
// Convert line loop to surface
Plane Surface(50) = {1};
Plane Surface(51) = {2};
Plane Surface(52) = {3};
// Place crack in surface
Line{12} In Surface {50};
// Combine planes in loop
Surface Loop(100) = {50,51,52};
//
Point{11} In Surface{51};
Point{12} In Surface{52};

// Convert surface to volume
//Volume(1000) = {50};
