// Embedded Crack for discrete crack simulations

Mesh.Algorithm = 8; // Delaunay for quads

// Dimension in r
r_d = 0.01;
// Length of crack
c_d = 0.01;
// rougher mesh to finer
ms_r = 5E-4;
ms = 1E-7;

// Points
Point(1) = {-r_d, 0.000001, 0, ms_r};
Point(2) = {-r_d, -0.000001, 0, ms_r};
Point(3) = {0, -r_d, 0, ms_r};
Point(4) = {r_d, 0, 0, ms_r};
Point(5) = {0, r_d, 0, ms_r};
Point(6) = {0, 0, 0, ms}; // < centerpoint

// Arcs
Circle(1) = {2,6,3};
Circle(2) = {3,6,4};
Circle(3) = {4,6,5};
Circle(4) = {5,6,1};

// Lines
Line(5) ={1,6};
Line(6) ={6,2};

//Make line loops and surfaces
Line Loop(1) = {1,2,3,4,5,6};
Plane Surface(1) = {1};
