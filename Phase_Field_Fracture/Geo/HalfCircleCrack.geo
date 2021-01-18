// Embedded Crack for discrete crack simulations

Mesh.Algorithm = 8; // Delaunay for quads

// Dimension in r
r_d = 0.1;
// Length of crack
c_d = 0.05;
// rougher mesh to finer
ms_r = 1E-5; // 1E-6;
ms = 1E-2; // 2E-3;

// Points
Point(1) = {0, 0, 0, ms_r}; // < centerpoint
Point(2) = {r_d, 0, 0, ms};
Point(3) = {0, r_d, 0, ms};
Point(4) = {-r_d, 0, 0, ms};

// Arcs
Circle(1) = {2,1,3};
Circle(2) = {3,1,4};

// Lines
Line(5) ={4,1};
Line(6) ={1,2};

//Make line loops and surfaces
Line Loop(1) = {1,2,5,6};
Plane Surface(1) = {1};
