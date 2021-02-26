// Embedded Crack for discrete crack simulations

Mesh.Algorithm = 8; // Delaunay for quads

// Dimensions in x y
x_d = 0.01;
y_d = 0.01;
// Length of crack
c_d = 0.0005;
ms_r = 0.0001;
ms = 0.0001;

// Points
Point(1) = {-x_d/2, -0.000001, 0, ms_r};
Point(2) = {-x_d/2, -y_d/2, 0, ms_r};
Point(3) = {0, -y_d/2, 0, ms_r};
Point(4) = {x_d/2, -y_d/2, 0, ms};
Point(5) = {x_d/2, 0, 0, ms};
Point(6) = {x_d/2, y_d/2, 0, ms};
Point(7) = {0, y_d/2, 0, ms_r};
Point(8) = {-x_d/2, y_d/2, 0, ms_r};
Point(9) = {-x_d/2, +0.000001, 0, ms_r};
Point(10) = {0, 0, 0, ms_r};

// Lines
Line(1) ={1,2}; Line(2) ={2,3};
Line(3) ={3,4}; Line(4) ={4,5};
Line(5) ={5,6}; Line(6) ={6,7};
Line(7) ={7,8}; Line(8) ={8,9};
// Center lines
Line(9) ={9,10}; Line(10) = {1,10};
Line(11) = {3,10};
Line(12) = {5,10};
Line(13) = {7,10};

//Make line loops and surfaces
Line Loop(1) = {1,2,11,-10};
Line Loop(2) = {3,4,12,-11};
Line Loop(3) = {5,6,13,-12};
Line Loop(4) = {7,8,9,-13};

Plane Surface(1) = {1};
Plane Surface(2) = {2};
Plane Surface(3) = {3};
Plane Surface(4) = {4};

//Assign physical IDs (numbers picked just so they will be obvious in .msh file)
// Physical Surface(100) = {1}; // bottom will have ID 100
// Physical Surface(101) = {2}; // top will have ID 101
// Physical Line(10) = {7}; // crack will have ID 10

// more physical ids on top and bottom for boundary conditions in a solver
// Physical Line(1) = {1};
// 101 Physical Line(2) = {4};
