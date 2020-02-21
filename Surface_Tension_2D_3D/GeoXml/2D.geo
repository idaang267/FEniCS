// Gmsh 3D Sphere

// Unit sphere
RSphere = 1.0;

// Mesh resolution
lcSphere = 0.1; //.25;

// Define Points of circle
Point(1) = {0,0,0,lcSphere};
Point(2) = {RSphere,0,0,lcSphere};
Point(3) = {0,RSphere,0,lcSphere};
Point(4) = {-RSphere,0,0.0,lcSphere};
Point(5) = {0,-RSphere,0.0,lcSphere};
// Circle command creates an arc, subdivide circle into four arcs
Circle(1) = {2,1,3};
Circle(2) = {3,1,4};
Circle(3) = {4,1,5};
Circle(4) = {5,1,2};
// Define all lines pointing outwards from center
Line(5) = {1,3};
Line(6) = {1,5};
Line(7) = {1,2};
Line(8) = {1,4};
// Connect all lines
Line Loop(7) = {1,-5,7};
Line Loop(8) = {5,2,-8};
Line Loop(9) = {8,3,-6};
Line Loop(10) = {6,4,-7};
// Create surfaces from each line loop
Surface(11) = {7};
Surface(12) = {8};
Surface(13) = {9};
Surface(14) = {10};
