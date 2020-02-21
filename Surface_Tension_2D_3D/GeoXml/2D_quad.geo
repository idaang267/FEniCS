// Gmsh 3D Sphere

Mesh.Algorithm = 8; // Delaunay for quads
// Mesh.RecombinationAlgorithm = 2; // or 3

// Unit sphere
RSphere = 1.0;
RCoord = 0.707106781188;

// Mesh resolution
lcSphere = 0.1;

// DEFINE POINTS
// Define inside points making up inner sphere
Point(1) = {0,0,0,lcSphere};
Point(2) = {RCoord,RCoord,0,lcSphere};
Point(3) = {-RCoord,RCoord,0,lcSphere};
Point(4) = {-RCoord,-RCoord,0,lcSphere};
Point(5) = {RCoord,-RCoord,0,lcSphere};

Point(6) = {0,RSphere,0,lcSphere};
Point(7) = {-RSphere,0,0,lcSphere};
Point(8) = {0,-RSphere,0,lcSphere};
Point(9) = {RSphere,0,0,lcSphere};

// Define arcs of interest
// xy-plane of central sphere
Circle(1) = {2,1,6};
Circle(2) = {6,1,3};
Circle(3) = {3,1,7};
Circle(4) = {7,1,4};
Circle(5) = {4,1,8};
Circle(6) = {8,1,5};
Circle(7) = {5,1,9};
Circle(8) = {9,1,2};

// Create box lines
Line(9) = {2,3};
Line(10) = {3,4};
Line(11) = {4,5};
Line(12) = {5,2};

// Create crossing lines
Line(13) = {1,2};
Line(14) = {1,3};
Line(15) = {1,4};
Line(16) = {1,5};

// Connect lines in two loops, a quarter circle each
Line Loop(50) = {13,9,-14};
Line Loop(52) = {14,10,-15};
Line Loop(54) = {15,11,-16};
Line Loop(56) = {16,12,-13};

// Connect edges of squares
Line Loop(58) = {1,2,-9};
Line Loop(60) = {3,4,-10};
Line Loop(62) = {5,6,-11};
Line Loop(64) = {7,8,-12};

// Surfaces
Surface(51) = {50};
Surface(53) = {52};
Surface(55) = {54};
Surface(57) = {56};
Surface(59) = {58};
Surface(61) = {60};
Surface(63) = {62};
Surface(65) = {64};

//
Surface Loop(100) = {51,53,55,57,59,61,63,65};

//
Recombine Surface{100};
// Recombine Surface{101};

Physical Surface(1) = {100};
//Volume(1000) = {100};
//Volume(2001) = {101};


//
//Recombine Volume{2000};
//Recombine Volume{2001};


// Physical Volume(1) = {2000,2001,2002,2003};
