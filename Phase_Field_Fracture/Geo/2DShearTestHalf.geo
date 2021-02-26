// Half Rectangle with left hand crack

Mesh.Algorithm = 8; // Delaunay for quads
Mesh.RecombinationAlgorithm = 2; // or 3

// Rectangle length and height
L = 5.0;  // Length scale along x axis
H = 1.0;  // Length scale along y axis
pf = 0.01; // Length scale along y axis of refinement

// Mesh resolution
res_r = 0.001; // 1E-7;
res = 0.1; // 0.1;

// DEFINE POINTS
Point(1) = {0   , 0,  0, res_r};
Point(2) = {L/2 , 0,  0, res_r};
Point(3) = {L/2 , pf,  0, res_r};
Point(4) = {L/2 , H,  0, res};
Point(5) = {0   , H,  0, res};
Point(6) = {-L/2, H,  0, res};
Point(7) = {-L/2, pf,  0, res_r};
Point(8) = {-L/2, 0,  0, res_r};

// DEFINE LINES
Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,5};
Line(5) = {5,6};
Line(6) = {6,7};
Line(7) = {7,8};
Line(8) = {8,1};
Line(9) = {3,7};

// Connect lines
Line Loop(1) = {1,2,9,7,8};
Line Loop(2) = {3,4,5,6,-9};

// Surfaces
Plane Surface(50) = {1};
Plane Surface(51) = {2};

Surface Loop(100) = {50,51};
