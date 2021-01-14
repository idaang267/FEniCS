// Half Rectangle with left hand crack

Mesh.Algorithm = 8; // Delaunay for quads
Mesh.RecombinationAlgorithm = 2; // or 3

// Rectangle length and height
L = 5.0;  // Length scale along x axis
H = 1.0;  // Length scale along y axis
C = L/2;  // Length scale of crack

// Mesh resolution
res_r = 0.0000001;
res = 0.1;

// DEFINE POINTS
Point(1) = {0   , 0, 0, res_r};
Point(2) = {L/2 , 0, 0, res};
Point(3) = {L/2 , H, 0, res};
Point(4) = {0   , H, 0, res};
Point(5) = {-L/2, H, 0, res};
Point(6) = {-L/2, 0, 0, res};

// DEFINE LINES

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,5};
Line(5) = {5,6};
Line(6) = {6,1};
Line(7) = {4,1};

// Connect lines
Line Loop(50) = {1,2,3,7};
Line Loop(51) = {4,5,6,-7};
// Surfaces
Plane Surface(100) = {50};
Plane Surface(101) = {51};
