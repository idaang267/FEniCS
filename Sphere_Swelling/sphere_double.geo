// Gmsh 3D Sphere

Mesh.Algorithm = 8; // Delaunay for quads
Mesh.RecombinationAlgorithm = 2; // or 3

// Unit sphere
RSphere = 0.9;
RSphere2 = 1.0;

// Mesh resolution
lcSphere = 0.1; //.25;
lcSphere2 = 0.1;

// Define inside points making up inner sphere
Point(1) = {0,0,0,lcSphere};
Point(2) = {RSphere,0,0,lcSphere};
Point(3) = {0,RSphere,0,lcSphere};
Point(4) = {-RSphere,0,0.0,lcSphere};
Point(5) = {0,-RSphere,0.0,lcSphere};
Point(6) = {0,0,-RSphere,lcSphere};
Point(7) = {0,0,RSphere,lcSphere};
// Define outside points making up outer sphere
Point(8)  = {RSphere2,0,0,lcSphere2};
Point(9)  = {0,RSphere2,0,lcSphere2};
Point(10) = {-RSphere2,0,0.0,lcSphere2};
Point(11) = {0,-RSphere2,0.0,lcSphere2};
Point(12) = {0,0,-RSphere2,lcSphere2};
Point(13) = {0,0,RSphere2,lcSphere2};

// Define arcs of interest
// xy-plane of central sphere
Circle(1) = {2,1,3};
Circle(2) = {3,1,4};
Circle(3) = {4,1,5};
Circle(4) = {5,1,2};
// yz-plane of central sphere
Circle(5) = {3,1,6};
Circle(6) = {6,1,5};
Circle(7) = {5,1,7};
Circle(8) = {7,1,3};
// xz-plane of central sphere
Circle(9) = {2,1,7};
Circle(10) = {7,1,4};
Circle(11) = {4,1,6};
Circle(12) = {6,1,2};
// xy-plane of outer sphere
Circle(13) = {8,1,9};
Circle(14) = {9,1,10};
Circle(15) = {10,1,11};
Circle(16) = {11,1,8};
// yz-plane of outer sphere
Circle(17) = {9,1,12};
Circle(18) = {12,1,11};
Circle(19) = {11,1,13};
Circle(20) = {13,1,9};
// xz-plane of outer sphere
Circle(21) = {8,1,13};
Circle(22) = {13,1,10};
Circle(23) = {10,1,12};
Circle(24) = {12,1,8};

// Create six extra lines between the inner and outer sphere
Line(25) = {2,8};
Line(26) = {3,9};
Line(27) = {4,10};
Line(28) = {5,11};
Line(29) = {6,12};
Line(30) = {7,13};

// Create two extra lines on x-axis for centerpoint visualization
Line(31) = {1,4};
Line(32) = {1,2};

// Connect all lines in a loop
// Central sphere
Line Loop(25) = {2,8,-10};
Line Loop(27) = {10,3,7};

Line Loop(29) = {-8,-9,1};
Line Loop(37) = {-7,4,9};

Line Loop(31) = {-11,-2,5};
Line Loop(39) = {-4,12,-6};

Line Loop(33) = {-5,-12,-1};
Line Loop(35) = {-3,11,6};

// Outer Sphere
Line Loop(41) = {14,20,-22};
Line Loop(43) = {22,15,19};

Line Loop(45) = {-20,-21,13};
Line Loop(47) = {-23,-14,17};

Line Loop(49) = {-17,-24,-13};
Line Loop(51) = {-15,23,18};

Line Loop(53) = {-19,16,21};
Line Loop(55) = {-16,24,-18};

// Flat planes in xy
Line Loop(59) = {1,26,-13,-25};
Line Loop(61) = {2,27,-14,-26};
Line Loop(63) = {3,28,-15,-27};
Line Loop(65) = {4,25,-16,-28};

// Inner sphere plane in XY
Line Loop(67) = {1,2,-31,32};
Line Loop(69) = {31,3,4,-32};

// Surfaces on hemisphere in front of xy plane for inner circle
Surface(26) = {25};
Surface(28) = {27};
Surface(30) = {29};
Surface(38) = {37};
// Surface on hemisphere behind xy plane for inner circle
Surface(32) = {31};
Surface(34) = {33};
Surface(36) = {35};
Surface(40) = {39};

// Surfaces on hemisphere in front of xy plane for outer circle
Surface(42) = {41};
Surface(44) = {43};
Surface(46) = {45};
Surface(54) = {53};
// Surface on hemisphere behind xy plane for outer circle
Surface(48) = {47};
Surface(50) = {49};
Surface(52) = {51};
Surface(56) = {55};

//
Surface(60) = {59};
Surface(62) = {61};
Surface(64) = {63};
Surface(66) = {65};

// Inner sphere plane in XY
Surface(68) = {67};
Surface(70) = {69};

//
Surface Loop(100) = {68,70,26,28,38,30};
Surface Loop(101) = {68,70,32,34,36,40};
Surface Loop(102) = {26,28,30,38,60,62,64,66,42,44,46,54};
Surface Loop(103) = {32,34,36,40,60,62,64,66,48,50,52,56};

//
Recombine Surface{100};
Recombine Surface{101};
Recombine Surface{102};
Recombine Surface{103};

Volume(2000) = {100};
Volume(2001) = {101};
Volume(2002) = {102};
Volume(2003) = {103};

//
Recombine Volume{2000};
Recombine Volume{2001};
Recombine Volume{2002};
Recombine Volume{2003};

Physical Volume(1) = {2000,2001,2002,2003};

//
// Surface Loop(100) = {40,38,28,26,32,36,34,30};
