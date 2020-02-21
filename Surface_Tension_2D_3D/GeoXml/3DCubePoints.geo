// Benchmark existing on GMSH tutorials

//Delete All;
Mesh.Algorithm = 8;
Mesh.Algorithm3D = 1;

RCube = 1.;
lcCube = 0.08;

// Define all points of the circle
Point(1) = {0,     0,     0,     lcCube};
Point(2) = {RCube, 0,     0,     lcCube};
Point(3) = {RCube, 0,     RCube, lcCube};
Point(4) = {0,     0,     RCube, lcCube};
Point(5) = {0,     RCube, 0,     lcCube};
Point(6) = {RCube, RCube, 0,     lcCube};
Point(7) = {RCube, RCube, RCube, lcCube};
Point(8) = {0,     RCube, RCube, lcCube};

// Bottom Lines
Line(10) = {1, 2};
Line(11) = {2, 3};
Line(12) = {3, 4};
Line(13) = {4, 1};
// Four Lines
Line(14) = {1, 5};
Line(15) = {2, 6};
Line(16) = {3, 7};
Line(17) = {4, 8};
// Top Lines
Line(18) = {5, 6};
Line(19) = {6, 7};
Line(20) = {7, 8};
Line(21) = {8, 5};

Line Loop(10) = {10, 11, 12, 13};
Surface(11) = {10};
Line Loop(12) = {10, 15, -18, -14};
Surface(13) = {12};
Line Loop(14) = {11, 16, -19, -15};
Surface(15) = {14};
Line Loop(16) = {12, 17, -20, -16};
Surface(17) = {16};
Line Loop(18) = {13, 14, -21, -17};
Surface(19) = {18};
Line Loop(20) = {18, 19, 20, 21};
Surface(21) = {20};

Surface Loop(29) = {11,13,15,17,19,21};
Volume(1111) = {29};

p0 = newp;
Point(p0) = {0.5, RCube, 0.5, lcCube};
Point{p0} In Volume {1111};

p1 = newp;
Point(p1) = {0.6, RCube, 0.5, lcCube};
Point{p1} In Volume {1111};
l1 = newl;
Line(l1) = {p0, p1};
Curve{l1} In Volume {1111};

p2 = newp;
Point(p2) = {0.7, RCube, 0.5, lcCube};
Point{p2} In Volume {1111};
l2 = newl;
Line(l2) = {p1, p2};
Curve{l2} In Volume {1111};

p3 = newp;
Point(p3) = {0.8, RCube, 0.5, lcCube};
Point{p3} In Volume {1111};
l3 = newl;
Line(l3) = {p2, p3};
Curve{l3} In Volume {1111};

p4 = newp;
Point(p4) = {0.9, RCube, 0.5, lcCube};
Point{p4} In Volume {1111};
l4 = newl;
Line(l4) = {p3, p4};
Curve{l4} In Volume {1111};

p5 = newp;
Point(p5) = {1.0, RCube, 0.5, lcCube};
Point{p5} In Volume {1111};
l5 = newl;
Line(l5) = {p4, p5};
Curve{l5} In Volume {1111};
