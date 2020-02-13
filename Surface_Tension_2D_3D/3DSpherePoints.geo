// Benchmark existing on GMSH tutorials

//Delete All;
Mesh.Algorithm = 2;
Mesh.Algorithm3D = 1;

RSphere = 1.;
lcSphere = 0.15;

// Define all points of the circle
Point(1) = {0,         0,       0, lcSphere};
Point(2) = {RSphere,   0,       0, lcSphere};
Point(3) = {0,         RSphere, 0, lcSphere};
Point(4) = {-RSphere,  0,       0, lcSphere};
Point(5) = {0,        -RSphere, 0, lcSphere};
Point(6) = {0,         0,       -RSphere,lcSphere};
Point(7) = {0,         0,       RSphere,lcSphere};

Circle(1) = {2,1,3};
Circle(2) = {3,1,4};
Circle(3) = {4,1,5};
Circle(4) = {5,1,2};
Circle(5) = {3,1,6};
Circle(6) = {6,1,5};
Circle(7) = {5,1,7};
Circle(8) = {7,1,3};
Circle(9) = {2,1,7};
Circle(10) = {7,1,4};
Circle(11) = {4,1,6};
Circle(12) = {6,1,2};

Line Loop(13) = {2,8,-10};
Surface(14) = {13};
Line Loop(15) = {10,3,7};
Surface(16) = {15};
Line Loop(17) = {-8,-9,1};
Surface(18) = {17};
Line Loop(19) = {-11,-2,5};
Surface(20) = {19};
Line Loop(21) = {-5,-12,-1};
Surface(22) = {21};
Line Loop(23) = {-3,11,6};
Surface(24) = {23};
Line Loop(25) = {-7,4,9};
Surface(26) = {25};
Line Loop(27) = {-4,12,-6};
Surface(28) = {27};

Surface Loop(29) = {28,26,16,14,20,24,22,18};
Volume(1111) = {29};

p0 = newp;
Point(p0) = {0.1, 0.0, 0.0, lcSphere};
Point{p0} In Volume {1111};
l0 = newl;
Line(l0) = {1, p0};
Curve{l0} In Volume {1111};

p1 = newp;
Point(p1) = {0.2, 0.0, 0.0, lcSphere};
Point{p1} In Volume {1111};
l1 = newl;
Line(l1) = {p0, p1};
Curve{l1} In Volume {1111};

p2 = newp;
Point(p2) = {0.3, 0.0, 0.0, lcSphere};
Point{p2} In Volume {1111};
l2 = newl;
Line(l2) = {p1, p2};
Curve{l2} In Volume {1111};

p3 = newp;
Point(p3) = {0.4, 0.0, 0.0, lcSphere};
Point{p3} In Volume {1111};
l3 = newl;
Line(l3) = {p2, p3};
Curve{l3} In Volume {1111};

p4 = newp;
Point(p4) = {0.5, 0.0, 0.0, lcSphere};
Point{p4} In Volume {1111};
l4 = newl;
Line(l4) = {p3, p4};
Curve{l4} In Volume {1111};

p5 = newp;
Point(p5) = {0.6, 0.0, 0.0, lcSphere};
Point{p5} In Volume {1111};
l5 = newl;
Line(l5) = {p4, p5};
Curve{l5} In Volume {1111};

p6 = newp;
Point(p6) = {0.7, 0.0, 0.0, lcSphere};
Point{p6} In Volume {1111};
l6 = newl;
Line(l6) = {p5, p6};
Curve{l6} In Volume {1111};

p7 = newp;
Point(p7) = {0.8, 0.0, 0.0, lcSphere};
Point{p7} In Volume {1111};
l7 = newl;
Line(l7) = {p6, p7};
Curve{l7} In Volume {1111};

p8 = newp;
Point(p8) = {0.9, 0.0, 0.0, lcSphere};
Point{p8} In Volume {1111};
l8 = newl;
Line(l8) = {p7, p8};
Curve{l8} In Volume {1111};

//pLoc = 0.0;
//subD = RSphere/lcSphere;
//For pLoc In {0:RSphere}
//  p[pLoc] = newp;
//  l[pLoc] = newl;
//  Point(p[pLoc]) = {pLoc, 0.0, 0.0, lcSphere};
//  Point{p[pLoc]} In Volume {1111};
//  Line(l[pLoc]) = {1, p[pLoc]};
//  Curve{l[pLoc]} In Volume {1111};
//  pLoc = pLoc + subD;

// Wall ID = 2100
//Physical Surface(2100) = {28,26,16,14,20,24,22,18};
// Farfield CBC ID = 2222
//Physical Surface(2222) = {128,126,116,114,120,124,122,118};
// Interior ID = 1000
//Physical Volume(1000) = {200};
