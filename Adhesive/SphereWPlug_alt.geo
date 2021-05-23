// Define hollow sphere with thickness and circular plug region with different
// material properties

Mesh.Algorithm = 8;
Mesh.RecombinationAlgorithm = 2;

res = 0.4;  // Mesh Resolution
R = 7.1; // Sphere Radius
Rp = 2.5; // Plug Radius
t = 0.2;  // Wall thickness
// Arc length equation
th = Rp/R;
h = R*(1-Cos(th));
c = (R^2 - (R-h)^2)^(1/2);
Rt = R + t;
ht = Rt*(1-Cos(th));
ct = (Rt^2 - (Rt-ht)^2)^(1/2);
// Determining center point of circle on circle
beta = 0; gamma = Pi/2; t = 0;
z = Sin(th)*Sin(beta)*Cos(t)+ Cos(th)*Cos(beta);

// plane points inner
Point(1) = {0, 0, 0, res};
Point(2) = {R, 0, 0, res}; Point(3) = {0, 0, -R, res};
Point(4) = {-R, 0, 0, res}; Point(5) = {0, 0, R, res};
Point(6) = {0, R, 0, res}; Point(7) = {0, -R, 0, res};
// Define points for plug
Point(8) = {c, 0, R-h, res}; Point(9) = {-c, 0, R-h, res};
Point(10) = {0, c, R-h, res}; Point(11) = {0, -c, R-h, res};
// plane points outer
Point(12) = {Rt, 0, 0, res}; Point(13) = {0, 0, -Rt, res};
Point(14) = {-Rt, 0, 0, res}; Point(15) = {0, 0, Rt, res};
Point(16) = {0, Rt, 0, res}; Point(17) = {0, -Rt, 0, res};
// Define points for plug
Point(18) = {ct, 0, Rt-ht, res}; Point(19) = {-ct, 0, Rt-ht, res};
Point(20) = {0, ct, Rt-ht, res}; Point(21) = {0, -ct, Rt-ht, res};
// Center points for plug
Point(22) = {0, 0, R*z}; Point(23) = {0, 0, Rt*z};

// Inner xz plane
Circle(1) = {5,1,2}; Circle(2) = {2,1,3};
Circle(3) = {3,1,4}; Circle(4) = {4,1,5};
// Inner yz plane
Circle(5) = {7,1,3}; Circle(6) = {3,1,6};
Circle(7) = {6,1,5}; Circle(8) = {5,1,7};
// Inner xy plane
Circle(9) = {7,1,2}; Circle(10) = {2,1,6};
Circle(11) = {6,1,4}; Circle(12) = {4,1,7};
// Outer xz plane
Circle(13) = {15,1,12}; Circle(14) = {12,1,13};
Circle(15) = {13,1,14}; Circle(16) = {14,1,15};
// Outer yz plane
Circle(17) = {17,1,13}; Circle(18) = {13,1,16};
Circle(19) = {16,1,15}; Circle(20) = {15,1,17};
// Outer xy plane
Circle(21) = {17,1,12}; Circle(22) = {12,1,16};
Circle(23) = {16,1,14}; Circle(24) = {14,1,17};
// Connector Lines
Line(25) = {2, 12}; Line(26) = {3, 13};
Line(27) = {4, 14}; Line(28) = {5, 15};
Line(29) = {6, 16}; Line(30) = {7, 17};

// Sphere upper inner surface
Line Loop(11) = {1,10,7}; Line Loop(12) = {2,6,-10};
Line Loop(13) = {3,-11,-6}; Line Loop(14) = {4,-7,11};
// Sphere lower inner surface
Line Loop(15) = {1,-9,-8}; Line Loop(16) = {2,-5,9};
Line Loop(17) = {3,12,5}; Line Loop(18) = {4,8,-12};
// Sphere upper outer surface
Line Loop(19) = {13,22,19}; Line Loop(20) = {14,18,-22};
Line Loop(21) = {15,-23,-18}; Line Loop(22) = {16,-19,23};
// Sphere lower outer surface
Line Loop(23) = {13,-21,-20}; Line Loop(24) = {14,-17,21};
Line Loop(25) = {15,24,17}; Line Loop(26) = {16,20,-24};
// Connections following the xz plane
Line Loop(30) = {28,13,-25,-1}; Line Loop(31) = {25,14,-26,-2};
Line Loop(32) = {26,15,-27,-3}; Line Loop(33) = {27,16,-28,-4};
// Connections following the yz plane
Line Loop(34) = {30,17,-26,-5}; Line Loop(35) = {26,18,-29,-6};
Line Loop(36) = {29,19,-28,-7}; Line Loop(37) = {28,20,-30,-8};
// Connections following the xy plane
Line Loop(38) = {30,21,-25,-9}; Line Loop(39) = {25,22,-29,-10};
Line Loop(40) = {29,23,-27,-11}; Line Loop(41) = {27,24,-30,-12};

// Convert line loops to surface
Surface(11) = {11}; Surface(12) = {12}; Surface(13) = {13}; Surface(14) = {14};
Surface(15) = {15}; Surface(16) = {16}; Surface(17) = {17}; Surface(18) = {18};
Surface(19) = {19}; Surface(20) = {20}; Surface(21) = {21}; Surface(22) = {22};
Surface(23) = {23}; Surface(24) = {24}; Surface(25) = {25}; Surface(26) = {26};
Surface(30) = {30}; Surface(31) = {31}; Surface(32) = {32}; Surface(33) = {33};
Surface(34) = {34}; Surface(35) = {35}; Surface(36) = {36}; Surface(37) = {37};
Surface(38) = {38}; Surface(39) = {39}; Surface(40) = {40}; Surface(41) = {41};

// Connector Lines within surfaces
Line(31) = {8,18}; Line(32) = {9,19};
Line(33) = {10,20}; Line(34) = {11,21};
// Inner Plug Lines to add into surface
Circle(35) = {11,22,8}; Circle(36) = {8,22,10};
Circle(37) = {10,22,9}; Circle(38) = {9,22,11};
// Outer Plug lines to add into surface
Circle(39) = {21,23,18}; Circle(40) = {18,23,20};
Circle(41) = {20,23,19}; Circle(42) = {19,23,21};

// Inside Points
Point{11} In Surface{15}; Point{8} In Surface{15};
Point{8} In Surface{11}; Point{10} In Surface{11};
Point{10} In Surface{14}; Point{9} In Surface{14};
Point{9} In Surface{18}; Point{11} In Surface{18};
// Outside Points
Point{21} In Surface{23}; Point{18} In Surface{23};
Point{18} In Surface{19}; Point{20} In Surface{19};
Point{20} In Surface{22}; Point{19} In Surface{22};
Point{19} In Surface{26}; Point{21} In Surface{26};

//
Line Loop(105) = {35,31,-39,-34};
Surface(105) = {105};

// Bottom volume of sphere
Surface Loop(101) = {15, 37, 30, 38, 23};
Volume(1001) = {101};

Surface{105} In Volume{1001};
