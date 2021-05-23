// Define hollow sphere with thickness and circular plug region with different
// material properties

// Mesh.Algorithm = 8;
// Mesh.RecombinationAlgorithm = 2;

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

// x-z plane
Circle(1) = {8,1,2}; Circle(2) = {2,1,3}; Circle(3) = {3,1,4};
Circle(4) = {4,1,9}; Circle(5) = {9,1,5}; Circle(6) = {5,1,8};
// yz plane
Circle(7) = {7,1,3}; Circle(8) = {3,1,6}; Circle(9) = {6,1,10};
Circle(10) = {10,1,5}; Circle(11) = {5,1,11}; Circle(12) = {11,1,7};
// xy plane
Circle(13) = {7,1,2}; Circle(14) = {2,1,6};
Circle(15) = {6,1,4}; Circle(16) = {4,1,7};
// Plug inner
Circle(17) = {11,22,8}; Circle(18) = {8,22,10};
Circle(19) = {10,22,9}; Circle(20) = {9,22,11};

// Outer x-z plane
Circle(21) = {18,1,12}; Circle(22) = {12,1,13}; Circle(23) = {13,1,14};
Circle(24) = {14,1,19}; Circle(25) = {19,1,15}; Circle(26) = {15,1,18};
// Outer yz plane
Circle(27) = {17,1,13}; Circle(28) = {13,1,16}; Circle(29) = {16,1,20};
Circle(30) = {20,1,15}; Circle(31) = {15,1,21}; Circle(32) = {21,1,17};
// Outer xy plane
Circle(33) = {17,1,12}; Circle(34) = {12,1,16};
Circle(35) = {16,1,14}; Circle(36) = {14,1,17};
// Outer Plug
Circle(37) = {21,23,18}; Circle(38) = {18,23,20};
Circle(39) = {20,23,19}; Circle(40) = {19,23,21};
// Lines
Line(41) = {2, 12}; Line(42) = {3, 13}; Line(43) = {4, 14}; Line(44) = {5, 15};
Line(45) = {6, 16}; Line(46) = {7, 17}; Line(47) = {8, 18}; Line(48) = {9, 19};
Line(49) = {10, 20}; Line(50) = {11, 21};

// Plug surface inner
Line Loop(11) = {10,6,18}; Line Loop(12) = {-10,19,5};
Line Loop(13) = {-5,20,-11}; Line Loop(14) = {11,17,-6};
// Sphere upper inner surface
Line Loop(15) = {1,14,9,-18}; Line Loop(16) = {2,8,-14};
Line Loop(17) = {3,-15,-8}; Line Loop(18) = {4,-19,-9,15};
// Sphere lower inner surface
Line Loop(19) = {1,-13,-12,17}; Line Loop(20) = {2,-7,13};
Line Loop(21) = {3,16,7}; Line Loop(22) = {4,20,12,-16};
// Plug connection between inner and outer surfaces
Line Loop(23) = {44,26,-47,-6}; Line Loop(24) = {44,-30,-49,10};
Line Loop(25) = {44,-25,-48,5}; Line Loop(26) = {44,31,-50,-11};
Line Loop(27) = {50,37,-47,-17}; Line Loop(28) = {47,38,-49,-18};
Line Loop(29) = {49,39,-48,-19}; Line Loop(30) = {48,40,-50,-20};
// Connections following the xz plane and yz plane respectively
Line Loop(31) = {47,21,-41,-1}; Line Loop(32) = {41,22,-42,-2};
Line Loop(33) = {42,23,-43,-3}; Line Loop(34) = {43,24,-48,-4};
Line Loop(35) = {46,27,-42,-7}; Line Loop(36) = {42,28,-45,-8};
Line Loop(37) = {45,29,-49,-9}; Line Loop(38) = {50,32,-46,-12};
// Plug surface outer
Line Loop(39) = {26,38,30}; Line Loop(40) = {-30,39,25};
Line Loop(41) = {-25,40,-31}; Line Loop(42) = {31,37,-26};
// Sphere upper outer surface
Line Loop(43) = {21,34,29,-38}; Line Loop(44) = {22,28,-34};
Line Loop(45) = {23,-35,-28}; Line Loop(46) = {24,-39,-29,35};
// Sphere lower outer surface
Line Loop(47) = {21,-33,-32,37}; Line Loop(48) = {22,-27,33};
Line Loop(49) = {23,36,27}; Line Loop(50) = {24,40,32,-36};
// Connections following the xy plane
Line Loop(51) = {46,33,-41,-13}; Line Loop(52) = {41,34,-45,-14};
Line Loop(53) = {45,35,-43,-15}; Line Loop(54) = {43,36,-46,-16};

// Convert line loops to surface
Surface(11) = {11}; Surface(12) = {12}; Surface(13) = {13}; Surface(14) = {14};
Surface(15) = {15}; Surface(16) = {16}; Surface(17) = {17}; Surface(18) = {18};
Surface(19) = {19}; Surface(20) = {20}; Surface(21) = {21}; Surface(22) = {22};
Surface(23) = {23}; Surface(24) = {24}; Surface(25) = {25}; Surface(26) = {26};
Surface(27) = {27}; Surface(28) = {28}; Surface(29) = {29}; Surface(30) = {30};
Surface(31) = {31}; Surface(32) = {32}; Surface(33) = {33}; Surface(34) = {34};
Surface(35) = {35}; Surface(36) = {36}; Surface(37) = {37}; Surface(38) = {38};
Surface(39) = {39}; Surface(40) = {40}; Surface(41) = {41}; Surface(42) = {42};
Surface(43) = {43}; Surface(44) = {44}; Surface(45) = {45}; Surface(46) = {46};
Surface(47) = {47}; Surface(48) = {48}; Surface(49) = {49}; Surface(50) = {50};
Surface(51) = {51}; Surface(52) = {52}; Surface(53) = {53}; Surface(54) = {54};

// Surface loops for plug
Surface Loop(101) = {11, 23, 28, 24, 39}; Surface Loop(102) = {12, 24, 29, 25, 40};
Surface Loop(103) = {13, 25, 30, 26, 41}; Surface Loop(104) = {14, 26, 27, 23, 42};
// Surface loops for top hemisphere
Surface Loop(105) = {15,31,52,37,28,43}; Surface Loop(106) = {16,32,52,36,44};
Surface Loop(107) = {17,33,36,53,45}; Surface Loop(108) = {18,34,53,37,29,46};
// Surface loops for bottom hemisphere
Surface Loop(109) = {19,31,51,38,27,47}; Surface Loop(110) = {20,32,51,35,48};
Surface Loop(111) = {21,33,35,54,49}; Surface Loop(112) = {22,34,54,38,30,50};

// Conversion to volume
Volume(1001) = {101}; Volume(1002) = {102}; Volume(1003) = {103}; Volume(1004) = {104};
Volume(1005) = {105}; Volume(1006) = {106}; Volume(1007) = {107}; Volume(1008) = {108};
Volume(1009) = {109}; Volume(1010) = {110}; Volume(1011) = {111}; Volume(1012) = {112};

Physical Surface(1) = {11,12,13,14};
Physical Surface(2) = {15,16,17,18,19,20,21,22};
Physical Surface(3) = {23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38};
Physical Surface(4) = {39,40,41,42,43,44,45,46,47,48,49,50};
Physical Surface(5) = {51,52,53,54};

Physical Volume(2001) = {1001,1002,1003,1004};
Physical Volume(2002) = {1005,1006,1007,1008,1009,1010,1011,1012};
