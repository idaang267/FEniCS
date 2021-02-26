
// Mesh.Algorithm = 8; // Delaunay for quads

// Define points on surface
// 0.05 0.1 and 0.3 Small Wrinkle
// 0.04 0.06 and 0.3 Wrinkling
// 0.03 0.04 and 0.1 Wrinkling
// 0.02 0.04 and 0.1 Diverges
// 0.01 0.02 and 0.1 Diverges
ms = 0.02;
ms_l = 0.04;
ms_r = 0.1;

// Dimensions in x y and z
x_d = 1.0;
y_d = 0.1;
z_d = 0.1;
// Length of crack
c_d = 0.5;

// Back plane
Point(1) = {0, 0, 0, ms_l};
Point(2) = {c_d, 0, 0, ms};
Point(3) = {x_d, 0, 0, ms_r};
Point(4) = {x_d, y_d, 0, ms_r};
Point(5) = {c_d, y_d, 0, ms};
Point(6) = {0, y_d, 0, ms_l};
Point(7) = {0, y_d/2, 0, ms_l};
Point(8) = {c_d, y_d/2, 0, ms};
// Front Plane
Point(9) = {0, 0, z_d, ms_l};
Point(10) = {c_d, 0, z_d, ms};
Point(11) = {x_d, 0, z_d, ms_r};
Point(12) = {x_d, y_d, z_d, ms_r};
Point(13) = {c_d, y_d, z_d, ms};
Point(14) = {0, y_d, z_d, ms_l};
Point(15) = {0, y_d/2, z_d, ms_l};
Point(16) = {c_d, y_d/2, z_d, ms};

// Back lines
Line(1) = {1,2}; Line(2) = {2,3};
Line(3) = {3,4}; Line(4) = {4,5};
Line(5) = {5,6}; Line(6) = {6,7};
Line(7) = {7,1};

// Front lines
Line(8) = {9,10}; Line(9) = {10,11};
Line(10) = {11,12}; Line(11) = {12,13};
Line(12) = {13,14}; Line(13) = {14,15};
Line(14) = {15,9};

// Lines making up crack edges not crack surface
Line(15) = {7,8}; Line(16) = {15,16};

// Connector lines
Line(17) = {1,9}; Line(18) = {3,11};
Line(19) = {4,12}; Line(20) = {6,14};

Line Loop(1) = {1,2,3,4,5,6,7};
Line Loop(2) = {8,9,10,11,12,13,14};
Line Loop(3) = {1,2,18,-9,-8,-17};
Line Loop(4) = {18,10,-19,-3};
Line Loop(5) = {4,5,20,-12,-11,-19};
Line Loop(6) = {6,7,17,-14,-13,-20};

Plane Surface(50) = {1};
Plane Surface(51) = {2};
Plane Surface(52) = {3};
Plane Surface(53) = {4};
Plane Surface(54) = {5};
Plane Surface(55) = {6};

Line{15} In Surface {50};
Line{16} In Surface {51};

Surface Loop(100) = {50,51,52,53,54,55};
// Convert surface to volume
Volume(1000) = {100};

/*
// Connector lines between front plane and back plane
Line(21) = {1,9}; Line(22) = {2,10};
Line(23) = {3,11}; Line(24) = {4,12};
Line(25) = {5,13}; Line(26) = {6,14};
Line(27) = {7,15}; Line(28) = {8,16};

Line Loop(1) = {1,8,9,5,6,7};
Line Loop(2) = {2,3,4,-9,-8};

Line Loop(3) = {11,18,19,15,16,17};
Line Loop(4) = {12,13,14,-19,-18};

Line Loop(5) = {1,2,23,-12,-11,-21};
Line Loop(6) = {23,13,-24,-3};
Line Loop(7) = {4,5,26,-15,-14,-24};
Line Loop(8) = {6,7,21,-17,-16,-26};

Plane Surface(50) = {1};
Plane Surface(51) = {2};
Plane Surface(52) = {3};
Plane Surface(53) = {4};
Plane Surface(54) = {5};
Plane Surface(55) = {6};
Plane Surface(56) = {7};
Plane Surface(57) = {8};

// Combine planes in loop
Surface Loop(100) = {50,51,52,53,54,55,56,57};

// Convert surface to volume
Volume(1000) = {100};
*/

/*
Physical Volume(1) = {1000};

x = 1;
p_count_f = 15;
p_count_b = 16;
inc = (x_d - c_d)/10;

For x In {1:18}
  //Printf("Printout Test %g ", p_count_f) ;
  Point(p_count_f) = {c_d+inc, y_d/2, 0, ms_l};
  Point{p_count_f} In Surface{9};
  Point(p_count_b) = {c_d+inc, y_d/2, z_d, ms_l};
  Point{p_count_b} In Surface{10};
  // Update p_count front and back for next step
  p_count_f = p_count_f + 2;
  p_count_b = p_count_b + 2;
  x = x + 1;
  inc = inc + (x_d - c_d)/10;
EndFor
*/
