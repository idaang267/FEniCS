// Gmsh project created on Wed Nov  6 15:07:48 2019
resolution2 = 0.1;
resolution1 = 0.1;
l0 = 1;
///////////////////// Points ///////////////////////
// Outer surface
Point(1)  = {0.00/l0, 0.00/l0, 0.00/l0, resolution1};
Point(2)  = {1.00/l0, 0.00/l0, 0.00/l0, resolution1};
Point(3)  = {1.00/l0, 1.00/l0, 0.00/l0, resolution1};
Point(4)  = {0.00/l0, 1.00/l0, 0.00/l0, resolution1};

// Cylinder
Point(10) = {0.50/l0, 0.50/l0, 0.00/l0, resolution2};
// Inner
Point(11) = {0.75/l0, 0.50/l0, 0.00/l0, resolution2};
Point(12) = {0.50/l0, 0.75/l0, 0.00/l0, resolution2};
Point(13) = {0.25/l0, 0.50/l0, 0.00/l0, resolution2};
Point(14) = {0.50/l0, 0.25/l0, 0.00/l0, resolution2};
// Outer
Point(15) = {0.80/l0, 0.50/l0, 0.00/l0, resolution2};
Point(16) = {0.50/l0, 0.80/l0, 0.00/l0, resolution2};
Point(17) = {0.20/l0, 0.50/l0, 0.00/l0, resolution2};
Point(18) = {0.50/l0, 0.20/l0, 0.00/l0, resolution2};

//////////////////////// Line ///////////////////////
// Tissue
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Curve Loop(5) = {1,2,3,4};

// Cylinders
Circle(10) = {11, 10, 12};
Circle(11) = {12, 10, 13};
Circle(12) = {13, 10, 14};
Circle(13) = {14, 10, 11};
Curve Loop(14) = {10,11,12,13};
Circle(15) = {15, 10, 16};
Circle(16) = {16, 10, 17};
Circle(17) = {17, 10, 18};
Circle(18) = {18, 10, 15};
Curve Loop(19) = {15,16,17,18};

////////////////////// Surface ///////////////////////
// Cylinder
Plane Surface(1) = {14, 19};
// Tissue
Plane Surface(2) = {5, 19};

////////////////////// Volume ///////////////////////
Extrude {0, 0, 1/l0} {
  Surface{1}; Surface{2}; 
}
Surface Loop(1) = {103, 74, 2, 78, 82, 86, 61, 32, 1, 36, 40, 44};
//Volume(3) = {1};

////////////////// Physical domain ///////////////////
Physical Volume(99) = {1, 2};
Physical Surface(1) = {61};                        // top of cylinder
Physical Surface(2) = {32, 44, 36, 40};            // side of inner cylinder
//Physical Surface(3) = {103};                       // Top surface except for cylinder
Physical Surface(3) = {103, 74, 78, 82, 86, 1, 2}; // the others of 1 and 2
//Physical Surface(3) = {2, 1};                      // bottom surface of cube
//Physical Surface(3) = {74, 78, 82, 86};            // four side surfaces of cube
//Physical Surface(6) = {56, 60, 48, 52};            // side of outer cylinder

