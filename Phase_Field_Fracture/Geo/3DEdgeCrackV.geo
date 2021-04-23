// 3D Strip with Notch Problem

// 1=MeshAdapt, 2=Automatic, 5=Delaunay, 6=Frontal, 7=BAMG, 8=DelQuad
Mesh.Algorithm = 8;

// mesh size
hr = 0.01;	// More refined
h  = 0.2;

// Set width, height, and thickness
W  = 1.0; H  = 3.0; T =-hr;
H0 = 15*hr; H1 = 5*hr; W0 = H1;

// Upper surface Points
// Define outer points anti-clockwise
Point(1) = {-W, -H1, T, hr}; Point(2) = {-W, -H0, T, hr};
Point(3) = {-W, -H,  T, h}; Point(4) = {W, -H,  T, h};
Point(5) = { W, -H0, T, hr}; Point(6) = {W, -H1, T, hr};
Point(7) = {W, 0, T, hr}; Point(8) = {W, H1, T, hr};
Point(9) = {W,  H0, T, hr}; Point(10) = {W, H, T, h};
Point(11) = {-W,  H,  T, h}; Point(12) = {-W,  H0,  T, hr};
Point(13) = {-W,  H1, T, hr};
// Define Points for semi-circle
Point(14) = {-W0, H1, T, hr};
Point(15) = {-W0, 0., T, hr};
Point(16) = {-W0,-H1, T, hr};

// Outer lines
Line(1) = {1, 2}; Line(2) = {2, 3};
Line(3) = {3, 4}; Line(4) = {4, 5};
Line(5) = {5, 6}; Line(6) = {6, 7};
Line(7) = {7, 8}; Line(8) = {8, 9};
Line(9) = {9, 10}; Line(10) = {10, 11};
Line(11) = {11, 12}; Line(12) = {12, 13};
// Lines for notch
Line(13) = {13, 14};
Circle(14) = {16,15,14};
Line(15) = {16, 1};
// Inner lines
Line(16) = {2, 5}; Line(17) = {9, 12};
Line(18) = {16,6}; Line(19) = {8, 14};

// Bottom, Middle, and top section
Line Loop(1) = {2, 3, 4, -16};
Line Loop(2) = {1, 16, 5, -18, 15};
Line Loop(3) = {-14, 18, 6, 7, 19};
Line Loop(4) = {12, 13, -19, 8, 17};
Line Loop(5) = {11, -17, 9, 10};

// Convert to plane surface
Plane Surface(50) = {1};
Plane Surface(51) = {2};
Plane Surface(52) = {3};
Plane Surface(53) = {4};
Plane Surface(54) = {5};

Extrude {0, 0, -T} {
	Surface{50};
	Surface{51};
	Surface{52};
	Surface{53};
	Surface{54};
}
