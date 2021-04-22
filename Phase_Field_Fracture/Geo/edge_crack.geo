
// mesh size
hr = 0.01;h  = 6*hr;
W  = 1.0;  H  = 3.0;
H0 = 10*hr;H1 = 5*hr;
W0 = H1;
T  =-2*hr;
//upper surface
Point(101) = {-W, -H,  T, h};
Point(102) = { W, -H,  T, h};
Point(103) = { W, -H0, T, hr};
Point(104) = { W,  H0, T, hr};
Point(105) = { W,  H,  T, h};
Point(106) = {-W,  H,  T, h};
Point(107) = {-W0,  H0, T, hr};
Point(108) = {-W0, -H0, T, hr};
Point(109) = {-W,  H1, T, h};
Point(110) = {-W, -H1, T, h};
Point(111) = {-W0,-H1, T, hr};
Point(112) = {-W0, H1, T, hr};
Point(113) = {-W0, 0., T, hr};

Line(101) = {101, 102};
Line(102) = {102, 103};
Line(103) = {103, 104};
Line(104) = {104, 105};
Line(105) = {105, 106};
Line(106) = {106, 109};
Line(107) = {107, 112};
Line(108) = {109, 112};
Circle(109)={111,113,112};
Line(110) = {111,110};
Line(111) = {111, 108};
Line(112) = {110, 101};
Line(113) = {108, 103};
Line(114) = {107, 104};

Line Loop(11) = {101, 102, -113, -111,  110, 112};
Line Loop(12) = {-113, -111, 109, -107, 114, -103};
Line Loop(13) = {104, 105, 106, 108, -107, 114};
Plane Surface(11) = {11};
Plane Surface(12) = {12};
Plane Surface(13) = {13};

//+
Extrude {0, 0, 0.04} {
	Surface{11};
	Surface{12};
	Surface{13};
}

// 1=MeshAdapt, 2=Automatic, 5=Delaunay, 6=Frontal, 7=BAMG, 8=DelQuad
//Mesh.Algorithm = 6;
//Mesh 3;
