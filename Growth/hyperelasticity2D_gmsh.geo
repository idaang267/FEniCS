/*********************************************************************
 *
 *  Gmsh tutorial 1
 *
 *  Variables, elementary entities (points, curves, surfaces), physical
 *  entities (points, curves, surfaces)
 *
 *********************************************************************/

// Define points on surface
ms = 0.0025;
ms2 = 0.025;
Point(1) = {0, 0, 0, ms2};
Point(2) = {1, 0, 0, ms2};
Point(3) = {0, 0.09, 0, ms};
Point(4) = {1, 0.09, 0, ms};
Point(5) = {0, 0.1, 0, ms};
Point(6) = {1, 0.1, 0, ms};

// Define straight lines
Line(1) = {1,2};
Line(2) = {2,4};
Line(3) = {1,3};
Line(4) = {3,4};
Line(5) = {4,6};
Line(6) = {3,5};

// Create sine curve
pList[0] = 5; // Start point
nPoints = 100; // Number of discretization points
For i In {1:nPoints}
  xL = 1; // Length of sine curve
  a = 0.0001;
  l = 0.00001;
  x =xL*i/(nPoints+1);
  pList[i] = newp;
  Point(pList[i]) = {x,0.1+(a*Sin(2*Pi*x/l)),0};
EndFor
pList[nPoints+1] = 6; // Last point
Spline(7) = pList[];

Transfinite Curve {4, 7} = Ceil(1/ms+1) Using Progression 1;
Transfinite Curve {6, 5} = Ceil(0.01/ms+1) Using Progression 1;


// Create final surfaces
Curve Loop(1) = {1, 2, -4, -3};
Curve Loop(2) = {4, 5, -7, -6};
Plane Surface(1) = {1};
Plane Surface(2) = {2};

Transfinite Surface {2};
Surface Loop(1) = {1,2};

