// Gmsh project created on Tue Apr 16 13:42:28 2019
//+
Point(1) = {0, 20, 0, 0};
//+
Point(2) = {0, 21.57, 0, 0};
//+
Point(3) = {20, 0, 0, 0};
//+
Point(4) = {21.57, 0, 0, 0};
//+
Point(5) = {0, 0, 0, 0};
//+
Circle(1) = {1, 5, 3};
//+
Circle(2) = {2, 5, 4};
//+
Line(3) = {1, 2};
//+
Line(4) = {4, 3};
//+
Line Loop(1) = {2, 4, -1, 3};
//+
Plane Surface(1) = {1};
//+
Physical Line("PRESSURE") = {1};
//+
Physical Line("DISPLACEMENT_1") = {4};
//+
Physical Line("DISPLACEMENT_2") = {3};
//+
Physical Surface("MATERIAL") = {1};
