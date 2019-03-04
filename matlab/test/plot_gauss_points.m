
close all
clear all





material_points = csvread('./../../build/bin/cylinder/materialpoints.csv');



plot(material_points(:,1),material_points(:,2),'r*');

axis equal 


nodes = csvread('./../../build/bin/cylinder/nodes.csv');

hold on
plot(nodes(:,1),nodes(:,2),'b.');


tri = csvread('./../../build/bin/cylinder/triangles.csv');

triplot(tri,nodes(:,1),nodes(:,2));
