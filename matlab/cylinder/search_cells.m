
close all
clear all





material_points = csvread('./../../build/bin/cylinder/materialpoints.csv');



plot(material_points(:,1),material_points(:,2),'r*');

axis equal 


nodes = csvread('./../../build/bin/cylinder/nodes.csv');

hold on
plot(nodes(:,1),nodes(:,2),'b.');


%tri = csvread('./../../build/bin/cylinder/triangles.csv');

%triplot(tri,nodes(:,1),nodes(:,2));


cells = csvread('./../../build/bin/cylinder/search_cells.csv');


active_cells =  csvread('./../../build/bin/cylinder/search_cells_nodes.csv');


for k = 1:length(active_cells)
hold on
    i = active_cells(k)+1;
rectangle('Position',[cells(i,1), cells(i,3), cells(i,2)-cells(i,1), cells(i,4)-cells(i,3)])
end


