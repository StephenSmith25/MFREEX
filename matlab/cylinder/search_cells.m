
close all
clear all





%material_points = csvread('./../../build/bin/cylinder/materialpoints.csv');



%plot(material_points(:,1),material_points(:,2),'r*');

axis equal 


nodes = csvread('./../../build/bin/cylinder/nodes.csv');

hold on
plot(nodes(:,1),nodes(:,2),'b.');


%tri = csvread('./../../build/bin/cylinder/triangles.csv');

%triplot(tri,nodes(:,1),nodes(:,2));


cells = csvread('./../../build/bin/cylinder/search_cells.csv');


cell_nodes =  csvread('./../../build/bin/cylinder/search_cells_nodes.csv');

ix = find(cell_nodes(:,1) == 1);

for k = 1:length(ix)
hold on
    i = ix(k)
rectangle('Position',[cells(i,1), cells(i,3), cells(i,2)-cells(i,1), cells(i,4)-cells(i,3)])
end


