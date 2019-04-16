
close all
clear all





material_points = csvread('./../../build/bin/preform_alt/materialpoints.csv');
figure
cells = csvread('./../../build/bin/preform_alt/search_cells.csv');

for i = 1:length(cells)
hold on

rectangle('Position',[cells(i,1), cells(i,3), cells(i,2)-cells(i,1), cells(i,4)-cells(i,3)])
end

active_cells =  csvread('./../../build/bin/preform_alt/active_cells.csv');



for k = 1:length(active_cells)
hold on
    i = active_cells(k)+1;
rectangle('Position',[cells(i,1), cells(i,3), cells(i,2)-cells(i,1), cells(i,4)-cells(i,3)],'FaceColor',[0 .5 .5])
end

axis equal

figure

plot(material_points(:,1),material_points(:,2),'r*');

%axis equal 


nodes = csvread('./../../build/bin/preform_alt/nodes.csv');

hold on
plot(nodes(:,1),nodes(:,2),'b.');


tri = csvread('./../../build/bin/preform_alt/triangles.csv');

triplot(tri,nodes(:,1),nodes(:,2));

axis equal



cells = csvread('./../../build/bin/preform_alt/search_cells.csv');
active_cells =  csvread('./../../build/bin/preform_alt/active_cells.csv');



for k = 1:length(active_cells)
hold on
    i = active_cells(k)+1;
rectangle('Position',[cells(i,1), cells(i,3), cells(i,2)-cells(i,1), cells(i,4)-cells(i,3)])
end