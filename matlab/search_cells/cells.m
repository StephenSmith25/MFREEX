

close all
clear all


search_cells= csvread('./../../build/bin/cell_search/search_cells.csv');



figure
for i = 1:length(search_cells)
hold on
rectangle('Position',[search_cells(i,1), search_cells(i,3), search_cells(i,2)-search_cells(i,1),....
    search_cells(i,4)-search_cells(i,3)])
end

hold on

r = 0.5;

theta = linspace(0,360,50)';


nodes = csvread('./../../build/bin/cell_search/nodes.csv');
hold on
plot(nodes(:,1),nodes(:,2),'b.');

hold on
plot(5,5,'r.');


X = [5+r.*cosd(theta),5+r.*sind(theta)];

hold on
plot(X(:,1),X(:,2),'r.');





xlim([0 10])
ylim([0 10])

