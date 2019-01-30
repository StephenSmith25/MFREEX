clc;
clear all;
close all;



% read cells and sites




[vorCells,nodes] = readCells("simple");

figure
for i = 1:length(vorCells)
    rand1 = rand(1,3);
    tempArray = vorCells{i};
    hold on
    fill(tempArray(:,1),tempArray(:,2),rand1);
end

nodes = csvread("simple.sites");
boundaryNodes = csvread("boundary.txt");
boundaryNodes = [boundaryNodes ; boundaryNodes(1)];
plot(nodes(boundaryNodes,1),nodes(boundaryNodes,2),'r-','linewidth',2)
hold on
hold on
plot(nodes(:,1),nodes(:,2),'k.','markersize',6)
axis off
axis equal
print(gcf,'cells.png','-dpng','-r600');         