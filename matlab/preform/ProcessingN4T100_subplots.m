clear all
close all

PLOT_GRAPHS = true;
TMAX = 0.3;


path = './../../build/bin/preform/Displacement/';
pathSR = './../../build/bin/preform/srRod/';
boundaryNodes = csvread('./../../build/bin/preform/boundary.txt');
boundaryNodes = [boundaryNodes;boundaryNodes(1)];

addpath(path)
displacementdir = path ;
d = dir(displacementdir);
d1 = dir([displacementdir,'*.csv']);
numFiles = size(d,1) -3 ;

plotFiles = ceil(linspace(1,numFiles,10));


figure
for i = 1:10
   
    subplot(1,10,i)
    
    
plot_point =171;
filename = strcat(path,'displacement_',num2str(plotFiles(i)),'.csv');
disp = csvread(filename,1);

plot(disp(:,1),disp(:,2),'k.','markersize',5)           % line plot
hold on
plot(-disp(:,1),disp(:,2),'k.','markersize',5)           % line plot
ymax = max(disp(:,2));
hold on
plot(disp(plot_point,1),disp(plot_point,2),'r*')
axis equal 
hold on
plot(disp(boundaryNodes,1),disp(boundaryNodes,2),'b--')
hold on
plot(-disp(boundaryNodes,1),disp(boundaryNodes,2),'b--')

hold on
filename = strcat(pathSR,'srRod_',num2str(plotFiles(i)),'.csv');
disp = csvread(filename,1);
hold on
plot(disp(:,1),disp(:,2),'r-','linewidth',1)           % line plot
hold on
plot(-disp(:,1),disp(:,2),'r-','linewidth',1)           % line plot
xlim([-50,50])
ylim([-170,ymax])

    
end


print -dpng2 Displacement.png
