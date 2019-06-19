
close all
clear all 

load ./Experimental/20190515_forStephen.mat

a = process_data(5);
b = a.outer_strain;
c = a.middle_strain;
time = a.time;
pressure = a.cavity_pressure;
dt = a.duration/size(b,3);

outer_coord = a.outer_coord;

t = 0;
T_OFFSET = 0.1005; %0.1005;

element = 84;





path = './../../build/bin/preform/Displacement/';
pathSR = './../../build/bin/preform/srRod/';
addpath(path)


displacementdir = path ;
d = dir(displacementdir);
d1 = dir([displacementdir,'*.csv']);
numFiles = size(d,1) -3 ;

plotFiles = ceil(linspace(1,numFiles,10));

plot_point =138; %%2241
filename = strcat(path,'displacement_',num2str(plotFiles(1)),'.csv');
disp = csvread(filename,1);


boundaryNodes = csvread('./../../build/bin/preform/boundary.txt');

boundaryNodes = [boundaryNodes;boundaryNodes(1)];



% -------------------------------------------------------------------------%
%                           PLOT 1
% -------------------------------------------------------------------------%




plot(disp(:,1),disp(:,2),'k.','markersize',3)           % line plot
hold on
ymax = 0;
hold on



axis equal 
hold on
plot(disp(boundaryNodes,1),disp(boundaryNodes,2),'k-')




outer_coords_i = outer_coord(:,:,1);
hold on 
plot(outer_coords_i(:,1),outer_coords_i(:,2),'b','linewidth',2);
hold on 


% -------------------------------------------------------------------------%
%                           PLOT 2
% -------------------------------------------------------------------------%




filename = strcat(path,'displacement_',num2str(plotFiles(5)),'.csv');
disp = csvread(filename,1);




plot(disp(:,1),disp(:,2),'k.','markersize',3)           % line plot
hold on
plot(disp(boundaryNodes,1),disp(boundaryNodes,2),'k-')


outer_coords_i = outer_coord(:,:,floor(18*end/32));
hold on 
plot(outer_coords_i(:,1),outer_coords_i(:,2),'b','linewidth',2);


% -------------------------------------------------------------------------%
%                           PLOT 3
% -------------------------------------------------------------------------%
filename = strcat(path,'displacement_',num2str(plotFiles(10)),'.csv');
disp = csvread(filename,1);




plot(disp(:,1),disp(:,2),'k.','markersize',5)           % line plot
hold on
plot(disp(boundaryNodes,1),disp(boundaryNodes,2),'k-')

hold on




outer_coords_i = outer_coord(:,:,floor(32*end/32));
hold on 
plot(outer_coords_i(:,1),outer_coords_i(:,2),'b','linewidth',2);
hold on 



axis equal
axis off
xlim([0,45])
ylim([-350,ymax])


