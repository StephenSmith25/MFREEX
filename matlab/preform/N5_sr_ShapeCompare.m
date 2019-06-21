
close all
clear all 

load ./Experimental/20190515_forStephen.mat

a = process_data(4);
b = a.outer_strain;
c = a.middle_strain;
time = a.time;
pressure = a.cavity_pressure;
dt = a.duration/size(b,3);

outer_coord = a.outer_Vic3D_coord;

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







axis equal 
hold on
plot(disp(boundaryNodes,1),disp(boundaryNodes,2),'k-')


outer_coords_i = outer_coord(:,:,1);
hold on 
plot(outer_coords_i(:,1),outer_coords_i(:,2),'g','linewidth',2);


% -------------------------------------------------------------------------%
%                           PLOT 2
% -------------------------------------------------------------------------%




filename = strcat(path,'displacement_',num2str(plotFiles(4)),'.csv');
disp = csvread(filename,1);




plot(disp(:,1),disp(:,2),'k.','markersize',3)           % line plot
hold on
plot(disp(boundaryNodes,1),disp(boundaryNodes,2),'k-')




outer_coords_i = outer_coord(:,:,floor(10*end/24));
hold on 
plot(outer_coords_i(:,1),outer_coords_i(:,2),'r','linewidth',2);
hold on 

% -------------------------------------------------------------------------%
%                           PLOT 3
% -------------------------------------------------------------------------%
filename = strcat(path,'displacement_',num2str(plotFiles(10)),'.csv');
disp = csvread(filename,1);




plot(disp(:,1),disp(:,2),'k.','markersize',5)           % line plot
hold on
plot(disp(boundaryNodes,1),disp(boundaryNodes,2),'k-')



outer_coords_i = outer_coord(:,:,floor(26*end/32));
hold on 
plot(outer_coords_i(:,1),outer_coords_i(:,2),'b','linewidth',2);



% 
axis equal
xlim([0,60])
ylim([-350,ymax])



%set axis
set(gca, 'FontName', 'cmr12')
% set x tics and y tics
set(gca,...
'Units','normalized',...
'FontWeight','normal',...
'FontSize',14,... % size ofiguref numbers on axis
'FontName','cmr14') % font name
set(gca,'TickLabelInterpreter', 'latex');
set(gcf, 'Color', 'w');


% Y label
ylabel({'z (mm)'},...
'interpreter','latex',...
'FontSize',14,... % font size
'FontName','cmr14')
% X label
xlabel('r (mm)',...
'interpreter','latex',...
'FontSize',14,... % font size
'FontName','cmr14')
% title



export_fig n5shape.png -m5


