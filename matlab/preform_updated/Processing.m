
clear all
close all

PLOT_GRAPHS = false;
PLOT_DOMAINS_INFLUENCE = true; 


WITH_MOULD = true;
TMAX = 0.6;





path = './../../build/bin/preform_alt/Displacement/';
pathSR = './../../build/bin/preform_alt/srRod/';
boundaryNodes = csvread('./../../build/bin/preform_alt/boundary.txt');

boundaryNodes = [boundaryNodes;boundaryNodes(1)];

addpath(path)
displacementdir = path ;
d = dir(displacementdir);
d1 = dir([displacementdir,'*.csv']);
numFiles = size(d,1) -3 ;

plotFiles = ceil(linspace(1,numFiles,10));

plot_point =177;
filename = strcat(path,'displacement_',num2str(plotFiles(1)),'.txt');
disp = csvread(filename);

% draw mould
Radius_mould = 24;
Height_mould = 30;
thickness_mould = 5;
mould_length = 180;
height = 74.2650;
Radius_out = 11.25;
top_point = 48.920;



mould_nodes = [];
mould_nodes = [mould_nodes ; [Radius_out,height+5]   ];
mould_nodes = [mould_nodes ; [Radius_out,68.92]   ];
mould_nodes = [mould_nodes ; [Radius_out+Radius_mould,top_point]   ];
mould_nodes = [mould_nodes ; [Radius_out+Radius_mould,height-mould_length]   ];
mould_nodes = [mould_nodes ; [0,height-mould_length]   ];


mould_nodes = [mould_nodes ; [0,height-mould_length-thickness_mould]   ];
mould_nodes = [mould_nodes ; [Radius_out+Radius_mould+thickness_mould,height-mould_length-thickness_mould]   ];
mould_nodes = [mould_nodes ; [Radius_out+Radius_mould+thickness_mould,top_point]   ];
mould_nodes = [mould_nodes ; [Radius_out+thickness_mould,68.92]   ];
mould_nodes = [mould_nodes ; [Radius_out+thickness_mould,height+5]   ];
mould_nodes = [mould_nodes ; [Radius_out,height+5]   ];

if ( WITH_MOULD == false)
    mould_nodes = zeros(1,2);
end




figure




subplot(1,3,1)       % add first plot in 2 x 2 grid

hold on 
fill(mould_nodes(:,1),mould_nodes(:,2),'w');
hold on
fill(-mould_nodes(:,1),mould_nodes(:,2),'w');
hold on
plot(disp(:,1),disp(:,2),'k.','markersize',3)           % line plot
hold on
plot(-disp(:,1),disp(:,2),'k.','markersize',3)           % line plot
ymax = max(disp(:,2));



axis equal 
hold on
%plot(disp(boundaryNodes,1),disp(boundaryNodes,2),'b-')
hold on
%plot(-disp(boundaryNodes,1),disp(boundaryNodes,2),'b-')


boundary_nodes_xy = disp(boundaryNodes,1:2);
Rout = max(disp(:,1));
height = max(disp(:,2));




hold on
filename = strcat(pathSR,'srRod_',num2str(plotFiles(1)),'.csv');
disp = csvread(filename,1);
hold on
plot(disp(:,1),disp(:,2),'r-','linewidth',1)           % line plot
hold on
plot(-disp(:,1),disp(:,2),'r-','linewidth',1)           % line plot
xlim([-50,50])
ylim([-170,ymax])

       % line plot




filename = strcat(path,'displacement_',num2str(plotFiles(6)),'.txt');
disp = csvread(filename);
subplot(1,3,2)       % add first plot in 2 x 2 grid
hold on 
fill(mould_nodes(:,1),mould_nodes(:,2),'w');
hold on
fill(-mould_nodes(:,1),mould_nodes(:,2),'w');
hold on



plot(disp(:,1),disp(:,2),'k.','markersize',3)           % line plot
hold on
plot(-disp(:,1),disp(:,2),'k.','markersize',3) % line plot
hold on
%plot(disp(plot_point,1),disp(plot_point,2),'r*')
axis equal
hold on
%plot(disp(boundaryNodes,1),disp(boundaryNodes,2),'b-')
hold on
%plot(-disp(boundaryNodes,1),disp(boundaryNodes,2),'b-')
%plot(disp(boundaryNodes,1),disp(boundaryNodes,2),'b-')
hold on
filename = strcat(pathSR,'srRod_',num2str(plotFiles(6)),'.csv');
disp = csvread(filename,1);
hold on
plot(disp(:,1),disp(:,2),'r-','linewidth',1)           % line plot
hold on
plot(-disp(:,1),disp(:,2),'r-','linewidth',1)           % line plot
xlim([-50,50])
ylim([-170,ymax])

filename = strcat(path,'displacement_',num2str(plotFiles(10)),'.txt');
disp = csvread(filename);


subplot(1,3,3)       % add first plot in 2 x 2 grid

hold on 
fill(mould_nodes(:,1),mould_nodes(:,2),'w');
hold on
fill(-mould_nodes(:,1),mould_nodes(:,2),'w');
hold on

plot(disp(:,1),disp(:,2),'k.','markersize',5)           % line plot
hold on
plot(-disp(:,1),disp(:,2),'k.','markersize',5)           % line plot
hold on
%plot(disp(boundaryNodes,1),disp(boundaryNodes,2),'b-')
hold on
%plot(-disp(boundaryNodes,1),disp(boundaryNodes,2),'b-')
hold on
%plot(disp(plot_point,1),disp(plot_point,2),'r*')


hold on
plot(disp(181,1),disp(181,2),'ro');
hold on
plot(disp(180,1),disp(180,2),'ro');
hold on
plot(disp(2,1),disp(2,2),'ro');
hold on
plot(disp(182,1),disp(182,2),'ro');

hold on
plot(disp(267,1),disp(267,2),'ro');
hold on
plot(disp(352,1),disp(352,2),'ro');
hold on
plot(disp(176,1),disp(176,2),'ro');
hold on
plot(0.724474,-2.774268,'m*');

axis equal
hold on
filename = strcat(pathSR,'srRod_',num2str(plotFiles(10)),'.csv');
disp = csvread(filename,1);
hold on
plot(disp(:,1),disp(:,2),'r-','linewidth',1)           % line plot
hold on
plot(-disp(:,1),disp(:,2),'r-','linewidth',1)           % line plot
xlim([-50,50])
ylim([-170,ymax])




print -dpng2 Displacement.png


figure

cells = csvread('./../../build/bin/preform_alt/search_cells.csv');
active_cells =  csvread('./../../build/bin/preform_alt/active_cells.csv');


filename = strcat(path,'displacement_',num2str(plotFiles(10)),'.txt');
disp = csvread(filename);


plot(disp(:,1),disp(:,2),'k.','markersize',5)           % line plot

hold on
plot(disp(266,1),disp(266,2),'yo');
hold on
plot(disp(265,1),disp(265,2),'bo');
hold on
plot(disp(366,1),disp(366,2),'go');
hold on
plot(disp(365,1),disp(365,2),'ko');

hold on
plot(disp(465,1),disp(465,2),'cyo');
hold on
plot(disp(153,1),disp(153,2),'ro');
hold on
%plot(disp(176,1),disp(176,2),'ro');
hold on
plot(17.166,30.97,'m*');

cells = csvread('./../../build/bin/preform_alt/search_cells.csv');
active_cells =  csvread('./../../build/bin/preform_alt/active_cells.csv');



for k = 1:length(active_cells)
hold on
    i = active_cells(k)+1;
rectangle('Position',[cells(i,1), cells(i,3), cells(i,2)-cells(i,1), cells(i,4)-cells(i,3)])
end

axis equal

m = dlmread('./../../build/bin/preform_alt/pressureTime.txt',' ');
%m(ceil(length(m)/4):end,:) = smoothdata(m(ceil(length(m)/4):end,:))
figure
plot(m(:,1),m(:,2),'k-','linewidth',3);
xlabel('time');
ylabel('Pressure');
xlim([0,TMAX])
ylim([0,1])

b = csvread("../preform/Experimental/N4T100_exp_pressure.csv");
hold on
plot(b(:,1),b(:,2)/1000,'b--')


c = csvread("../preform/Experimental/N4T100_fe_pressure.csv");

hold on
plot(c(:,1),c(:,2)/1000,'r--')

legend('Mesh-free','Experiemntal','FE (Abaqus)')
ylabel('Pressure (MPa)')
xlabel('Time (s)')


%set axis
set(gca, 'FontName', 'cmr12')
% set x tics and y tics
set(gca,...
'Units','normalized',...
'FontWeight','normal',...
'FontSize',14,... % size ofiguref numbers on axis
'FontName','cmr14') % font name
set(gca,'TickLabelInterpreter', 'latex');

% Y label
ylabel({'Pressure (MPa) '},...
'interpreter','latex',...
'FontSize',14,... % font size
'FontName','cmr14')
% X label
xlabel('Time (s)',...
'interpreter','latex',...
'FontWeight','normal',...
'FontSize',14,... % font size
'FontName','cmr14')
% legend

legend({'Meshfree','Experimental','FE (Abaqus )'},... % { 'legend1', 'legend2',...}
'interpreter','latex',...
'FontSize',12,...
'FontName','cmr14',...
'Location','SouthEast')


set(gcf, 'Color', 'w');

print -dpng2 pressure.png



set(gcf, 'Color', 'w');

