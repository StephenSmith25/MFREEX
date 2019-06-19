
clear all
close all

PLOT_GRAPHS = true;
PLOT_DOMAINS_INFLUENCE = true; 


WITH_MOULD = false;
TMAX = 0.30;




path = './../../build/bin/preform/Displacement/';
pathSR = './../../build/bin/preform/srRod/';
boundaryNodes = csvread('./../../build/bin/preform/boundary.txt');
path_base = './../../build/bin/preform';

boundaryNodes = [boundaryNodes;boundaryNodes(1)];

addpath(path)
displacementdir = path ;
d = dir(displacementdir);
d1 = dir([displacementdir,'*.csv']);
numFiles = size(d,1) -3 ;

plotFiles = ceil(linspace(1,numFiles,10));

plot_point =155; %%2241
filename = strcat(path,'displacement_',num2str(plotFiles(1)),'.csv');
disp = csvread(filename,1);

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



% -------------------------------------------------------------------------%
%                           PLOT 1
% -------------------------------------------------------------------------%



subplot(1,3,1)       % add first plot in 2 x 2 grid

hold on 
fill(mould_nodes(:,1),mould_nodes(:,2),'w');
hold on
fill(-mould_nodes(:,1),mould_nodes(:,2),'w');
hold on
plot(disp(:,1),disp(:,2),'k.','markersize',3)           % line plot
hold on
plot(-disp(:,1),disp(:,2),'k.','markersize',3)           % line plot
ymax = 0;
hold on
plot(-disp(plot_point,1),disp(plot_point,2),'r*')
hold on
plot(disp(plot_point,1),disp(plot_point,2),'r*')




c = [[disp(:,1);-disp(:,1)],[disp(:,2);disp(:,2)]]








save('disp_tstart.dat', 'c', '-ascii', '-double', '-tabs')



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
plot(disp(:,1),disp(:,2),'r');
hold on
plot(-disp(:,1),disp(:,2),'r');


c = [[disp(1:end-1,1);-disp(end:-1:1,1)],[disp(1:end-1,2);disp(end:-1:1,2)]]


save('disp_rod_tstart.dat', 'c', '-ascii', '-double', '-tabs')





hold on
xlim([-50,50])
ylim([-350,ymax])



% -------------------------------------------------------------------------%
%                           PLOT 2
% -------------------------------------------------------------------------%




filename = strcat(path,'displacement_',num2str(plotFiles(5)),'.csv');
disp = csvread(filename,1);
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
hold on
plot(-disp(plot_point,1),disp(plot_point,2),'r*')
hold on
plot(disp(plot_point,1),disp(plot_point,2),'r*')


c = [[disp(:,1);-disp(:,1)],[disp(:,2);disp(:,2)]]
save('disp_tmid.dat', 'c', '-ascii', '-double', '-tabs')




%plot(disp(plot_point,1),disp(plot_point,2),'r*')
axis equal
hold on
%plot(disp(boundaryNodes,1),disp(boundaryNodes,2),'b-')
hold on
%plot(-disp(boundaryNodes,1),disp(boundaryNodes,2),'b-')
%plot(disp(boundaryNodes,1),disp(boundaryNodes,2),'b-')
hold on

hold on
filename = strcat(pathSR,'srRod_',num2str(plotFiles(5)),'.csv');
disp = csvread(filename,1);
hold on
plot(disp(:,1),disp(:,2),'r');
hold on
plot(-disp(:,1),disp(:,2),'r');


c = [[disp(1:end-1,1);-disp(end:-1:1,1)],[disp(1:end-1,2);disp(end:-1:1,2)]]


save('disp_rod_tmid.dat', 'c', '-ascii', '-double', '-tabs')





xlim([-50,50])
ylim([-350,ymax])



print -dpng2 Displacement.png

% -------------------------------------------------------------------------%
%                           PLOT 3
% -------------------------------------------------------------------------%
filename = strcat(path,'displacement_',num2str(plotFiles(10)),'.csv');
disp = csvread(filename,1);




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
plot(-disp(plot_point,1),disp(plot_point,2),'r*')
hold on
plot(disp(plot_point,1),disp(plot_point,2),'r*')


c = [[disp(:,1);-disp(:,1)],[disp(:,2);disp(:,2)]];
save('disp_tend.dat', 'c', '-ascii', '-double', '-tabs')



axis equal
hold on


hold on
filename = strcat(pathSR,'srRod_',num2str(plotFiles(10)),'.csv');
disp = csvread(filename,1);
hold on
plot(disp(:,1),disp(:,2),'r');
hold on
plot(-disp(:,1),disp(:,2),'r');


c = [[disp(1:end-1,1);-disp(end:-1:1,1)],[disp(1:end-1,2);disp(end:-1:1,2)]]


save('disp_rod_tend.dat', 'c', '-ascii', '-double', '-tabs')




xlim([-60,60])
ylim([-350,ymax])



%% PRESSURE

m = dlmread('./../../build/bin/preform/pressureTime.txt',' ');
%m(ceil(length(m)/4):end,:) = smoothdata(m(ceil(length(m)/4):end,:))
figure
plot(m(:,1),m(:,2),'k-','linewidth',3);
xlabel('time');
ylabel('Pressure');
xlim([0,TMAX])
ylim([0,1])

m = m(1:10:end,:);

save('pressure_num.dat', 'm', '-ascii', '-double', '-tabs')

b = csvread("../preform/Experimental/N5_sr_pressure.csv");
hold on
plot(b(:,1)-0.008,b(:,2)/10,'b--')
b(:,2) = b(:,2)/10;


save('pressure_exp.dat', 'b', '-ascii', '-double', '-tabs')




legend('Mesh-free','Experiemntal')
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

legend({'Meshfree','Experimental'},... % { 'legend1', 'legend2',...}
'interpreter','latex',...
'FontSize',12,...
'FontName','cmr14',...
'Location','SouthEast')


set(gcf, 'Color', 'w');

print -dpng2 pressure.png
