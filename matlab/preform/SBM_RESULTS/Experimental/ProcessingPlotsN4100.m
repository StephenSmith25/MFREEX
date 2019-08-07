close all
clear all 

set(0,'defaulttextinterpreter','latex')
displacementdir = './Displacement';
d = dir(displacementdir);
d1 = dir([displacementdir,'*.txt']);
numFiles = size(d,1) - 2;


%
srdir = './Displacement';
s = dir(srdir);
s1 = dir([srdir,'*.txt']);


boundaryNodes = csvread('boundary.txt');
boundaryNodes = [boundaryNodes;boundaryNodes(1)];



plot_node = 68:78;


plotFiles = ceil(linspace(1,numFiles,10));


filename = strcat('displacement_',num2str(plotFiles(1)),'.txt');
disp = csvread(filename);
filename = strcat('srRod_',num2str(plotFiles(1)),'.txt');
dispRod = csvread(filename);


figure('position', [150, 150, 1200, 1200])

p = plot(disp(:,1),disp(:,2),disp(boundaryNodes,1),disp(boundaryNodes,2),dispRod(:,1),dispRod(:,2))   ;        % line plot
p(1).Marker = 'o';
p(1).LineStyle = 'none';
p(1).MarkerSize = 6;
p(2).Color = 'k';
p(3).LineWidth = 2;
hold on
plot(disp(plot_node,1),disp(plot_node,2),'r+');
axis equal

%set axis
set(gca, 'FontName', 'cmr14')
% set x tics and y tics
set(gca,...
'Units','normalized',...
'FontUnits','points',...
'FontWeight','normal',...
'FontSize',14,... % size ofiguref numbers on axis
'FontName','cmr14') % font name
set(gca,'TickLabelInterpreter', 'latex');
xticks([0 10 20 30 40 50])
grid on
print -depsc2 model.eps


figure('position', [150, 150, 1200, 1200])
subplot(1,3,1)       % add first plot in 2 x 2 grid

p = plot(disp(:,1),disp(:,2),disp(boundaryNodes,1),disp(boundaryNodes,2),dispRod(:,1),dispRod(:,2))   ;        % line plot
p(1).Marker = '.';
p(1).LineStyle = 'none';
p(1).MarkerSize = 5;
p(2).Color = 'k';
p(3).LineWidth = 2;
hold on
plot(disp(plot_node,1),disp(plot_node,2),'r+');
axis equal
xlim([0,50])
ylim([-150,77])
%set axis
set(gca, 'FontName', 'cmr14')
% set x tics and y tics
set(gca,...
'Units','normalized',...
'FontUnits','points',...
'FontWeight','normal',...
'FontSize',14,... % size ofiguref numbers on axis
'FontName','cmr14') % font name
set(gca,'TickLabelInterpreter', 'latex');
xticks([0 10 20 30 40 50])
grid on


filename = strcat('displacement_',num2str(plotFiles(5)),'.txt');
disp = csvread(filename);
filename = strcat('srRod_',num2str(plotFiles(5)),'.txt');
dispRod = csvread(filename);


subplot(1,3,2)       % add first plot in 2 x 2 grid
p = plot(disp(:,1),disp(:,2),disp(boundaryNodes,1),disp(boundaryNodes,2),dispRod(:,1),dispRod(:,2))   ;        % line plot
p(1).Marker = '.';
p(1).LineStyle = 'none';
p(1).MarkerSize = 5;
p(2).Color = 'k';
p(3).LineWidth = 2;
hold on 
plot(disp(plot_node,1),disp(plot_node,2),'r+');

axis equal
xlim([0,50])
ylim([-150,77])
%set axis
set(gca, 'FontName', 'cmr14')
% set x tics and y tics
set(gca,...
'Units','normalized',...
'FontUnits','points',...
'FontWeight','normal',...
'FontSize',14,... % size ofiguref numbers on axis
'FontName','cmr14') % font name
set(gca,'TickLabelInterpreter', 'latex');
xticks([0 10 20 30 40 50])
grid on


filename = strcat('displacement_',num2str(plotFiles(10)),'.txt');
disp = csvread(filename);
filename = strcat('srRod_',num2str(plotFiles(10)),'.txt');
dispRod = csvread(filename);


subplot(1,3,3)       % add first plot in 2 x 2 grid
p = plot(disp(:,1),disp(:,2),disp(boundaryNodes,1),disp(boundaryNodes,2),dispRod(:,1),dispRod(:,2))   ;        % line plot
p(1).Marker = '.';
p(1).LineStyle = 'none';
p(1).MarkerSize = 5;
p(2).Color = 'k';
p(3).LineWidth = 2;
hold on 
plot(disp(plot_node,1),disp(plot_node,2),'r+');

axis equal
xlim([0,50])
ylim([-150,77])

%set axis
set(gca, 'FontName', 'cmr14')
% set x tics and y tics
set(gca,...
'Units','normalized',...
'FontUnits','points',...
'FontWeight','normal',...
'FontSize',14,... % size ofiguref numbers on axis
'FontName','cmr14') % font name
set(gca,'TickLabelInterpreter', 'latex');
xticks([0 10 20 30 40 50])
grid on
set(gcf, 'Color', 'w');

print -depsc2 deformation.eps

%export_fig('./Deformation.eps','-zbuffer','-r300')

%% sweep 
% 
for i = 1:length(disp)-1
    bf(i,:) = [i,i+1] ;
    
end

m = dlmread('pressureTime.txt',' ');
%m(ceil(length(m)/4):end,:) = smoothdata(m(ceil(length(m)/4):end,:))
figure
plot(m(:,1),m(:,2),'k-','linewidth',3);
xlabel('time');
ylabel('Pressure');
xlim([0,0.5])
ylim([0,0.9])

b = csvread("N4T100_exp_pressure.csv")
hold on
plot(b(:,1),b(:,2)/1000,'b--')


c = csvread("N4T100_fe_pressure.csv")

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
'FontUnits','points',...
'FontWeight','normal',...
'FontSize',14,... % size ofiguref numbers on axis
'FontName','cmr14') % font name
set(gca,'TickLabelInterpreter', 'latex');

% Y label
ylabel({'Pressure (MPa) '},...
'FontUnits','points',...
'interpreter','latex',...
'FontSize',14,... % font size
'FontName','cmr14')
% X label
xlabel('Time (s)',...
'FontUnits','points',...
'interpreter','latex',...
'FontWeight','normal',...
'FontSize',14,... % font size
'FontName','cmr14')
% legend

legend({'Meshfree','Experimental','FE (Abaqus )'},... % { 'legend1', 'legend2',...}
'FontUnits','points',...
'interpreter','latex',...
'FontSize',12,...
'FontName','cmr14',...
'Location','SouthEast')


set(gcf, 'Color', 'w');

print -depsc2 pressure.eps
set(gcf, 'Color', 'w');
%export_fig('./Pressure.png','-zbuffer','-r600')
%legend('nodes','ux = 0', 'uy = 0','Prescribed Displacement')

% %% tempreature
% d = csvread('IRprofile.csv');
% 
% figure
% plot(d(:,1),d(:,2),'k--','linewidth',2);
% hold on
% plot(d(:,1),d(:,3),'k.-','linewidth',2);
% 
% d = csvread('IRprofile.csv');
% 
% xlim([0,85])
% ylim([80,110])
% 
% % set axis
%  set(gca, 'FontName', 'cmr12')
% % set x tics and y tics
% set(gca,...
% 'Units','normalized',...
% 'FontUnits','points',...
% 'FontWeight','normal',...
% 'FontSize',14,... % size ofiguref numbers on axis
% 'FontName','cmr14') % font name
% set(gca,'TickLabelInterpreter', 'latex');
% 
% % Y label
% ylabel({'Temperature ($^{\circ}$C) '},...
% 'FontUnits','points',...
% 'interpreter','latex',...
% 'FontSize',14,... % font size
% 'FontName','cmr14')
% % X label
% xlabel('Distance from the neck (mm)',...
% 'FontUnits','points',...
% 'interpreter','latex',...
% 'FontWeight','normal',...
% 'FontSize',14,... % font size
% 'FontName','cmr14')
% % legend
% 
% legend({'Outer','Inner'},... % { 'legend1', 'legend2',...}
% 'FontUnits','points',...
% 'interpreter','latex',...
% 'FontSize',12,...
% 'FontName','cmr14',...
% 'Location','NorthEast')
% 
% 
% print -depsc2 Temperature.eps
% 
% set(gcf, 'Color', 'w');
% export_fig('./Temperature.png','-zbuffer','-r600')
% 
% 
% filename = strcat('displacement_',num2str(plotFiles(10)),'.txt');
% 
% 
% disp = csvread(filename);
% % sweep the profile
% [tri, xyz] = bfRevolve(bf, disp, 20);
%  
% % display the surface
% figure
% h= trisurf(tri, xyz(:,1),xyz(:,2),xyz(:,3));
% set(h,'edgecolor','k','facecolor',[0 0.8 0.6])
% hold on
% plot3(xyz(:,1),xyz(:,2),xyz(:,3),'k.','markersize',3)
% axis equal
% axis([-50 50 -50 50 -200 max(disp(:,2))])
% view(45,30);
% % colormap summer
% shading interp
% yourColorMap(1, :) = [0.9,0.9,0.9];
% colormap(yourColorMap);
% set(gcf, 'Color', 'w');
% export_fig('./Bottle.eps','-zbuffer','-r600')
% print -depsc2 bottle.eps
% 
Historydir = './History';
d = dir(Historydir);
d1 = dir([Historydir,'*.txt']);
numFiles = size(d,1)/3 - 2;
plotFiles = ceil(linspace(1,numFiles,100));
for i = 1:length(plotFiles)
    
    

    filename = strcat('strain_',num2str(plotFiles(i)),'.txt');

    strain = csvread(filename);
    
 
    F =  [strain(1,1) strain(1,2) 0 ;
        strain(2,1) strain(2,2) 0;
        0 0 strain(3,3)]  ;
    
    
    [R U V] = poldecomp(F);
    
    
    true_strain = logm(U);
    
    
    hoop_strain(i) =true_strain(3,3);
    axial_strain(i) = true_strain(2,2);
    radial_strain(i) = true_strain(1,1);
    shear_strain(i) = true_strain(1,2);
    jacobian(i) = det(F);

    time(i) = strain(3,2);
  
    
end



figure
subplot(1,2,1)
plot(time,(axial_strain),'k','linewidth',2);
% hold on
% plot(time,log(1+shear_strain),'k');
xlim([0,0.5])

c = csvread("N4100_exp_strain_axial.csv");

hold on
plot(c(:,1),c(:,2),'b--')

c = csvread("N4100_fe_strain_axial.csv");

hold on
plot(c(:,1),c(:,2),'r--')


hold on
plot(time,jacobian,'g--')

hold on
plot(time,(shear_strain),'g','linewidth',2);


%set axis
set(gca, 'FontName', 'cmr12')
% set x tics and y tics
set(gca,...
'Units','normalized',...
'FontUnits','points',...
'FontWeight','normal',...
'FontSize',14,... % size ofiguref numbers on axis
'FontName','cmr14') % font name
set(gca,'TickLabelInterpreter', 'latex');

% Y label
ylabel({'True Strain '},...
'FontUnits','points',...
'interpreter','latex',...
'FontSize',14,... % font size
'FontName','cmr14')
% X label
xlabel('Time',...
'FontUnits','points',...
'interpreter','latex',...
'FontWeight','normal',...
'FontSize',14,... % font size
'FontName','cmr14')
% legend

legend({'Axial (Meshfree)', 'Axial (exp)','Axial (FE)', 'Jacobian', 'Shear' },... % { 'legend1', 'legend2',...}
'FontUnits','points',...
'interpreter','latex',...
'FontSize',12,...
'FontName','cmr14',...
'Location','SouthEast')

grid on 



print -depsc2 strain_axial.eps

subplot(1,2,2)
plot(time,(hoop_strain),'k','linewidth',2);


c = csvread("N4100_exp_strain_hoop.csv");

hold on
plot(c(:,1),c(:,2),'b--')


c = csvread("N4100_fe_strain_hoop.csv");

hold on
plot(c(:,1),c(:,2),'r--')

xlim([0,0.5])


%set axis
set(gca, 'FontName', 'cmr12')
% set x tics and y tics
set(gca,...
'Units','normalized',...
'FontUnits','points',...
'FontWeight','normal',...
'FontSize',14,... % size ofiguref numbers on axis
'FontName','cmr14') % font name
set(gca,'TickLabelInterpreter', 'latex');

% Y label
ylabel({'True Strain '},...
'FontUnits','points',...
'interpreter','latex',...
'FontSize',14,... % font size
'FontName','cmr14')
% X label
xlabel('Time',...
'FontUnits','points',...
'interpreter','latex',...
'FontWeight','normal',...
'FontSize',14,... % font size
'FontName','cmr14')
% legend

legend({'Hoop (Meshfree)','Hoop (exp)','Hoop (FE)' },... % { 'legend1', 'legend2',...}
'FontUnits','points',...
'interpreter','latex',...
'FontSize',12,...
'FontName','cmr14',...
'Location','SouthEast')

grid on 



print -depsc2 strain_hoop.eps




for i = 1:length(plotFiles)
    
    

    filename = strcat('stress_bond_',num2str(plotFiles(i)),'.txt');

    stress = csvread(filename);
    hoop_stress_bond(i) = stress(3,3);
    axial_stress_bond(i) = stress(2,2);
    radial_stress_bond(i) = stress(1,1);
  
    filename = strcat('stress_conf_',num2str(plotFiles(i)),'.txt');

    stress = csvread(filename);
    hoop_stress_conf(i) = stress(3,3);
    axial_stress_conf(i) = stress(2,2);
    radial_stress_conf(i) = stress(1,1);
    
  
end

figure

subplot(1,2,1)
hold on 
plot(time,hoop_stress_conf,'r');
hold on
plot(time,hoop_stress_bond,'b');
hold on
plot(time,log(jacobian)*1.8e9,'m')

%set axis
set(gca, 'FontName', 'cmr12')
% set x tics and y tics
set(gca,...
'Units','normalized',...
'FontUnits','points',...
'FontWeight','normal',...
'FontSize',14,... % size ofiguref numbers on axis
'FontName','cmr14') % font name
set(gca,'TickLabelInterpreter', 'latex');

% Y label
ylabel({'Cauchy Stress '},...
'FontUnits','points',...
'interpreter','latex',...
'FontSize',14,... % font size
'FontName','cmr14')
% X label
xlabel('Time',...
'FontUnits','points',...
'interpreter','latex',...
'FontWeight','normal',...
'FontSize',14,... % font size
'FontName','cmr14')
% legend

legend({ 'Hoop (Conf)','Hoop (Bond)','Volumetric' },... % { 'legend1', 'legend2',...}
'FontUnits','points',...
'interpreter','latex',...
'FontSize',12,...
'FontName','cmr14',...
'Location','SouthEast')

grid on 

subplot(1,2,2)
plot(time,axial_stress_conf,'r');
hold on
plot(time,axial_stress_bond,'b');

%set axis
set(gca, 'FontName', 'cmr12')
% set x tics and y tics
set(gca,...
'Units','normalized',...
'FontUnits','points',...
'FontWeight','normal',...
'FontSize',14,... % size ofiguref numbers on axis
'FontName','cmr14') % font name
set(gca,'TickLabelInterpreter', 'latex');

% Y label
ylabel({'Cauchy Stress '},...
'FontUnits','points',...
'interpreter','latex',...
'FontSize',14,... % font size
'FontName','cmr14')
% X label
xlabel('Time',...
'FontUnits','points',...
'interpreter','latex',...
'FontWeight','normal',...
'FontSize',14,... % font size
'FontName','cmr14')
% legend

legend({ 'Axial (Conf)','Axial (Bond)' },... % { 'legend1', 'legend2',...}
'FontUnits','points',...
'interpreter','latex',...
'FontSize',12,...
'FontName','cmr14',...
'Location','SouthEast')

grid on 



print -depsc2 stress.eps


