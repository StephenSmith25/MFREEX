close all

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

plot_node = 107;



plotFiles = ceil(linspace(1,numFiles,10));


filename = strcat('displacement_',num2str(plotFiles(1)),'.txt');
disp = csvread(filename);
filename = strcat('srRod_',num2str(plotFiles(1)),'.txt');
dispRod = csvread(filename);

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
xlim([0,0.8])
ylim([0,0.9])

b = csvread("N2100_exp_pressure.csv")
hold on
plot(b(:,1)-0.014,b(:,2),'r--')
 


c = csvread("N2100_fe_pressure.csv")

hold on
plot(c(:,1)-0.014,c(:,2),'b--')

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
'Location','NorthEast')


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
numFiles = size(d,1)/2 - 2;
plotFiles = ceil(linspace(1,numFiles,100));
for i = 1:length(plotFiles)
    
    

    filename = strcat('strain_',num2str(plotFiles(i)),'.txt');

    strain = csvread(filename);
    
    F =  [strain(1,1) strain(1,2) 0 ;
        strain(2,1) strain(2,2) 0;
        0 0 strain(3,3)]  ;
    
    Fbar = det(F)^(-1/3)*(F);
    
    [R, U ,V] = poldecomp(Fbar);
    
    
    true_strain = logm(U);
    
    hoop_strain(i) =true_strain(3,3);
    axial_strain(i) = true_strain(2,2);
    time(i) = strain(3,2);
  
  
    
end



figure
hold on 
hold on
plot(time,(axial_strain),'k','linewidth',2);
% hold on
% plot(time,log(1+shear_strain),'k');

    
c = csvread("N2100_exp_strain_axial.csv")

hold on
plot(c(:,1),c(:,2),'b--')

c = csvread("N2100_fe_strain_axial.csv")

hold on
plot(c(:,1),c(:,2),'r*')

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

legend({'Axial (Meshfree)', 'Axial (exp)','Axial (FE)', },... % { 'legend1', 'legend2',...}
'FontUnits','points',...
'interpreter','latex',...
'FontSize',12,...
'FontName','cmr14',...
'Location','NorthWest')

grid on 



print -depsc2 strain_axial.eps

figure
plot(time,(hoop_strain),'k','linewidth',2);


c = csvread("N2100_exp_strain_hoop.csv")

hold on
plot(c(:,1),c(:,2),'b--')


c = csvread("N2100_fe_strain_hoop.csv")

hold on
plot(c(:,1),c(:,2),'r*')



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
'Location','NorthWest')

grid on 



print -depsc2 strain_hoop.eps

% 
% 
% for i = 1:length(plotFiles)
%     
%     
% 
%     filename = strcat('stress_',num2str(plotFiles(i)),'.txt');
% 
%     stress = csvread(filename);
%     hoop_stress(i) = stress(3,3);
%     axial_stress(i) = stress(2,2);
%     shear_stress(i) = stress(1,2);
%   
%     
% end
% 
% figure
% hold on 
% plot(1:length(hoop_stress),hoop_stress,'b');
% hold on
% plot(1:length(hoop_strain),axial_stress,'r');
% hold on
% plot(1:length(hoop_strain),shear_stress,'k');
% 
% 
% %set axis
% set(gca, 'FontName', 'cmr12')
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
% ylabel({'Cauchy Stress '},...
% 'FontUnits','points',...
% 'interpreter','latex',...
% 'FontSize',14,... % font size
% 'FontName','cmr14')
% % X label
% xlabel('Time',...
% 'FontUnits','points',...
% 'interpreter','latex',...
% 'FontWeight','normal',...
% 'FontSize',14,... % font size
% 'FontName','cmr14')
% % legend
% 
% legend({'Hoop', 'Axial', 'Shear' },... % { 'legend1', 'legend2',...}
% 'FontUnits','points',...
% 'interpreter','latex',...
% 'FontSize',12,...
% 'FontName','cmr14',...
% 'Location','NorthWest')
% 
% grid on 
% 
% 
% 
% print -depsc2 stress.eps


