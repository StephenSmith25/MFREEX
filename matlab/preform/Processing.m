clear all
close all

path = './../../build/bin/preform/Displacement/';
pathSR = './../../build/bin/preform/srRod/';
boundaryNodes = csvread('./../../build/bin/preform/boundary.txt');
boundaryNodes = [boundaryNodes;boundaryNodes(1)];

addpath(path)
displacementdir = path ;
d = dir(displacementdir);
d1 = dir([displacementdir,'*.txt']);
numFiles = size(d,1) -3 ;

plotFiles = ceil(linspace(1,numFiles,10));

plot_point =59;
filename = strcat(path,'displacement_',num2str(plotFiles(1)),'.txt');
disp = csvread(filename);
figure
subplot(2,3,1)       % add first plot in 2 x 2 grid
plot(disp(:,1),disp(:,2),'k.')           % line plot
hold on
plot(disp(plot_point,1),disp(plot_point,2),'r*')
axis equal 
hold on
plot(disp(boundaryNodes,1),disp(boundaryNodes,2),'b-')

hold on
filename = strcat(pathSR,'srRod_',num2str(plotFiles(1)),'.txt');
disp = csvread(filename);
hold on
plot(disp(:,1),disp(:,2),'r-','linewidth',1)           % line plot




filename= strcat(path,'displacement_',num2str(plotFiles(2)),'.txt');
disp = csvread(filename);

subplot(2,3,2)       % add first plot in 2 x 2 grid
plot(disp(:,1),disp(:,2),'k.')           % line plot
hold on
plot(disp(plot_point,1),disp(plot_point,2),'r*')
axis equal 
hold on
plot(disp(boundaryNodes,1),disp(boundaryNodes,2),'b-')
hold on
filename = strcat(pathSR,'srRod_',num2str(plotFiles(2)),'.txt');
disp = csvread(filename);
hold on
plot(disp(:,1),disp(:,2),'r-','linewidth',1)           % line plot


filename = strcat(path,'displacement_',num2str(plotFiles(4)),'.txt');
disp = csvread(filename);


subplot(2,3,3)       % add first plot in 2 x 2 grid
plot(disp(:,1),disp(:,2),'k.')           % line plot
hold on
plot(disp(plot_point,1),disp(plot_point,2),'r*')
axis equal
hold on
plot(disp(boundaryNodes,1),disp(boundaryNodes,2),'b-')
hold on
filename = strcat(pathSR,'srRod_',num2str(plotFiles(4)),'.txt');
disp = csvread(filename);
hold on
plot(disp(:,1),disp(:,2),'r-','linewidth',1)           % line plot





filename = strcat(path,'displacement_',num2str(plotFiles(7)),'.txt');
disp = csvread(filename);
subplot(2,3,4)       % add first plot in 2 x 2 grid
plot(disp(:,1),disp(:,2),'k.')           % line plot
hold on
plot(disp(plot_point,1),disp(plot_point,2),'r*')
axis equal
hold on
plot(disp(boundaryNodes,1),disp(boundaryNodes,2),'b-')
hold on
filename = strcat(pathSR,'srRod_',num2str(plotFiles(7)),'.txt');
disp = csvread(filename);
hold on
plot(disp(:,1),disp(:,2),'r-','linewidth',1)           % line plot




filename = strcat(path,'displacement_',num2str(plotFiles(8)),'.txt');
disp = csvread(filename);
subplot(2,3,5)       % add first plot in 2 x 2 grid
plot(disp(:,1),disp(:,2),'k.')           % line plot
hold on
plot(disp(plot_point,1),disp(plot_point,2),'r*')
axis equal
hold on
plot(disp(boundaryNodes,1),disp(boundaryNodes,2),'b-')
hold on
filename = strcat(pathSR,'srRod_',num2str(plotFiles(8)),'.txt');
disp = csvread(filename);
hold on
plot(disp(:,1),disp(:,2),'r-','linewidth',1)           % line plot


filename = strcat(path,'displacement_',num2str(plotFiles(10)),'.txt');
disp = csvread(filename);


subplot(2,3,6)       % add first plot in 2 x 2 grid
plot(disp(:,1),disp(:,2),'k.')           % line plot
hold on
plot(disp(plot_point,1),disp(plot_point,2),'r*')
axis equal
hold on
plot(disp(boundaryNodes,1),disp(boundaryNodes,2),'b-')
hold on
filename = strcat(pathSR,'srRod_',num2str(plotFiles(10)),'.txt');
disp = csvread(filename);
hold on
plot(disp(:,1),disp(:,2),'r-','linewidth',1)           % line plot



m = dlmread('./../../build/bin/preform/pressureTime.txt',' ');
%m(ceil(length(m)/4):end,:) = smoothdata(m(ceil(length(m)/4):end,:))
figure
plot(m(:,1),m(:,2),'k-','linewidth',3);
xlabel('time');
ylabel('Pressure');
xlim([0,0.5])
ylim([0,0.9])

b = csvread("Experimental/N4T100_exp_pressure.csv");
hold on
plot(b(:,1),b(:,2)/1000,'b--')


c = csvread("Experimental/N4T100_fe_pressure.csv");

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

path = './../../build/bin/preform/History/Strain';

addpath(path)
displacementdir = path ;
d = dir(displacementdir);
d1 = dir([displacementdir,'*.txt']);
numFiles = size(d,1) -3 ;

plotFiles = ceil(linspace(1,numFiles,100));
for i = 1:length(plotFiles)
    
    

    filename = strcat('./../../build/bin/preform/History/Strain/strain_',num2str(plotFiles(i)),'.txt');

    strain = csvread(filename);
    
 
    F =  [strain(1,1) strain(1,2) 0 ;
        strain(2,1) strain(2,2) 0;
        0 0 strain(3,3)]  ;
    
    
    [R U V] = poldecomp(F);
    
    
    true_strain = logm(V);
    
    
    hoop_strain(i) =true_strain(3,3);
    axial_strain(i) = true_strain(2,2);
    radial_strain(i) = true_strain(1,1);
    shear_strain(i) = true_strain(1,2);
    jacobian(i) = det(F);
    crit_lambda(i) = strain(2,3);
    gamma(i) = strain(3,1);
    max_lambda_n(i) = strain(1,3);
    time(i) = strain(3,2);
  
    
end



figure
subplot(1,2,1)
plot(time,(axial_strain),'k','linewidth',2);
% hold on
% plot(time,log(1+shear_strain),'k');
xlim([0,0.5])

c = csvread("Experimental/N4100_exp_strain_axial.csv");

hold on
plot(c(:,1),c(:,2),'b--')

c = csvread("Experimental/N4100_fe_strain_axial.csv");

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


c = csvread("Experimental/N4100_exp_strain_hoop.csv");

hold on
plot(c(:,1),c(:,2),'b--')


c = csvread("Experimental/N4100_fe_strain_hoop.csv");

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

figure 
subplot(1,2,1);
plot(time,crit_lambda,'k-');
hold on
plot(time,max_lambda_n,'b-');

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
ylabel({'Critical Network Stretch ($\lambda$)'},...
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


 
subplot(1,2,2);
plot(time,gamma,'k-');
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
ylabel({'Viscosity ($\gamma$) Pa'},...
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
path = './../../build/bin/preform/History/Stress';





addpath(path)
displacementdir = path ;
d = dir(displacementdir);
d1 = dir([displacementdir,'*.txt']);
numFiles = size(d,1)/2 -3 ;

plotFiles = ceil(linspace(1,numFiles,100)); 

for i = 1:length(plotFiles)
    
    

    filename = strcat('Bond_Stress_',num2str(plotFiles(i)),'.txt');

    stress = csvread(filename);
    hoop_stress_bond(i) = stress(3,3);
    axial_stress_bond(i) = stress(2,2);
    radial_stress_bond(i) = stress(1,1);
    shear_stress_bond(i) = stress(1,2);

    filename = strcat('Conformational_Stress_',num2str(plotFiles(i)),'.txt');

    stress = csvread(filename);
    hoop_stress_conf(i) = stress(3,3);
    axial_stress_conf(i) = stress(2,2);
    radial_stress_conf(i) = stress(1,1);
     shear_stress_conf(i) = stress(1,2);

  
end

figure

subplot(1,2,1)
hold on 
plot(time,hoop_stress_conf,'r');
hold on
plot(time,hoop_stress_bond,'b');
hold on
plot(time,shear_stress_bond,'m')
hold on
plot(time,shear_stress_conf,'g')
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

legend({ 'Hoop (Conf)','Hoop (Bond)','Shear (Bond)', 'Shear (Conf)' },... % { 'legend1', 'legend2',...}
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




