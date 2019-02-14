clear all
close all

PLOT_GRAPHS = true;
JOB_ID = string(296756);
TMAX = 0.25;

path = strcat('/home/stephen/Documents/remote_belinda/Meshfree/',JOB_ID,'/Displacement/');
pathSR = strcat('/home/stephen/Documents/remote_belinda/Meshfree/',JOB_ID,'/srRod/');

boundaryNodesfile = strcat('/home/stephen/Documents/remote_belinda/Meshfree/',JOB_ID,'/boundary.txt');

boundaryNodes = csvread(boundaryNodesfile);
boundaryNodes = [boundaryNodes;boundaryNodes(1)];

addpath(path)
displacementdir = path ;
d = dir(displacementdir);
numFiles = size(d,1) -3 ;

plotFiles = ceil(linspace(1,numFiles,10));

plot_point =191;
filename = strcat(path,'displacement_',num2str(plotFiles(1)),'.csv');
disp = csvread(filename,1);




figure
subplot(1,3,1)       % add first plot in 2 x 2 grid
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
filename = strcat(pathSR,'srRod_',num2str(plotFiles(1)),'.csv');
disp = csvread(filename,1);
hold on
plot(disp(:,1),disp(:,2),'r-','linewidth',1)           % line plot
hold on
plot(-disp(:,1),disp(:,2),'r-','linewidth',1)           % line plot
xlim([-50,50])
ylim([-170,ymax])

       % line plot




filename = strcat(path,'displacement_',num2str(plotFiles(5)),'.csv');
disp = csvread(filename,1);
subplot(1,3,2)       % add first plot in 2 x 2 grid
plot(disp(:,1),disp(:,2),'k.','markersize',3)           % line plot
hold on
plot(-disp(:,1),disp(:,2),'k.','markersize',3) % line plot
hold on
plot(disp(plot_point,1),disp(plot_point,2),'r*')
axis equal
hold on
plot(disp(boundaryNodes,1),disp(boundaryNodes,2),'b--')
hold on
%plot(-disp(boundaryNodes,1),disp(boundaryNodes,2),'b--')
%plot(disp(boundaryNodes,1),disp(boundaryNodes,2),'b-')
hold on
filename = strcat(pathSR,'srRod_',num2str(plotFiles(5)),'.csv');
disp = csvread(filename,1);
hold on
plot(disp(:,1),disp(:,2),'r-','linewidth',1)           % line plot
hold on
plot(-disp(:,1),disp(:,2),'r-','linewidth',1)           % line plot
xlim([-50,50])
ylim([-170,ymax])

filename = strcat(path,'displacement_',num2str(plotFiles(10)),'.csv');
disp = csvread(filename,1);


subplot(1,3,3)       % add first plot in 2 x 2 grid
plot(disp(:,1),disp(:,2),'k.','markersize',3)           % line plot
hold on
plot(-disp(:,1),disp(:,2),'k.','markersize',3)           % line plot

hold on
plot(disp(plot_point,1),disp(plot_point,2),'r*')
axis equal
hold on
%plot(disp(boundaryNodes,1),disp(boundaryNodes,2),'b-')
hold on
%plot(-disp(boundaryNodes,1),disp(boundaryNodes,2),'b-')
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

file_pressure = strcat('/home/stephen/Documents/remote_belinda/Meshfree/',JOB_ID,'/pressureTime.txt');

m = dlmread(file_pressure,' ');
%m(ceil(length(m)/4):end,:) = smoothdata(m(ceil(length(m)/4):end,:))
figure
plot(m(:,1),m(:,2),'k-','linewidth',3);
xlabel('time');
ylabel('Pressure');
xlim([0,TMAX])
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

if ( PLOT_GRAPHS)
    
path = strcat('/home/stephen/Documents/remote_belinda/Meshfree/',JOB_ID,'/History/');
path_strain = strcat(path,'Strain/');

addpath(path)
displacementdir = path_strain ;
d = dir(displacementdir);
numFiles = size(d,1) -3 ;

plotFiles = ceil(linspace(1,numFiles,100));
for i = 1:length(plotFiles)
    
    

    filename = strcat(path_strain,'strain_',num2str(plotFiles(i)),'.txt');

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


strain = figure('Position', get(0, 'Screensize'));
subplot(1,2,1)
plot(time,(axial_strain),'k','linewidth',2);
% hold on
% plot(time,log(1+shear_strain),'k');
xlim([0,TMAX])

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
'FontWeight','normal',...
'FontSize',14,... % size ofiguref numbers on axis
'FontName','cmr14') % font name
set(gca,'TickLabelInterpreter', 'latex');

% Y label
ylabel({'True Strain '},...
'interpreter','latex',...
'FontSize',14,... % font size
'FontName','cmr14')
% X label
xlabel('Time',...
'interpreter','latex',...
'FontWeight','normal',...
'FontSize',14,... % font size
'FontName','cmr14')
% legend

legend({'Axial (Meshfree)', 'Axial (exp)','Axial (FE)', 'Jacobian', 'Shear' },... % { 'legend1', 'legend2',...}
'interpreter','latex',...
'FontSize',12,...
'FontName','cmr14',...
'Location','SouthEast')

grid on 




subplot(1,2,2)
plot(time,(hoop_strain),'k','linewidth',2);


c = csvread("Experimental/N4100_exp_strain_hoop.csv");

hold on
plot(c(:,1),c(:,2),'b--')


c = csvread("Experimental/N4100_fe_strain_hoop.csv");

hold on
plot(c(:,1),c(:,2),'r--')

xlim([0,TMAX])


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
ylabel({'True Strain '},...
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
'interpreter','latex',...
'FontSize',12,...
'FontName','cmr14',...
'Location','SouthEast')

grid on 



print -dpng2 strain.png

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
'FontWeight','normal',...
'FontSize',14,... % size ofiguref numbers on axis
'FontName','cmr14') % font name
set(gca,'TickLabelInterpreter', 'latex');

% Y label
ylabel({'Cauchy Stress '},...
'interpreter','latex',...
'FontSize',14,... % font size
'FontName','cmr14')
% X label
xlabel('Time',...
'interpreter','latex',...
'FontWeight','normal',...
'FontSize',14,... % font size
'FontName','cmr14')
% legend

legend({ 'Hoop (Conf)','Hoop (Bond)','Shear (Bond)', 'Shear (Conf)' },... % { 'legend1', 'legend2',...}
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
'FontWeight','normal',...
'FontSize',14,... % size ofiguref numbers on axis
'FontName','cmr14') % font name
set(gca,'TickLabelInterpreter', 'latex');

% Y label
ylabel({'Cauchy Stress '},...
'interpreter','latex',...
'FontSize',14,... % font size
'FontName','cmr14')
% X label
xlabel('Time',...
'interpreter','latex',...
'FontWeight','normal',...
'FontSize',14,... % font size
'FontName','cmr14')
% legend

legend({ 'Axial (Conf)','Axial (Bond)' },... % { 'legend1', 'legend2',...}
'interpreter','latex',...
'FontSize',12,...
'FontName','cmr14',...
'Location','SouthEast')

grid on 



print -dpng2 stress.png



 end