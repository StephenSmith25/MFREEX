clear all
close all 

path = './../../build/bin/taylor_bar/Displacement';

addpath(path)
displacementdir = path ;
d = dir(displacementdir);
d1 = dir([displacementdir,'*.txt']);
numFiles = size(d,1) -3 ;

plotFiles = ceil(linspace(1,numFiles,10));

% first plot
filename = strcat(path,'/displacement_',num2str(plotFiles(1)),'.txt');
disp = csvread(filename);
figure
subplot(1,3,1)       % add first plot in 2 x 2 grid

plot(disp(:,1),disp(:,2),'k.','markersize',3)           % line plot

hold on
plot(-disp(:,1), disp(:,2),'k.','markersize',3);
plate_width = 3*max(disp(:,1)) ;
hold on
axis equal

plot([-plate_width,plate_width],[0,0],'r');

max_Y = max(disp(:,2));
min_Y = min(disp(:,2));

max_X = max(disp(:,2));
min_X = -max(disp(:,2));
ylim([min_Y - 0.1*max_Y,1.2*max_Y])
xlim([-plate_width - 0.1*plate_width,plate_width + 0.1*plate_width])



% second plot
filename = strcat(path,'/displacement_',num2str(plotFiles(4)),'.txt');
disp = csvread(filename);


subplot(1,3,2)       % add first plot in 2 x 2 grid
plot(disp(:,1),disp(:,2),'k.','markersize',3)           % line plot

hold on
plot(-disp(:,1), disp(:,2),'k.','markersize',3);
hold on
plot([-plate_width,plate_width],[0,0],'r');

axis equal


ylim([min_Y - 0.1*max_Y,1.2*max_Y])
xlim([-plate_width - 0.1*plate_width,plate_width + 0.1*plate_width])

    

% third plot
filename = strcat(path,'/displacement_',num2str(plotFiles(10)),'.txt');
disp = csvread(filename);


subplot(1,3,3)       % add first plot in 2 x 2 grid
plot(disp(:,1),disp(:,2),'k.','markersize',3)           % line plot

hold on
plot(-disp(:,1), disp(:,2),'k.','markersize',3);

hold on
plot([-plate_width,plate_width],[0,0],'r-');
axis equal
ylim([min_Y - 0.1*max_Y,1.2*max_Y])
xlim([-plate_width - 0.1*plate_width,plate_width + 0.1*plate_width])


print(gcf,'Displacement_taylor.png','-dpng','-r300');        

path = './../../build/bin/taylor_bar/History/Stress';




addpath(path)
displacementdir = path ;
d = dir(displacementdir);
d1 = dir([displacementdir,'*.txt']);
numFiles = size(d,1)/2 -3 ;


for i = 1:length(plotFiles)
    
    

    filename = strcat('Stress_',num2str(plotFiles(i)),'.txt');

    stress = csvread(filename);
    hoop_stress(i) = stress(3,3)-(1/3)*trace(stress);
    axial_stress(i) = stress(2,2) - (1/3)*trace(stress);
    radial_stress(i) = stress(1,1) - (1/3)*trace(stress);
    shear_stress(i) = stress(1,2);

 

  
end

time = [1:1:length(hoop_stress)];

figure
hold on 
plot(time,hoop_stress,'g');
hold on 
plot(time,radial_stress,'m');
hold on

plot(time,axial_stress,'r');
hold on
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

legend({ 'Hoop', 'Radial', 'Axial' },... % { 'legend1', 'legend2',...}
'interpreter','latex',...
'FontSize',12,...
'FontName','cmr14',...
'Location','SouthEast')
%