clear all
close all 

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
subplot(2,3,1)       % add first plot in 2 x 2 grid
plot(disp(:,1),disp(:,2),'k.')           % line plot
hold on
plot(-disp(:,1), disp(:,2),'k.');

plate_width = 3*max(disp(:,1)) ;
hold on
plot([-plate_width,plate_width],[0,0],'r');

max_Y = max(disp(:,2));
min_Y = min(disp(:,2));

max_X = max(disp(:,2));
min_X = -max(disp(:,2));
ylim([min_Y - 0.1*max_Y,2*max_Y])
xlim([-plate_width - 0.1*plate_width,plate_width + 0.1*plate_width])




% second plot
filename = strcat(path,'/displacement_',num2str(plotFiles(2)),'.txt');
disp = csvread(filename);


subplot(2,3,2)       % add first plot in 2 x 2 grid
plot(disp(:,1),disp(:,2),'k.')           % line plot
 

hold on
plot(-disp(:,1), disp(:,2),'k.');
hold on
plot([-plate_width,plate_width],[0,0],'r');
ylim([min_Y - 0.1*max_Y,2*max_Y])
xlim([-plate_width - 0.1*plate_width,plate_width + 0.1*plate_width])

    

% third plot
filename = strcat(path,'/displacement_',num2str(plotFiles(4)),'.txt');
disp = csvread(filename);


subplot(2,3,3)       % add first plot in 2 x 2 grid
plot(disp(:,1),disp(:,2),'k.')           % line plot

hold on
plot(-disp(:,1), disp(:,2),'k.');

hold on
plot([-plate_width,plate_width],[0,0],'r');

ylim([min_Y - 0.1*max_Y,2*max_Y])
xlim([-plate_width - 0.1*plate_width,plate_width + 0.1*plate_width])



% fourth plot
filename = strcat(path,'/displacement_',num2str(plotFiles(7)),'.txt');
disp = csvread(filename);

subplot(2,3,4)       % add first plot in 2 x 2 grid
plot(disp(:,1),disp(:,2),'k.')           % line plot

hold on
plot(-disp(:,1), disp(:,2),'k.');

hold on
plot([-plate_width,plate_width],[0,0],'r');
ylim([min_Y - 0.1*max_Y,2*max_Y])
xlim([-plate_width - 0.1*plate_width,plate_width + 0.1*plate_width])



filename = strcat(path,'/displacement_',num2str(plotFiles(8)),'.txt');
disp = csvread(filename);

% 5th plot
subplot(2,3,5)       % add first plot in 2 x 2 grid
plot(disp(:,1),disp(:,2),'k.')           % line plot


hold on
plot(-disp(:,1), disp(:,2),'k.');

hold on
plot([-plate_width,plate_width],[0,0],'r');

ylim([min_Y - 0.1*max_Y,2*max_Y])
xlim([-plate_width - 0.1*plate_width,plate_width + 0.1*plate_width])



% 6th plot
filename = strcat(path,'/displacement_',num2str(plotFiles(9)),'.txt');
disp = csvread(filename);


subplot(2,3,6)       % add first plot in 2 x 2 grid
plot(disp(:,1),disp(:,2),'k.')           % line plot

hold on
plot(-disp(:,1), disp(:,2),'k.');

hold on
plot([-plate_width,plate_width],[0,0],'r');
xlim([-plate_width - 0.1*plate_width,plate_width + 0.1*plate_width])

ylim([min_Y - 0.1*max_Y,2*max_Y])
