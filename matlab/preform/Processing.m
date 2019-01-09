clear all
close all

path = './../../build/bin/preform/Displacement';

addpath(path)
displacementdir = path ;
d = dir(displacementdir);
d1 = dir([displacementdir,'*.txt']);
numFiles = size(d,1) -3 ;

plotFiles = ceil(linspace(1,numFiles,10));


filename = strcat(path,'/displacement_',num2str(plotFiles(1)),'.txt');
disp = csvread(filename);
figure
subplot(2,3,1)       % add first plot in 2 x 2 grid
plot(disp(:,1),disp(:,2),'k.')           % line plot
axis equal 

filename = strcat(path,'/displacement_',num2str(plotFiles(2)),'.txt');
disp = csvread(filename);


subplot(2,3,2)       % add first plot in 2 x 2 grid
plot(disp(:,1),disp(:,2),'k.')           % line plot
 

axis equal 



filename = strcat(path,'/displacement_',num2str(plotFiles(4)),'.txt');
disp = csvread(filename);


subplot(2,3,3)       % add first plot in 2 x 2 grid
plot(disp(:,1),disp(:,2),'k.')           % line plot
axis equal


filename = strcat(path,'/displacement_',num2str(plotFiles(7)),'.txt');
disp = csvread(filename);


subplot(2,3,4)       % add first plot in 2 x 2 grid
plot(disp(:,1),disp(:,2),'k.')           % line plot
axis equal



filename = strcat(path,'/displacement_',num2str(plotFiles(8)),'.txt');
disp = csvread(filename);


subplot(2,3,5)       % add first plot in 2 x 2 grid
plot(disp(:,1),disp(:,2),'k.')           % line plot
axis equal

filename = strcat(path,'/displacement_',num2str(plotFiles(10)),'.txt');
disp = csvread(filename);


subplot(2,3,6)       % add first plot in 2 x 2 grid
plot(disp(:,1),disp(:,2),'k.')           % line plot
axis equal


