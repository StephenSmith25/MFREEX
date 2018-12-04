close all
clear all 


path = './../../build/bin/beamUL/Displacement';

addpath(path)
displacementdir = path ;
d = dir(displacementdir);
d1 = dir([displacementdir,'*.txt']);
numFiles = size(d,1) -3 ;

plotFiles = ceil(linspace(1,numFiles,10));


filename = strcat('displacement_',num2str(plotFiles(1)),'.txt');
disp = csvread(filename);
figure
subplot(2,3,1)       % add first plot in 2 x 2 grid
plot(disp(:,1),disp(:,2),'k.')           % line plot
axis equal 

filename = strcat('displacement_',num2str(plotFiles(2)),'.txt');
disp = csvread(filename);


subplot(2,3,2)       % add first plot in 2 x 2 grid
plot(disp(:,1),disp(:,2),'k.')           % line plot
 

axis equal 



filename = strcat('displacement_',num2str(plotFiles(4)),'.txt');
disp = csvread(filename);


subplot(2,3,3)       % add first plot in 2 x 2 grid
plot(disp(:,1),disp(:,2),'k.')           % line plot
axis equal


filename = strcat('displacement_',num2str(plotFiles(7)),'.txt');
disp = csvread(filename);


subplot(2,3,4)       % add first plot in 2 x 2 grid
plot(disp(:,1),disp(:,2),'k.')           % line plot
axis equal



filename = strcat('displacement_',num2str(plotFiles(8)),'.txt');
disp = csvread(filename);


subplot(2,3,5)       % add first plot in 2 x 2 grid
plot(disp(:,1),disp(:,2),'k.')           % line plot
axis equal

filename = strcat('displacement_',num2str(plotFiles(9)),'.txt');
disp = csvread(filename);


subplot(2,3,6)       % add first plot in 2 x 2 grid
plot(disp(:,1),disp(:,2),'k.')           % line plot
axis equal
%plot(d_nodes1(:,1),d_nodes1(:,2),'g+','markersize',8);
%plot(d_nodes2(:,1),d_nodes2(:,2),'y+','markersize',8);
%plot(p_nodes(:,1),p_nodes(:,2),'b+','markersize',8);
%axis off

%legend('nodes','ux = 0', 'uy = 0','Prescribed Displacement')



sizeA = [2 inf];
formatSpec = '%f %f';
fileID = fopen('./../../build/bin/beamUL/loadDisp.txt','r');
A = fscanf(fileID,formatSpec,sizeA);

A = A'

figure

plot(A(:,1),A(:,2),'k-');

sizeA = [2 inf];
formatSpec = '%f %f';
fileID = fopen('./exactSol.txt','r');
B = fscanf(fileID,formatSpec,sizeA);

B = B'

hold on 
plot(B(:,1),B(:,2),'r.');
