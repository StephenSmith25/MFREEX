clear all

close all 
path = './../../build/bin/cylinder/Displacement';

addpath(path)
displacementdir = path ;
d = dir(displacementdir);
d1 = dir([displacementdir,'*.txt']);
numFiles = size(d,1) -3 ;

plotFiles = ceil(linspace(1,numFiles,10));
material_points = csvread('./../../build/bin/cylinder/materialpoints.csv');


filename = strcat(path,'/displacement_',num2str(plotFiles(1)),'.txt');
disp = csvread(filename);
figure
subplot(1,3,1)       % add first plot in 2 x 2 grid
plot(disp(:,1),disp(:,2),'k.')           % line plot
xlim([0,30])
ylim([0,30])


axis equal 

filename = strcat(path,'/displacement_',num2str(plotFiles(5)),'.txt');
disp = csvread(filename);


subplot(1,3,2)       % add first plot in 2 x 2 grid
plot(disp(:,1),disp(:,2),'k.')           % line plot
 
xlim([0,30])
ylim([0,30])

axis equal 


filename = strcat(path,'/displacement_',num2str(plotFiles(10)),'.txt');
disp = csvread(filename);


subplot(1,3,3)       % add first plot in 2 x 2 grid
plot(disp(:,1),disp(:,2),'k.')           % line plot
axis equal
hold on
plot(material_points(:,1),material_points(:,2),'r*');
xlim([0,40])
ylim([0,30])

axis equal

%plot(d_nodes1(:,1),d_nodes1(:,2),'g+','markersize',8);
%plot(d_nodes2(:,1),d_nodes2(:,2),'y+','markersize',8);
%plot(p_nodes(:,1),p_nodes(:,2),'b+','markersize',bu8);
%axis off

% %legend('nodes','ux = 0', 'uy = 0','Prescribed Displacement')



sizeA = [2 inf];
formatSpec = '%f %f';
fileID = fopen('./../../build/bin/cylinder/loadDisp.txt','r');
A = fscanf(fileID,formatSpec,sizeA);

A = A';

figure

plot(A(1:1:end,1),A(1:1:end,2),'kx','markersize',6);
% 
sizeA = [2 inf];
formatSpec = '%f %f';
fileID = fopen('./exactSol.txt','r');

B = fscanf(fileID,formatSpec,sizeA);

B = B';

hold on 
plot(B(:,1),B(:,2),'r-');
