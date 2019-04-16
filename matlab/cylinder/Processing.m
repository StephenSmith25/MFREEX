clear all

close all 
path = './../../build/bin/cylinder/Displacement';

addpath(path)
displacementdir = path ;
d = dir(displacementdir);
d1 = dir([displacementdir,'*.txt']);
numFiles = size(d,1) -3 ;

plotFiles = ceil(linspace(1,numFiles,10));



figure





%material_points = csvread('./../../build/bin/cylinder/materialpoints.csv');

filename = strcat(path,'/displacement_',num2str(plotFiles(10)),'.txt');
disp = csvread(filename);
plot(disp(:,1),disp(:,2),'k.')           % line plot
axis equal
hold on
%plot(material_points(:,1),material_points(:,2),'r.');
xlim([0,40])
ylim([0,30])




%material_points = csvread('./../../build/bin/cylinder/materialpoints_0.csv');

filename = strcat(path,'/displacement_',num2str(plotFiles(1)),'.txt');
disp = csvread(filename);
plot(disp(:,1),disp(:,2),'k.')           % line plot
axis equal
hold on
%plot(material_points(:,1),material_points(:,2),'r.');
xlim([0,40])
ylim([0,30])



saveas(gcf,'Displacement_cylinder','epsc')



sizeA = [2 inf];
formatSpec = '%f %f';
fileID = fopen('./../../build/bin/cylinder/loadDisp.txt','r');
A = fscanf(fileID,formatSpec,sizeA);

A = A';

figure



plot(A(1:10:end,1),A(1:10:end,2),'kx','markersize',6);
% 
sizeA = [2 inf];
formatSpec = '%f %f';
fileID = fopen('./exactSol.txt','r');

B = fscanf(fileID,formatSpec,sizeA);

B = B';

hold on 
plot(B(:,1),B(:,2)/1000,'r-');
xlabel('dr')
ylabel('Pressure')
legend('Meshfree','Exact','Location','northwest');

saveas(gcf,'Solution_cylinder','epsc')
