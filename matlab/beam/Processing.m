clear all

close all 
path = './../../build/bin/beam/Displacement';



cells = csvread('./../../build/bin/beam/search_cells.csv');

figure
for k = 1:length(cells)
hold on
    i = k;
rectangle('Position',[cells(i,1), cells(i,3), cells(i,2)-cells(i,1), cells(i,4)-cells(i,3)])
end

axis equal

nodes = csvread('./../../build/bin/beam/nodes.csv');

hold on
plot(nodes(:,1),nodes(:,2),'r.');

boundary_nodes = csvread('./../../build/bin/beam/boundary.txt');

boundary_nodes = [boundary_nodes;boundary_nodes(1)];

hold on
plot(nodes(boundary_nodes,1),nodes(boundary_nodes,2),'r-');

addpath(path)
displacementdir = path ;
d = dir(displacementdir);
d1 = dir([displacementdir,'*.txt']);
numFiles = size(d,1) -3 ;

plotFiles = ceil(linspace(1,numFiles,10));


filename = strcat(path,'/displacement_',num2str(plotFiles(1)),'.txt');
disp = csvread(filename);
figure
subplot(1,3,1)       % add first plot in 2 x 2 grid
plot(disp(:,1),disp(:,2),'k.')           % line plot
axis equal 
hold on
plot(disp(boundary_nodes,1),disp(boundary_nodes,2),'r-');
filename = strcat(path,'/displacement_',num2str(plotFiles(4)),'.txt');
disp = csvread(filename);
disp_temp = [disp(:,1)-20,disp(:,2)];
disp_temp = abs(disp_temp);

mid_boundary = disp_temp(boundary_nodes,:);

subplot(1,3,2)       % add first plot in 2 x 2 grid
plot(disp(:,1),disp(:,2),'k.')           % line plot
hold on
plot(disp(boundary_nodes,1),disp(boundary_nodes,2),'r-');



axis equal 



filename = strcat(path,'/displacement_',num2str(plotFiles(10)),'.txt');
disp = csvread(filename);


subplot(1,3,3)       % add first plot in 2 x 2 grid
plot(disp(:,1),disp(:,2),'k.')           % line plot
axis equal
hold on
plot(disp(boundary_nodes,1),disp(boundary_nodes,2),'r-');


final_boundary = [disp(boundary_nodes,1),disp(boundary_nodes,2)];

sizeA = [2 inf];
formatSpec = '%f %f';
fileID = fopen('./../../build/bin/beam/loadDisp.txt','r');
A = fscanf(fileID,formatSpec,sizeA);

A = A';


A_mod = A(1:100:end,:)
A_mod = [A_mod;A(end,:)];
figure

plot(A_mod(:,1),A_mod(:,2),'bo');

sizeA = [2 inf];
formatSpec = '%f %f';
fileID = fopen('./exactSol.txt','r');
B = fscanf(fileID,formatSpec,sizeA);

B = B';

hold on 
plot(B(:,1),B(:,2),'k-','markersize',6);






ix = find(A(:,1) < 0.8);
P_num = A(ix,2);
u_num = A(ix,1);


P_exact = interp1(B(:,1),B(:,2),u_num(:,1));

figure
plot(u_num,P_exact,'r.');
hold on
plot(u_num,P_num,'bo');

R2 = calculateR2(P_num,P_exact)

function R2 = calculateR2(z,z_est)
% calcuateR2 Cacluate R-squared
% R2 = calcuateR2(z,z_est) takes two inputs - The observed data x and its
% estimation z_est (may be from a regression or other model), and then
% compute the R-squared value a measure of goodness of fit. R2 = 0
% corresponds to the worst fit, whereas R2 = 1 corresponds to the best fit.
% 
% Copyright @ Md Shoaibur Rahman (shaoibur@bcm.edu)
r = z-z_est;
normr = norm(r);
SSE = normr.^2;
SST = norm(z-mean(z))^2;
R2 = 1 - SSE/SST;
end




