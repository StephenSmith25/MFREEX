close all
clear all 


displacementdir = './Displacement';
d = dir(displacementdir);
d1 = dir([displacementdir,'*.txt']);
numFiles = size(d,1) - 2;

plotFiles = ceil(linspace(1,numFiles,10));


filename = strcat('displacement_',num2str(plotFiles(1)),'.txt');
disp = csvread(filename);
figure
subplot(1,3,1)       % add first plot in 2 x 2 grid
plot(disp(:,1),disp(:,2),'k.')           % line plot
axis equal 

filename = strcat('displacement_',num2str(plotFiles(5)),'.txt');
disp = csvread(filename);


subplot(1,3,2)       % add first plot in 2 x 2 grid
plot(disp(:,1),disp(:,2),'k.')           % line plot
 

axis equal 



filename = strcat('displacement_',num2str(plotFiles(10)),'.txt');
disp = csvread(filename);


subplot(1,3,3)       % add first plot in 2 x 2 grid
plot(disp(:,1),disp(:,2),'k.')           % line plot
axis equal



%plot(d_nodes1(:,1),d_nodes1(:,2),'g+','markersize',8);
%plot(d_nodes2(:,1),d_nodes2(:,2),'y+','markersize',8);
%plot(p_nodes(:,1),p_nodes(:,2),'b+','markersize',8);
%axis off

%legend('nodes','ux = 0', 'uy = 0','Prescribed Displacement')


