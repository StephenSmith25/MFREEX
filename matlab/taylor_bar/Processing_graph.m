clear all
close all 

path = './../../build/bin/taylor_bar/Displacement';

addpath(path)
displacementdir = path ;
d = dir(displacementdir);
d1 = dir([displacementdir,'*.txt']);
numFiles = size(d,1) -3 ;

plotFiles = ceil(linspace(1,numFiles,30));
A = [];
t_max = 0.000051*10^6;
figure
for i = 1:length(plotFiles)
    
    
% second plot
filename = strcat(path,'/displacement_',num2str(plotFiles(i)),'.txt');
disp = csvread(filename);

y = max(disp(:,1))*100;


t = (t_max/(length(plotFiles)-1))*(i-1);

A = [A;t,y];




hold on
plot(t,y,'b.');


  
end

xlim([0,t_max])
