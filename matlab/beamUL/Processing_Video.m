close all
clear all

path = './../../build/bin/beamUL/Displacement';
displacementdir = path;
d = dir(displacementdir);
d1 = dir([displacementdir,'*.txt']);
numFiles = size(d,1) - 2;






figure

set(gcf,'color','w')
v = VideoWriter('Displacement.avi');
v.Quality = 95;
open(v)
ax = gca();

xlim([0,20])
ylim([-1,20])
axis equal


for i = 1:1:numFiles
filename = strcat(path,'/displacement_',num2str(i),'.txt');
disp = csvread(filename);
p = plot(disp(:,1),disp(:,2),'k.')   ;        % line plot


xlim([0,20])
ylim([-1,20])
 
drawnow();
writeVideo(v,getframe(ax));

end
close(gca)
close(v)


