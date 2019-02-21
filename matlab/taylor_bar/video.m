close all
path = './../../build/bin/taylor_bar/Displacement/';

addpath(path)
displacementdir = path ;
d = dir(displacementdir);
d1 = dir([displacementdir,'*.txt']);
numFiles = size(d,1) -3 ;


num_loops = 150;

F(num_loops) = struct('cdata',[],'colormap',[]);



plotFiles = ceil(linspace(1,numFiles,num_loops));


figure
for i = 1:num_loops
    
    filename = strcat(path,'displacement_',num2str(plotFiles(i)),'.txt');
    
    disp = csvread(filename,1);

   plot(disp(:,1),disp(:,2),'k.',-disp(:,1),disp(:,2),'k.',[-0.012,0.012],[0,0],'g','markersize',3);

  

    set(gcf,'color','w');
    
    ylim([-0.001,0.0250])
    axis equal
    drawnow
    F(i) = getframe(gcf);

    
end

v = VideoWriter('newfile','MPEG-4');
open(v)
writeVideo(v,F)
close(v)

