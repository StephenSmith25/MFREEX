clear all
close all 

close all 
path = './../../build/bin/cylinder/Displacement';

addpath(path)
displacementdir = path ;
d = dir(displacementdir);
d1 = dir([displacementdir,'*.txt']);
numFiles = size(d,1) -2 ;

num_loops = 50;




plotFiles = ceil(linspace(1,numFiles,num_loops));


F(num_loops) = struct('cdata',[],'colormap',[]);

for i = 1:num_loops
    
    filename = strcat(path,'/displacement_',num2str(plotFiles(i)),'.txt');
    
    disp = csvread(filename);
    
    plot(disp(:,1),disp(:,2),'r.')
  
    set(gcf,'color','w');
    
    xlim([0,50])
ylim([0,50])
    
    
    axis off
    axis equal
    drawnow
    F(i) = getframe(gcf);

    
end