clear all
close all


path = './../../build/bin/preform/Displacement/';
pathSR = './../../build/bin/preform/srRod/';
boundaryNodes = csvread('./../../build/bin/preform/boundary.txt');
boundaryNodes = [boundaryNodes;boundaryNodes(end)];

boundaryNodes = [boundaryNodes;boundaryNodes(1)];

addpath(path)
displacementdir = path ;
d = dir(displacementdir);
d1 = dir([displacementdir,'*.csv']);
numFiles = size(d,1) -3 ;

plotFiles = ceil(linspace(1,numFiles,200));

filename = strcat(path,'displacement_',num2str(plotFiles(1)),'.csv');
disp = csvread(filename,1);

% draw mould
Radius_mould = 24;
Height_mould = 30;
thickness_mould = 5;
mould_length = 180;
height = 74.2650;
Radius_out = 11.25;

top_point = 48.920;



mould_nodes = [];
mould_nodes = [mould_nodes ; [Radius_out,height+5]   ];
mould_nodes = [mould_nodes ; [Radius_out,68.92]   ];
mould_nodes = [mould_nodes ; [Radius_out+Radius_mould,top_point]   ];
mould_nodes = [mould_nodes ; [Radius_out+Radius_mould,height-mould_length]   ];
mould_nodes = [mould_nodes ; [0,height-mould_length]   ];


mould_nodes = [mould_nodes ; [0,height-mould_length-thickness_mould]   ];
mould_nodes = [mould_nodes ; [Radius_out+Radius_mould+thickness_mould,height-mould_length-thickness_mould]   ];
mould_nodes = [mould_nodes ; [Radius_out+Radius_mould+thickness_mould,top_point]   ];
mould_nodes = [mould_nodes ; [Radius_out+thickness_mould,68.92]   ];
mould_nodes = [mould_nodes ; [Radius_out+thickness_mould,height+5]   ];
mould_nodes = [mould_nodes ; [Radius_out,height+5]   ];

mould_nodes =zeros(3,2);

num_loops = length(plotFiles);

F(num_loops) = struct('cdata',[],'colormap',[]);

for i = 1:num_loops
    
    plot(i,i+1,'r.');
    filename = strcat(path,'displacement_',num2str(plotFiles(i)),'.csv');
    
    disp = csvread(filename,1);
    filename = strcat(pathSR,'srRod_',num2str(plotFiles(i)),'.csv');
disp_1 = csvread(filename,1);

   fill(disp(boundaryNodes,1),disp(boundaryNodes,2),'r',-disp(boundaryNodes,1),disp(boundaryNodes,2),'r', ....
      mould_nodes(:,1),mould_nodes(:,2),[0.5,0.5,0.5],-mould_nodes(:,1),mould_nodes(:,2),[0.5,0.5,0.5], ....
      disp_1(:,1),disp_1(:,2),'g',-disp_1(:,1),disp_1(:,2),'g')   ;

  

    set(gcf,'color','w');
    
   ylim([-250,0])
   xlim([-40,40])
    axis off
    axis equal
    drawnow
    F(i) = getframe(gcf);

    
end

v = VideoWriter('newfile.avi','Motion JPEG AVI');
v.Quality = 100;
open(v)
writeVideo(v,F)
close(v)

