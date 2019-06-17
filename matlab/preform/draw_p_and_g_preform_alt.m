% Draw Planar-Straight-Line-Graph of a preform
clc;
clear all;
close all;

%% RADII






[nodes,segments,X_GRID,Y_GRID,TEMP_GRID] = get_node_temperature(105,16);
%Vq = interp2(X_GRID,Y_GRID,TEMP_GRID,-30,10)

Y_GRID(:) = Y_GRID(:) - max(Y_GRID(:));
Y_GRID(:) = Y_GRID(:) ./ min(Y_GRID(:));

Y_GRID(:) = Y_GRID(:) - max(Y_GRID(:));

F = TriScatteredInterp(X_GRID(:),Y_GRID(:),TEMP_GRID(:));



nodes(:,2) = nodes(:,2) - max(nodes(:,2));

Y = nodes(:,2)./min(nodes(:,2));
[~,ix] = find(segments(:,3) == 2);
inner_X = nodes(5:length(ix)+4,:);
count = length(inner_X);

[~,ix] = find(segments(:,3) == 6);
outer_X = nodes(count+4:count+length(ix)+7,:);
outer_X = [outer_X ;nodes(1:4,:) ];

figure


plot(inner_X(:,2),inner_X(:,3),'b.');

hold on
plot(outer_X(:,2),outer_X(:,3),'r.');
%plot(nodes(:,1),nodes(:,2),'b.');


% 
% axis off
% axis equal
% 
% %boundary nodes
% 
% 
tri = csvread('./../../build/bin/preform/triangles.csv');




path = './../../build/bin/preform/Displacement/';
pathSR = './../../build/bin/preform/srRod/';
addpath(path)


displacementdir = path ;
d = dir(displacementdir);
d1 = dir([displacementdir,'*.csv']);
numFiles = size(d,1) -3 ;

plotFiles = ceil(linspace(1,numFiles,10));

plot_point =138; %%2241
filename = strcat(path,'displacement_',num2str(plotFiles(1)),'.csv');
nodes = csvread(filename,1);

figure
triplot(tri,nodes(:,1),nodes(:,2));


%F = TriScatteredInterp(X_GRID(:),Y_GRID(:),TEMP_GRID(:));

boundaryNodes = csvread('./../../build/bin/preform/boundary.txt');

plot(nodes(boundaryNodes,1),nodes(boundaryNodes,2),'r-');


axis off
axis equal

% figure
% p = patch('Faces',tri,'Vertices',nodes(:,1:2),'FaceVertexCData',nodes(:,3));
% p.FaceColor = 'flat';
% colorbar
% axis equal
% figure 
% 
% plot(nodes(:,1),nodes(:,2),'k.');
% hold on
% %plot(nodes(boundaryNodes(:,1),1),nodes(boundaryNodes(:,1),2),'b-')
% axis off
% axis equal
% hold on
% C = [];
% 
% 
%