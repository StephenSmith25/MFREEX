%% generate a random profile, then revolve it and display
clear all
close all 
% random points
pts = round(5*rand(10,2)+5);
pts = unique(pts, 'rows');
 
% turn it into a bounding profile
shp = alphaShape(pts(:,1), pts(:,2), 3.5);
[bf, bp] = boundaryFacets(shp);
 

num_revolves = 8;


% sweep the profile
[tri, xyz] = bfRevolve(bf, bp, num_revolves);
 

% display the surface
figure


for i = 1:length(bf) 
    for j = 1:num_revolves - 4
        tri_mod((i-1)*num_revolves/2 + j,:) = tri((i-1)*num_revolves + j,:);  
    end
    
    end

trisurf(tri_mod, xyz(:,1),xyz(:,2),xyz(:,3));

axis([-10 10 -10 10 0 10])