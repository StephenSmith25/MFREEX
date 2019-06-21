function [mould_coords] = draw_mould_ribs()


ylim = -300;
thickness = 5;
radius = 40;
radius_wall = 12.3;
wall_length = 270;
taper_length = 20;

fillet_top = 5
fillet_bot = 10;

% ribs 
rib_radii = 12;
rib_major_radii = 15;
rib_minor_radii = 5;

first_location = (3/8) * wall_length;

num_node_rib = 18;

second_location = (5/8) * wall_length;


% draw inside shape


mould_coords = [];

% bottom
mould_coords = [mould_coords ; [0,ylim]];

%bottom corner
mould_coords = [mould_coords ; [radius,ylim]];

%

mould_coords = [mould_coords ; [radius,ylim+first_location-rib_major_radii]];

theta = linspace(-120,-240,num_node_rib);
for i = 1:num_node_rib
   mould_coords = [mould_coords ; [radius + rib_minor_radii*cosd(theta(i)),ylim+first_location+rib_major_radii*sind(theta(i))]]; 
    
end
mould_coords = [mould_coords ; [radius,ylim+first_location+rib_major_radii]];
mould_coords = [mould_coords ; [radius,ylim+second_location-rib_major_radii]];

for i = 1:num_node_rib
   mould_coords = [mould_coords ; [radius + rib_minor_radii*cosd(theta(i)),ylim+second_location+rib_major_radii*sind(theta(i))]]; 
    
end
mould_coords = [mould_coords ; [radius,ylim+rib_major_radii+second_location]];


% top corner
mould_coords = [mould_coords ; [radius,ylim+wall_length]];

%inside taper
mould_coords = [mould_coords ; [radius_wall,ylim + wall_length+taper_length]];

%inside wall
mould_coords = [mould_coords ; [radius_wall,0]];


mould_coords = flip(mould_coords);

%export mould to file



% drwa mould

figure
hold on
plot(mould_coords(:,1),mould_coords(:,2),'linewidth',3);
axis off
axis equal 


path = './../../build/bin/preform/Displacement/';
addpath(path)


displacementdir = path ;
d = dir(displacementdir);
d1 = dir([displacementdir,'*.csv']);
numFiles = size(d,1) -3 ;

plotFiles = ceil(linspace(1,numFiles,10));

plot_point =138; %%2241
filename = strcat(path,'displacement_',num2str(plotFiles(1)),'.csv');
disp = csvread(filename,1);


boundaryNodes = csvread('./../../build/bin/preform/boundary.txt');

boundaryNodes = [boundaryNodes;boundaryNodes(1)];



hold on
plot(disp(:,1),disp(:,2),'k.','markersize',4)           % line plot
hold on
ymax = 0;
hold on





%segments
dlmwrite('../../problems/preform/mould.coords',length(mould_coords),'delimiter',' ')
dlmwrite('../../problems/preform/mould.coords',mould_coords,'-append',.........
    'delimiter',' ')









end







