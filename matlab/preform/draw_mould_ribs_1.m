
clear all
close all

ylim = -300;
thickness = 10;
radius = 40;
radius_wall = 13.5;
wall_length = 270;
taper_length =25;

fillet_top = 5;
fillet_bot = 10;

% ribs 
rib_radii = 12;
rib_major_radii = 15;
rib_minor_radii = 5;

first_location = (3/8) * wall_length;

num_node_rib = 18;

second_location = (5/8) * wall_length;
extra_height = 10;

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
mould_coords = [mould_coords ; [radius_wall,0+extra_height]];


%mould_coords = flip(mould_coords);

%export mould to file



% drwa mould

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







%inside wall
mould_coords = [mould_coords ; [radius_wall+thickness,0+extra_height]];
%inside taper
mould_coords = [mould_coords ; [radius_wall+thickness,ylim + wall_length+taper_length]];


% top corner
mould_coords = [mould_coords ; [radius+thickness,ylim+wall_length]];

%bottom corner
mould_coords = [mould_coords ; [radius+thickness,ylim-thickness]];

%bottom
mould_coords = [mould_coords ; [0,ylim-thickness]];



% mould_coords = flip(mould_coords);

hold on
fill(mould_coords(:,1),mould_coords(:,2),'b');
axis equal 


for i = 1:length(mould_coords)
   segments(i,:) = [i,i+1];
   
   if ( i == length(mould_coords))
          segments(i,:) = [i,1];

       
   end
   
end

% write the triangle file

fileID = fopen('mould.poly','w');

fprintf(fileID,'%d %d %d \n',length(mould_coords),2,0)

for i = 1:length(mould_coords)
   fprintf(fileID,'%d %f %f\n',i,mould_coords(i,1),mould_coords(i,2))
    
end
fprintf(fileID,'%d \n',length(segments))

for i = 1:length(mould_coords)
   fprintf(fileID,'%d %d %d\n',i,segments(i,1),segments(i,2))
    
end
fprintf(fileID,'%d \n',0)

fclose(fileID);






% write a vtk file 
fileID = fopen('mould.1.node','r');

line= fgetl(fileID);
line = split(line,' ');
numnodes = str2num(line{1});

for i = 1:numnodes
    line = fgetl(fileID);
    line = split(line,' ');

    
    if ( i < 10 )
    nodes(i,:) = [str2num(line{8}),str2num(line{10})];
    
    
    elseif ( i < 100)
     nodes(i,:) = [str2num(line{7}),str2num(line{9})];

    
    else    
     nodes(i,:) = [str2num(line{6}),str2num(line{8})];
      
    end    
end
fileID = fopen('mould.1.ele','r');


tri = [];
line= fgetl(fileID);
line = split(line,' ');
numele = str2num(line{1});

tri = zeros(numele,3);

for i = 1:numele
    line = fgetl(fileID);
    line = split(line,' ');

    count = 1;
    
        
    for k = 1:length(line)
        
       if ( ~isempty(line{k}))
           
           if ( count > 1)
               
              tri(i,count-1) = str2num(line{k}); 
           end
           
           
           
           count = count+1;
           
       end
       % else move on to next line
    end
    
end

fclose(fileID);
figure
plot(nodes(:,1),nodes(:,2),'r.');

triplot(tri,nodes(:,1),nodes(:,2));
z = zeros(length(nodes),1);

vtkwrite('mould.vtk','polydata','triangle',nodes(:,1),nodes(:,2),z,tri);    




axis off
axis equal











