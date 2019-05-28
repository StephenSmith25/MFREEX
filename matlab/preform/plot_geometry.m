
close all
clear all

nodes = csvread('./../../build/bin/preform/nodes.csv');

tri = csvread('./../../build/bin/preform/triangles.csv');

stress_points = csvread('./../../build/bin/preform/stress_points.csv');

%plot_cells();
hold on
plot(stress_points(:,1),stress_points(:,2),'r.');
hold on
plot(nodes(:,1),nodes(:,2),'b.');
axis equal
hold on
triplot(tri,nodes(:,1),nodes(:,2));


function [] = plot_cells()
fileID = fopen("./../../build/bin/preform/cells.txt");


tline = fgetl(fileID);
C = strsplit(tline);
num_verticies = str2double(C{1});
dim = str2double(C{2});

verticies = zeros(num_verticies,dim);

for i = 1:num_verticies


    tline = fgetl(fileID);
    C = strsplit(tline);
    
    for k = 1:dim
        verticies(i,k) = str2double(C{k});
    end
    
        
    
end





tline = fgetl(fileID);
C = strsplit(tline);
num_cells = str2double(C{1});

figure
for i = 1:num_cells

    tline = fgetl(fileID);
    C = strsplit(tline);
    
    num_cell_verticies = str2double(C{1});
    
    poly = zeros(num_cell_verticies,2);
    for k = 1:num_cell_verticies
        poly(k,1) = verticies(str2double(C{1+k})+1,1);
        poly(k,2) = verticies(str2double(C{1+k})+1,2);

    end
    area(i) = polyarea(poly(:,1),poly(:,2));
    
    hold on 
    fill(poly(:,1),poly(:,2),rand(1,3));
end
sum(area)
max(area)
axis equal

end
