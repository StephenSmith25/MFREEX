clear all
close all 



path = './../../build/bin/cylinder/';

figure



% read nodes
nodes = csvread(strcat(path,'Geometry/nodes.csv'));


plot(nodes(:,1),nodes(:,2),'r.');


% read domain elements
triangles = csvread(strcat(path,'Geometry/triangles.csv')) + 1;
hold on
triplot(triangles,nodes(:,1),nodes(:,2));
% read nodesets



% read SIDESETS
pressure_set = (csvread(strcat(path,'Geometry/pressureset.csv')) + 1)';
% for i = 1:length(pressure_set)
%     hold on
%    plot(nodes(pressure_set(i,:),1),nodes(pressure_set(i,:),2),'k-','linewidth',5); 
%     
% end
plot(nodes(pressure_set(:),1),nodes(pressure_set(:),2),'k-','linewidth',5); 

node_set_1 = csvread(strcat(path,'Geometry/nodeset_1.csv')) + 1;
plot(nodes(node_set_1(:,:),1),nodes(node_set_1(:,:),2),'r-','linewidth',5); 


node_set_2 = csvread(strcat(path,'Geometry/nodeset_2.csv')) + 1;
plot(nodes(node_set_2(:,:),1),nodes(node_set_2(:,:),2),'m-','linewidth',5); 


axis equal 
