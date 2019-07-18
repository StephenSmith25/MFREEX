clear all
close all 



path = './../../build/bin/rigid_punch/';

figure



% read nodes
nodes = csvread(strcat(path,'nodes.csv'));


plot(nodes(:,1),nodes(:,2),'r.');


% read domain elements
triangles = csvread(strcat(path,'triangles.csv')) ;
hold on
triplot(triangles,nodes(:,1),nodes(:,2));
% read nodesets

axis equal

