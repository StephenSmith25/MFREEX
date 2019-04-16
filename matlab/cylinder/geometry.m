clear all

close all 
path = './../../build/bin/cylinder/';

figure



% read nodes
nodes = csvread(strcat(path,'Geometry/nodes.csv'));


plot(nodes(:,1),nodes(:,2),'b.');
