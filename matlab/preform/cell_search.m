close all
clear all
clc


% generate random points in in a domain
L = 2;
W = 10;
num_width = 20;
num_length = 5;
x = linspace(0,W,num_width);
y = linspace(0,L,num_length);

[X,Y] = meshgrid(x,y);


nodes(:,1) = X(:);
nodes(:,2) = Y(:);

plot(nodes(:,1),nodes(:,2),'r.');
axis off
axis equal

% generate materail point


% find domain of influence of material points


% find max domain of influence


% generate cells