close all
clear all 

set(0,'DefaultTextFontname', 'latex')
set(0,'DefaultAxesFontName', 'latex')

path = './../../build/bin/beamUL/Cells';

addpath(path)

files = dir('./../../build/bin/beamUL/Cells/*.txt');
for file = files'
    [poly, ~, ~ ] = read_polygon(file.name);
    hold on
    fill(poly(:,1),poly(:,2),rand(1,3));
    % Do some stuff
end


ylim([0,1])
axis equal



