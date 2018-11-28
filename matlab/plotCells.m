close all
clear all 

set(0,'DefaultTextFontname', 'latex')
set(0,'DefaultAxesFontName', 'latex')




% figure
% filename = "cells.txt";
% [lShape, ~, ~ ] = read_polygon(filename);
% fill(lShape(:,1),lShape(:,2),'b');
% xlim([min(lShape(:,1)) - 0.5, max(lShape(:,1)) + 0.5])
% ylim([min(lShape(:,2)) - 0.5, max(lShape(:,2)) + 0.5])
% axis equal
% area = polyarea(lShape(:,1),lShape(:,2))
% 
% 
% 
% 
% figure
% filename = "cells.txt";
% subplot(1,2,1)
% [lShape, ~, ~ ] = read_polygon(filename);
% fill(lShape(:,1),lShape(:,2),'y');
% xlim([min(lShape(:,1)) - 0.5, max(lShape(:,1)) + 0.5])
% ylim([min(lShape(:,2)) - 0.5, max(lShape(:,2)) + 0.5])
% axis equal
% area = polyarea(lShape(:,1),lShape(:,2));
% 
% subplot(1,2,2)
% 
% ix = find ( lShape(:,2) == max(lShape(:,2)));
% 
% insideSurface = lShape(2:ix(1),:);
% insideSurface = [insideSurface ; [0, max(insideSurface(:,2))]];
% 
% 
% hold on
% % 
% filename = "result1.txt";

figure
files = dir('Cells/*.txt');

for file = files'
    [polya, ~, ~ ] = read_polygon(file.name);
    hold on
    fill(polya(:,1),polya(:,2),rand(1,3));
    % Do some stuff
end


 axis equal
% 
% figure
% fill(insideSurface(:,1),insideSurface(:,2),'y')
% axis off 
% axis equal
% 
% xq = linspace(0,max(insideSurface(:,1)),20);
% yq = linspace(0,max(insideSurface(:,2)),20);
% 
% [XQ,YQ] = meshgrid(xq,yq);
% XQ = XQ(:);
% YQ = YQ(:);
% 
% [in,on] = inpolygon(XQ,YQ,insideSurface(:,1),insideSurface(:,2));
% 
% points_in = in
% 
% points = [XQ(points_in),YQ(points_in)];
% 
% hold on
% plot(points(:,1),points(:,2),'k*')
% points = [points ; insideSurface];
% points = unique(points,'rows')
% TRI = delaunay(points(:,1),points(:,2));
% 
% 
% hold on
% triplot(TRI,points(:,1),points(:,2));

