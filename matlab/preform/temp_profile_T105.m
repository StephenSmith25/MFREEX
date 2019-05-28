function [temperatures] = get_node_temperature(OIL_TEMPERATURE,COOLING_TIME,x,y)


close all
clear all

% if less lose the through thickness temperature profile



% if in neck region we have to modify this problem 
[X,Y,TEMPERATURE] = temperature_profile(105,18);


Y = Y + 77.44;




h2 = figure(2);
h2 = pcolor(X,Y,TEMPERATURE);
set(h2, 'LineStyle','none')
%axis off;
axis([0,80,0,80]);
%%caxis([101 105])
colorbar;



end