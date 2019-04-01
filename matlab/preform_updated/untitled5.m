clear all
close all

figure

x = linspace(0,25);
y = sin(x/2);
yyaxis left
plot(x,y);

ylabel('left')

r = x.^2/2;
yyaxis right
plot(x,r);
% LABEL FOR RIGHT AXIS
% CHANGE POSITON BETWEEN 1 AND 1.1 ISH
y = ylabel('right')
set(y, 'Units', 'Normalized', 'Position', [1.1, 0.5, 0]);


% LEGEND
 columnlegend(2, {'Test1','Test2'}, 'location', 'northwest'); 