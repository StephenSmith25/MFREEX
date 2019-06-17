

close all
clear all

N = csvread('./../../build/bin/biaxial/Stress_11_nominal.txt');


figure


M = csvread('./BiaxialData/simul_sr4_t85_3.5x3.5.csv');
nominalStrain_E = M(:,3);
trueStrain_E = log( 1 + nominalStrain_E);
trueStress_E = M(:,4);
nominalStress_E = trueStress_E./trueStrain_E;
plot(nominalStrain_E,trueStress_E,'b+')

% 
M = csvread('./BiaxialData/simul_sr16_t105_4x4.txt');
nominalStrain_E = M(:,3);
trueStrain_E = M(:,2);
trueStress_E = M(:,3);
nominalStrain_E = exp(trueStrain_E) -1 ;
nominalStress = trueStress_E./exp(trueStrain_E) ; 
plot(nominalStrain_E,trueStress_E,'b-')





nominalStrain_S = N(:,1);
trueStress_S = N(:,2);

hold on

plot(nominalStrain_S(1:10:end),trueStress_S(1:10:end),'k')


%set axis
set(gca, 'FontName', 'cmr12')
% set x tics and y tics
set(gca,...
'Units','normalized',...
'FontWeight','normal',...
'FontSize',14,... % size ofiguref numbers on axis
'FontName','cmr14') % font name
set(gca,'TickLabelInterpreter', 'latex');

% Y label
ylabel({'True Stress (MPa) '},...
'interpreter','latex',...
'FontSize',14,... % font size
'FontName','cmr14')
% X label
xlabel('Nominal Strain',...
'interpreter','latex',...
'FontWeight','normal',...
'FontSize',14,... % font size
'FontName','cmr14')
% legend

legend({'Experimental','Buckley Model'},... % { 'legend1', 'legend2',...}
'interpreter','latex',...
'FontSize',12,...
'FontName','cmr14',...
'Location','NorthWest')


set(gcf, 'Color', 'w');

print -dpng2 stress_strain.png



