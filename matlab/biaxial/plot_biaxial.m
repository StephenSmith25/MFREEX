

close all
clear all 
M = csvread('./BiaxialData/simul_sr4_t105_4.5x4.5.txt');
N = csvread('./../../build/bin/biaxial/Stress_11_nominal.txt');


figure


M = csvread('./BiaxialData/simul_sr4_t85_3.5x3.5.csv');
% Experiemntal
nominalStrain_E = M(:,3);
trueStrain_E = log( 1 + nominalStrain_E);
trueStress_E = M(:,4);
nominalStress_E = trueStress_E./trueStrain_E;
plot(nominalStrain_E,trueStress_E,'b+')


% M = csvread('./BiaxialData/simul_sr16_t105_4x4.txt');
% Experiemntal
% nominalStrain_E = M(:,3);
% trueStrain_E = M(:,2);
% trueStress_E = M(:,3);
% nominalStrain_E = exp(trueStrain_E) -1 ;
% nominalStress = trueStress./exp(trueStrain) ; 
% plot(nominalStrain_E,trueStress_E,'b+')





nominalStrain_S = N(:,1);
trueStress_S = N(:,2);

hold on

plot(nominalStrain_S,trueStress_S,'k-')
