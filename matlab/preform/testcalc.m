
%(18.63+20)-(20*cosd(theta)) =18.44 



clear all
close all

Cv = 67;
Tinf = 342.68;
T = 273+103;
star_T = 358.15;


Tinf = linspace(320,350,100);
alpha_s = exp(Cv./(T - Tinf) - Cv./(star_T - Tinf));

plot(Tinf,alpha_s);



Vs = 2.814e-3;
Vp = 0.526e-3;
R = 8.314;

sigma_m = linspace(0,1e8,1000);

alpha_sig = (Vs/(2*R*T))*exp(-Vp.*sigma_m/(R*T))/sinh(Vs/(2*R*T));


figure
plot(sigma_m,alpha_sig);








