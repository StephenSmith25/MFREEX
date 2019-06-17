close all
clear all





sr = linspace(1,110,64);


N = 0.334;
K = 1*(1.6322e6);


gamma = K.*sr.^(N-1);


plot(sr,gamma);


hold on
N = 0.4;
K = 1*(1.6322e6);


gamma = K.*sr.^(N-1);


plot(sr,gamma);
legend('0.334','0.5')

% 
% sigma_m = linspace(-1e7,1e7,100)
% a = exp(-sigma_m*5.262e-4/(283*8413));
% 
% plot(sigma_m,a);



% 
% 
% C1 = -1.125e-2;
% C2 = 4.5850;
% beta = 0.97884;
% k = -2.5322e-2;
% b = 1.1552e1;
% temperature = 97.45+273.15;
% 
% sr = linspace(1,6,64);
% shifted_temperature = temperature .* 10.^(C1.*(sr-1)./(C2+sr-1));
% critLambda = k.*shifted_temperature + b;
% 
% 
% plot(sr,critLambda);