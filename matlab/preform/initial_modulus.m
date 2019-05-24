clear all
close all

eta = 0.001;
alpha = 0.1553;
Ns = 1.81e17;
temperature = 97.45+273.15;

params = [alpha;Ns;eta;temperature];

lambda_1 = 1.0001;
lambda_2 = 1.0002;

F = @(lambda) [lambda, 0, 0 ; 0, lambda  0 ; 0, 0, 1/(lambda*lambda)];
lambda_vec = @(lambda) [lambda,lambda,1/lambda^2];
[stress] = edwards_vilgis(params,lambda_vec(1.0001));

h = 1e-3;

% edwards vilgis model
for i = 1:10
   lambda = lambda_vec(1+i*h); 
   [stress] = edwards_vilgis(params,lambda);
   
   f(i) = stress(1,1);

end

first_derivative = diff(f)/h;


mu = first_derivative(1)/3;

nu = 0.495;

kappa = 2*mu * ( 1+nu)/(3*(1-2*nu));

%%  edwards-vilgis model
function [output_stress] = edwards_vilgis(params,lamda)

    %  Material model parameters  
    alpha = params(1);
    Ns = params(2);   
    Eta = params(3);       
    temperature = params(4);
    
    % boltzmann constant
    kB = 1.38e-17;
    Nc = 0;

    I1 = 0;
    for k = 1:1:3
        s(k) = lamda(k)^2;
        I1 = I1 + s(k);
    end
    ta = Ns * kB * temperature;
    tb = Nc * kB * temperature;
    tc = 1 + Eta;
    td = alpha^2;
    te = 1 - td;
    tf = 1 - td * I1;
    tg = 1 - 2 * td;
    th = 0;
    for k = 1:1:3
        th = th + s(k) / (1 + Eta * s(k));
    end
 
    for k = 1:1:3
    %   sigCL part
        if (tf ~= 0)
            tempa = tg / tf;
            tempb = td * te * I1 / tf^2;
        end
        tempc = s(k) * tb;
        sigCL = tempc * (tempa + tempb);
    %   sigSL part
        tempa = s(k) * ta;
        tempb = tc * te;
        tempc = 1 + Eta * s(k);
        if (tf ~= 0)
            tempe = tempb * td * th / tf^2;
            tempg = td / tf;
            if (tempc ~= 0)
                tempd = tempb / tf / tempc^2;
            end
        end
        if (tempc ~= 0)
            tempf = Eta / tempc;
        end
        sigSL = tempa * (tempd + tempe + tempf - tempg);
        stress(k) = (sigCL + sigSL);
    end
    %   turn to deviatoric stress
    for k = 1:1:2
        output_stress(k,k) = (stress(k) - stress(3))/1000;
    end
    output_stress(3,3) = 0;
    
end