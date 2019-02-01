close all
clear all

F = [ 1 3 0 ; 0 1 0 ; 0 0 1];


C= F'*F;

[Vec D] = eig(C);

 


lambda_1 = sqrt(D(1,1));
lambda_2 = sqrt(D(2,2));
lambda_3 = sqrt(D(3,3));

Csqr = Vec*D^2*Vec';

i1 = lambda_1 + lambda_2 + lambda_3;
i2 = lambda_1*lambda_2 + lambda_1*lambda_3 + lambda_2*lambda_3;
i3 = lambda_1*lambda_2 * lambda_3;

D1 = i1*i2 - i3;

V = (1/D1)* [ -Csqr + (i1^2 - i2)*C + i1*i3*eye(3)];

