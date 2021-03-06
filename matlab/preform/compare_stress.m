

tmax = 10;
D = 1;
eta = 3;
E = 1+eta;
alpha = 1;
lambda = [1,1,1];
lambda1 = 1;
lambda2 = 1;
lambda3 = 1;


zeta = 1-alpha^2*(lambda1 + lambda2 + lambda3);
xi = lambda1/(1+eta*lambda1) + lambda2/(1+eta*lambda2) + lambda3/(1+eta*lambda3) ;



sigma = 2*D*alpha * lambda * ( E * xi/zeta^2 - 1/zeta) + D*E *(2*lambda)./(zeta*(1+eta*lambda).^2)