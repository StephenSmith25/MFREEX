close all
clear all

lambda = 1;
rad = 1;
M_start = [ 1/rad^2 0 ; 0 1/rad^2];
F = [1+lambda 0 ; 0 1+0.5*lambda];
C = F'*F;


M_end = [M_start(1,1)*(1/C(1,1)) 0;
0 M_start(2,2)*(1/C(2,2))];


a = sqrt(1/M_end(1,1));
b = sqrt(1/M_end(2,2));


%

k = 1.5;
rad = 1;
M_start = [ 1/rad^2 0 ; 0 1/rad^2];
F = [1 k ; 0 1];
C = F'*F;


M_end = [M_start(1,1)*(1/C(1,1)) 0;
0 M_start(2,2)*(1/C(2,2))];


a = sqrt(1/M_end(1,1));
b = sqrt(1/M_end(2,2));


k_total =2; 
k = 0.5*k_total;
rad = 1;
M_start = [ 1/rad^2 0 ; 0 1/rad^2];
F = [1 k ; 0 1];

C = F'*F;

M_end = [M_start(1,1)*(1/C(1,1)) 0;
0 M_start(2,2)*(1/C(2,2))];


a_1 = sqrt(1/M_end(1,1));
b_1 = sqrt(1/M_end(2,2));






U = sqrtm(C);
[V,D] = eig(U);


N_1 = [V(1,1),V(1,2)]'
N_2 = [V(2,1),V(2,2)]';


M = [1 0 ; 0 1];


N_1 = [0.8507,0.5257]';
N_2 = [-0.5257,0.8507]';
M_test = (1/D(2,2)^2)*N_1*N_1' + (1/D(1,1)^2)*N_2*N_2'



a = D(1,1)*rad;
b = D(2,2)*rad;



theta=acosd(V(1,2));

theta = -theta;



theta = 31.7175;
R = [cosd(theta) -sind(theta) ; sind(theta) cosd(theta)];

S = [1/D(2,2)^2 0 ; 0 1/D(1,1)^2];


M_elli = R*S*M*R'

alpha = linspace(0,360,50);
theta = 31.7175;

x = b.*cosd(alpha).*cosd(theta) - a.*sind(alpha).*sind(theta);
y = b.*cosd(alpha).*sind(theta) + a.*sind(alpha).*cosd(theta);

plot(x,y,'b.');

hold on

rectangle('Position',[-2 -1.5 4 3])


function [output] = draw_general_ellipse_alt(M,theta,Cx,Cy)

    t = linspace(0,360,50);
    
    MI = [M(1),M(2); M(3), M(4)];
    
    [V,D] = eig(MI);
    
    
    a = 1.00/sqrt(D(1,1));
    b = 1.00/sqrt(D(2,2));
    x = a*cosd(t)*cos(theta) - b*sind(t)*sin(theta);
    y = a*cosd(t)*sin(theta) + b*sind(t)*cos(theta);
    X = [x;y];
    
        
   
    hold on
    plot(x+Cx,y+Cy,'r-');
    
    output = [x'+Cx,y'+Cy];
    
    
    


end
