close all
clear all


k_total =0.5; 
k = k_total;
rad = 1;
M_start = [ 1/rad^2 0 ; 0 1/rad^2];
F = [1 k ; k 1];

C = F'*F;


M_end = [M_start(1,1)*(1/C(1,1)) 0;
0 M_start(2,2)*(1/C(2,2))];


a_1 = sqrt(1/M_end(1,1));
b_1 = sqrt(1/M_end(2,2));

% anayltical formula
delta = k^4 -2*k^2 + 1;
tau = 2+2*k^2;

U_anal = 1/(sqrt(tau+2*sqrt(delta))) .* [1+k^2 + sqrt(delta) , 2*k ; 2*k , 1+k^2 + sqrt(delta)]






t = 1;
U_2 = 1/(sqrt(4+t^2)) .* [2 t ; t 2+t^2];
M_p2 = [(4+t^2)/4 0 ; 0 (4+t^2)/(2+t^2)^2]; 
a_p2 = sqrt(1/M_p2(1,1));
b_p2 = sqrt(1/M_p2(2,2));


aa = (2+t^2)/2 + sqrt((2+t^2)^2/4 -1);
bb = (2+t^2)/2 - sqrt((2+t^2)^2/4 -1);

B = [aa 0 ; 0 bb];
A = C;
[Da, eigA] = eig(A); % Find eigenvectors and eigenvalues of A
[Db, eigB] = eig(B); % Find eigenvectors and eigenvalues of B
Q = Db'*Da; % This relation can be derived if we substitute B = Db^T*eigB*Db and A = Da^T*eigA*Da in above equation and noting that eigA and eigB are same because it is a similarity transform
rotation_angle = acosd(0.5*trace(Q));




N1_p2 = [t, aa-1];
N1_p2 = N1_p2 ./ norm(N1_p2)

N2_p2 = [t, bb-1];
N2_p2 = N2_p2 ./ norm(N2_p2)


U = sqrtm(C);
R_1 = F*inv(U);
[V,D] = eig(U);
U_1 = U;

N_1 = [V(1,1),V(1,2)]';
N_2 = [V(2,1),V(2,2)]';


M = [1 0 ; 0 1];


N_1 = [0.8507,0.5257]';
N_2 = [-0.5257,0.8507]';
M_test = (1/D(2,2)^2)*N_1*N_1' + (1/D(1,1)^2)*N_2*N_2';



a = D(1,1)*rad;
b = D(2,2)*rad;



theta=acosd(V(1,2));

theta = -theta;
theta = 90+theta;


%theta = 31.7175;
R = [cosd(theta) -sind(theta) ; sind(theta) cosd(theta)];

S = [1/D(2,2)^2 0 ; 0 1/D(1,1)^2];


M_elli = R*S*M*R';

alpha = linspace(0,360,50);
%theta = 31.7175;

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
