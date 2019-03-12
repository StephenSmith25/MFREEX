
close all
clear all
clc

MI = [9 3 ; 3 9];


[V,D] = eig(MI);


num_points = 100;

t = linspace(-pi,pi,num_points);
a = sqrt(D(1,1));
b = sqrt(D(2,2));
n1 = V(:,1);
n2 = V(:,2);

figure
plot(0,0,'bo')

hold on
plot([0,5],[0,0],'k-');
hold on
plot([0,0],[0,5],'b-');



test_point = [3,0];

hold on
plot(test_point(1),test_point(2),'g.','markersize',10)

hI1 = (a*n1);
hI2 = b*n2;

hold on
plot([0,hI1(1)],[0,hI1(2)],'k--','linewidth',2);

hold on
plot([0,hI2(1)],[0,hI2(2)],'b--','linewidth',2);

h1 = a;
h2 = b;

for i = 1:num_points
   
   X = [a*cos(t(i)),b*sin(t(i))]';
   
   x_trans = V'*X;
   

   
   hold on
   plot(x_trans(1),x_trans(2),'r.');
   
   
       
       
end

axis square

num_tests = 50;
test_point_x = linspace(-3,3,num_tests);
test_point_y = linspace(0,0,num_tests);

test_points = [test_point_x',test_point_y'];

invMI = inv(MI)

for i = 1:length(test_points)
   
    xS = test_points(i,:);
 
    dI = xS * invMI * xS';
    
    
    r = sqrt(dI);
    drdx(i) =  ( 1 / (2*r)) * (2*xS(1)*invMI(1,1) + invMI(1,2)*xS(2))
   % drdx(i) = xS(1)/(r*9);
    
    
    [w(i),dw(i)] = cubic_spline(r);
    dwdx(i) = dw(i)*drdx(i);
    %dwdx(i) = dw(i) * sign(xS(1));
    
    
    
end
figure
plot(test_points(:,1),w)
hold on
plot(test_points(:,1),dwdx)


function [w,dw] = cubic_spline(r)


    if ( r <= 0.5)
   
        
        w = 2/3 - 4*r*r + 4*r*r*r;
        dw = -8*r + 12*r*r;
        
    elseif ( r > 0.5 ) && ( r <= 1) 
        
        w = 4/3 - 4*r + 4*r*r - (4/3)* r*r*r;
        dw = -4 + 8*r - (4)*r*r;
    else
        w = 0;
        dw = 0;
    end


end



