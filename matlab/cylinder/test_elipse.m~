
clear all 
close all 



MI = [ 9 0;0 9];

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
test_point = [linspace(-1.118,1.118,num_tests)',linspace(1.118,-1.118,num_tests)'];
x = linspace(-3,3,num_tests);
y = linspace(-3,3,num_tests);

invMI = inv(MI);


[X,Y] = meshgrid(x,y);
for i = 1:num_tests
    
    for j = 1:num_tests
  
    x_proj = [X(i,j),Y(i,j)]*hI1/(h1);
    y_proj = [X(i,j),Y(i,j)]*hI2/(h2);
   
    r_x(i,j) = (x_proj/h1);
    r_y(i,j) = (y_proj/h2);  
    
    dr_dx(i,j) = sign(-x_proj);
    dr_dy(i,j) = sign(y_proj);
    
    
    
    [w_x(i,j),dw_X(i,j)] = cubic_spline(abs(r_x(i,j)));
    [w_y(i,j),dw_Y(i,j)] = cubic_spline(abs(r_y(i,j))); 
    
    
    dwdX(i,j) = dr_dx(i,j) * dw_X(i,j);
    dwdY(i,j) = dr_dy(i,j) * dw_Y(i,j);
    
    dw_dx(i,j) = n1(1) * dwdX(i,j)*w_y(i,j) + n1(2) * w_x(i,j) * dwdY(i,j);
    
    dw_dy(i,j) = n2(1) * dwdX(i,j)*w_y(i,j) + n2(2) * w_x(i,j) * dwdY(i,j);

      
    end

    
    
end
figure
surf(X,Y,w_x.*w_y);

figure
surf(X,Y,dw_dx);

% figure
% subplot(1,2,1)
% plot(r_x,w_x.*w_y);
% hold on
% plot(r_x,dw_dx);
% 
% subplot(1,2,2)
% plot(r_y,w_y);
% hold on
% plot(r_y,dw_dy);



% second method
invMI = inv(MI);



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





