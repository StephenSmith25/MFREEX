clear all
close all

figure

% test general ellipse
MI = [3,2,2,4];
draw_general_ellipse(MI,0,0);
axis equal

function [] = draw_general_ellipse(M,x0,y0)

    t = linspace(0,360,50);
    
    MI = [M(1),M(2); M(3), M(4)];
    
    [V,D] = eig(MI);
    
    
    
    
    
    
    a = sqrt(D(1,1));
    b = sqrt(D(2,2));
    x = x0 + a*cosd(t);
    y = y0 + b*sind(t);
    X = [x;y];
    
    
    X_prime = V*X;
    
    hold on
    plot(X_prime(1,:),X_prime(2,:),'r-');
    
    
    
    


end