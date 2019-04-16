
clear all
close all 

X = linspace(-1,1,10);
Y = linspace(-1,1,10);


[x,y] = meshgrid(X,Y);



x = x(:);
y = y(:);

for ( i = 1:length(x))
    
    
    r = sqrt(x(i)^2 + y(i)^2);
    
    
    w(i) = cubic_spline(r);
    
end

a = [x,y,w'];
    


function [w,dw] = cubic_spline(r)


    if ( r <= 0.5)
   
        
        w = 2/3 - 4.*r.*r + 4.*r.*r.*r;
        dw = -8*r + 12*r*r;
        
    elseif ( r > 0.5 ) && ( r <= 1) 
        
        w = 4/3 - 4.*r + 4.*r.*r - (4/3).* r.*r.*r;
        dw = -4 + 8.*r - (4).*r.*r;
    else
        w = 0;
        dw = 0;
    end


end
