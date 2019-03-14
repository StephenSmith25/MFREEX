close all
clear all 


a0 = -2.2264e7;
a1 = 2.3704e7;
a2 = -9.3769e6;
a3 = 1.6212e6;
a4 = -9.73804e4;
a5 = -1.8801e3;
a6 = 559.3131;
a7 = 0.2565;


p = @(x) a0*x.^7 + a1 * x.^6 + a2*x.^5 + a3*x.^4 + a4*x.^3 + a5*x.^2 + a6*x + a7;


t = linspace(0,0.31,50);

y = p(t);

y= [y,y(end)];
t =[t,1];

plot(t,y)

data = [t',y'];