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



% should last about 0.4 seconds
% v at t=0 = 559.3131;
% v at 0.1 approx 180
% v at 0.25 = 600
% v at end = 0.2


data = csvread('N2T100_SR_PROFILE.csv')
x = data(:,1);
y = data(:,2);
p = polyfit(x,y,7);

a0 = p(1);
a1 = p(2);
a2 = p(3);
a3 = p(4);
a4 = p(5);
a5 = p(6);
a6 = p(7);
a7 = 0;


p = @(x) a0*x.^7 + a1 * x.^6 + a2*x.^5 + a3*x.^4 + a4*x.^3 + a5*x.^2 + a6*x + a7;
p1 = @(x) 7*a0*x.^6 + 6*a1*x.^5 + 5*a2 *x.^4 + 4*a3*x.^3 + a4*3*x.^2 + 2*a5*x + a6;
t = linspace(0,0.5,50);

figure
plot(x,y,'ro');
hold on
plot(t,p(t),'b');


figure
plot(t,p1(t));