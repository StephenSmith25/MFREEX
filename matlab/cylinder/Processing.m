clear all

close all 
path = './../../build/bin/cylinder/Displacement';
path_base = './../../build/bin/cylinder/';
addpath(path)
displacementdir = path ;
d = dir(displacementdir);
d1 = dir([displacementdir,'*.txt']);
numFiles = size(d,1) -3 ;

plotFiles = ceil(linspace(1,numFiles,10));



figure



%-------------------------------------------------------------------------%
%                          Plot 1 
%-------------------------------------------------------------------------%


filename = strcat(path_base,'/MaterialPoints/materialpoints_',num2str(plotFiles(1)),'.txt');
material_points = csvread(filename,1);


filename = strcat(path,'/displacement_',num2str(plotFiles(1)),'.txt');
disp = csvread(filename);
plot(disp(:,1),disp(:,2),'k.')           % line plot
axis equal
hold on
%plot(material_points(:,1),material_points(:,2),'b*');
xlim([0,40])
ylim([0,40])

filename = strcat(path_base,'/MaterialPoints/Domains/domains_',num2str(plotFiles(1)),'.txt');
domains = csvread(filename);
general_ellipse_drawer = @(t) draw_general_ellipse(domains(t,:),material_points(t,1),material_points(t,2));

general_ellipse_drawer(200);

ellipses = [];

for i = 1:200:length(domains)
    
   ellipses = [ellipses ; general_ellipse_drawer(i)]; 
end

%-------------------------------------------------------------------------%
%                          Plot 2 
%-------------------------------------------------------------------------%





%material_points = csvread('./../../build/bin/cylinder/Domains/domainmaterialpoints.csv');


filename = strcat(path_base,'/MaterialPoints/materialpoints_',num2str(plotFiles(10)),'.txt');
material_points = csvread(filename,1);





filename = strcat(path,'/displacement_',num2str(plotFiles(10)),'.txt');
disp = csvread(filename);
plot(disp(:,1),disp(:,2),'k.')           % line plot
axis equal
hold on
%plot(material_points(:,1),material_points(:,2),'r.');
xlim([0,40])
ylim([0,40])

filename = strcat(path_base,'/MaterialPoints/Domains/domains_',num2str(plotFiles(10)),'.txt');
domains = csvread(filename);
general_ellipse_drawer = @(t) draw_general_ellipse_alt(domains(t,1:4),domains(t,5),material_points(t,1),material_points(t,2));
ellipses = [];
for i = 1:200:length(domains)
    
   ellipses = [ellipses ; general_ellipse_drawer(i)]; 
    
end


saveas(gcf,'Displacement_cylinder','epsc')



figure

filename = strcat(path,'/displacement_',num2str(plotFiles(10)),'.txt');
disp = csvread(filename);
plot(disp(:,1),disp(:,2),'k.')           % line plot
axis equal
hold on
%plot(material_points(:,1),material_points(:,2),'r.');
xlim([0,40])
ylim([0,40])





saveas(gcf,'Displacement_cylinder','epsc')
%-------------------------------------------------------------------------%
%                          Load Displacement
%-------------------------------------------------------------------------%
sizeA = [2 inf];
formatSpec = '%f %f';
fileID = fopen('./../../build/bin/cylinder/loadDisp.txt','r');
A = fscanf(fileID,formatSpec,sizeA);

A = A';

figure



plot(A(1:25:end,1),A(1:25:end,2),'kx','markersize',6);
% 
sizeA = [2 inf];
formatSpec = '%f %f';
fileID = fopen('./exactSol.txt','r');

B = fscanf(fileID,formatSpec,sizeA);

B = B';

hold on 
plot(B(:,1),B(:,2)/1000,'r-');
xlabel('dr')
ylabel('Pressure')
legend('Meshfree','Exact','Location','northwest');

saveas(gcf,'Solution_cylinder','epsc')

latex_var_1 = [B(:,1),B(:,2)/1000];
latex_var_2 = A(1:25:end,:);


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



function [X_prime] = draw_general_ellipse(M,x0,y0)

    t = linspace(0,360,50);
    
    MI = [M(1),M(2); M(3), M(4)];
    
    [V,D] = eig(MI);
    
    
    
    
    
    
    a = 1.00/sqrt(D(1,1));
    b = 1.00/sqrt(D(2,2));
    x = a*cosd(t);
    y = b*sind(t);
    X = [x;y];
    
    X_prime = V'*X;
    
   
    
    hold on
    plot(X_prime(1,:)+x0,X_prime(2,:)+y0,'r-');
    X_prime = [X_prime(1,:)'+x0, X_prime(2,:)'+y0];
    
    
    
    


end
