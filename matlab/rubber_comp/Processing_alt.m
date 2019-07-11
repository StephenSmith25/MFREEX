clear all

close all 
path = './../../build/bin/rubber_comp/Displacement';
path_base = './../../build/bin/rubber_comp/';
addpath(path)
displacementdir = path ;
d = dir(displacementdir);
d1 = dir([displacementdir,'*.txt']);
numFiles = size(d,1) -3 ;

plotFiles = ceil(linspace(1,numFiles,100));



figure
subplot(1,3,1)
boundary_nodes = csvread('./../../build/bin/rubber_comp/boundary.txt');
boundary_nodes = [boundary_nodes;boundary_nodes(1)];

nodes = csvread('./../../build/bin/rubber_comp/nodes.csv');

top_boundary = find(nodes(:,2) == 0.1);

bot_boundary = find(nodes(:,2) == 0);



%-------------------------------------------------------------------------%
%                          Plot 1 
%-------------------------------------------------------------------------%


filename = strcat(path_base,'/MaterialPoints/materialpoints_',num2str(plotFiles(1)),'.txt');
%material_points = csvread(filename,1);


filename = strcat(path,'/displacement_',num2str(plotFiles(2)),'.txt');
disp = nodes;
plot(disp(:,1),disp(:,2),'k.')           % line plot
axis equal
hold on
%plot(material_points(:,1),material_points(:,2),'b*');

% 
% filename = strcat(path_base,'/MaterialPoints/Domains/domains_',num2str(plotFiles(1)),'.txt');
% domains = csvread(filename);
% general_ellipse_drawer = @(t) draw_general_ellipse(domains(t,:),material_points(t,1),material_points(t,2));
% 
% general_ellipse_drawer(200);
% 
% ellipses = [];
% 
% for i = 1:200:length(domains)
%     
%    ellipses = [ellipses ; general_ellipse_drawer(i)]; 
% end


start_disp = [disp];
start_top_boundary = disp(top_boundary,:);
start_bot_boundary = disp(bot_boundary,:);

plot(start_top_boundary(:,1),start_top_boundary(:,2),'b','linewidth',3);
plot(start_bot_boundary(:,1),start_bot_boundary(:,2),'b','linewidth',3);hold on



%-------------------------------------------------------------------------%
%                          Plot 2 
%-------------------------------------------------------------------------%




subplot(1,3,2)


filename = strcat(path,'/displacement_',num2str(plotFiles(38)),'.txt');
disp = csvread(filename);
hold on 
plot(disp(:,1),disp(:,2),'k.')           % line plot
axis equal
hold on
%plot(material_points(:,1),material_points(:,2),'r.');




mid_disp = [disp];
mid_top_boundary = disp(top_boundary,:);
mid_bot_boundary = disp(bot_boundary,:);
hold on
plot(mid_top_boundary(:,1),mid_top_boundary(:,2),'b','linewidth',3);
hold on
plot(mid_bot_boundary(:,1),mid_bot_boundary(:,2),'b','linewidth',3);
%material_points = csvread('./../../build/bin/cylinder/Domains/domainmaterialpoints.csv');

%material_points = csvread(filename,1);


subplot(1,3,3)

filename = strcat(path,'/displacement_',num2str(plotFiles(100)),'.txt');
disp = csvread(filename);
hold on 
plot(disp(:,1),disp(:,2),'k.')           % line plot

axis equal



%plot(material_points(:,1),material_points(:,2),'r.');




final_disp = [disp];
final_top_boundary = disp(top_boundary,:);
final_bot_boundary = disp(bot_boundary,:);

plot(final_top_boundary(:,1),final_top_boundary(:,2),'b','linewidth',3);
plot(final_bot_boundary(:,1),final_bot_boundary(:,2),'b','linewidth',3);


% %filename = strcat(path_base,'/MaterialPoints/Domains/domains_',num2str(plotFiles(10)),'.txt');
% domains = csvread(filename);
% general_ellipse_drawer = @(t) draw_general_ellipse_alt(domains(t,1:4),domains(t,5),material_points(t,1),material_points(t,2));
% ellipses = [];
% for i = 1:200:length(domains)
%     
%    ellipses = [ellipses ; general_ellipse_drawer(i)]; 
%     
% end


saveas(gcf,'Displacement_rubberr','epsc')

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
