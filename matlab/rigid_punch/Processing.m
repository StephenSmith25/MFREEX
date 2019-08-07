clear all

close all 

cells = csvread('./../../build/bin/rigid_punch/search_cells.csv');





path = './../../build/bin/rigid_punch/Displacement';
path_base = './../../build/bin/rigid_punch/';
addpath(path)
displacementdir = path ;
d = dir(displacementdir);
d1 = dir([displacementdir,'*.txt']);
numFiles = size(d,1) -3 ;

plotFiles = ceil(linspace(1,numFiles,100));

boundaryNodes = csvread('/home/stephen/Documents/Meshless/build/bin/rigid_punch/boundary.txt');
boundaryNodes = [boundaryNodes ; boundaryNodes(1)]
figure



%-------------------------------------------------------------------------%
%                          Plot 1 
%-------------------------------------------------------------------------%
subplot(1,2,1)

filename = strcat(path_base,'/MaterialPoints/materialpoints_',num2str(plotFiles(1)),'.txt');
material_points = csvread(filename,1);


ix = find(material_points(:,2) > 9.5);








filename = strcat(path,'/displacement_',num2str(plotFiles(1)),'.txt');
disp = csvread(filename);
plot(disp(:,1),disp(:,2),'k.')           % line plot
hold on
plot(disp(boundaryNodes,1),disp(boundaryNodes,2),'b-')

intitial_boundary = disp(boundaryNodes,:);
intitial_nodes = disp;

hold on
axis equal
hold on
%plot(material_points(:,1),material_points(:,2),'b*');
xlim([0,40])
ylim([0,15])

filename = strcat(path_base,'/MaterialPoints/Domains/domains_',num2str(plotFiles(1)),'.txt');
domains = csvread(filename);
general_ellipse_drawer = @(t) draw_general_ellipse(domains(t,:),material_points(t,1),material_points(t,2));

ix = find(material_points(:,2) > 8);
ia = find(material_points(ix,1) > 8);
ix = ix(ia);
ia = find(material_points(ix,1) <12);
ix = ix(ia);


iy = find(disp(:,2) > 8);
ia = find(disp(iy,1) > 7);
iy = iy(ia);
ia = find(disp(iy,1) <12);
iy = iy(ia);

disp_1 = disp(iy,:);
mat_1 = material_points(ix,:); 




ellipses_1 = [];
% for i = 1:3:length(ix)
%     j  = ix(i);
%    ellipses_1 = [ellipses_1 ; general_ellipse_drawer(j)]; 
% end

%-------------------------------------------------------------------------%
%                          Plot 2 
%-------------------------------------------------------------------------%





%material_points = csvread('./../../build/bin/cylinder/Domains/domainmaterialpoints.csv');

subplot(1,2,2)
filename = strcat(path_base,'/MaterialPoints/materialpoints_',num2str(plotFiles(98)),'.txt');
material_points = csvread(filename,1);





filename = strcat(path,'/displacement_',num2str(plotFiles(97)),'.txt');
disp = csvread(filename);
plot(disp(:,1),disp(:,2),'k.')           % line plot
axis equal
hold on
final_boundary = disp(boundaryNodes,:);
final_nodes = disp;

%plot(material_points(:,1),material_points(:,2),'r.');
xlim([0,40])
ylim([0,15])
hold on
plot(disp(boundaryNodes,1),disp(boundaryNodes,2),'b-')
filename = strcat(path_base,'/MaterialPoints/Domains/domains_',num2str(plotFiles(98)),'.txt');
domains = csvread(filename);
general_ellipse_drawer = @(t) draw_general_ellipse_alt(domains(t,1:4),domains(t,5),material_points(t,1),material_points(t,2));
ellipses = [];
% for i = 1:3:length(ix)
%     j  = ix(i);
%    ellipses = [ellipses ; general_ellipse_drawer(j)];  
% end

disp_2 = disp(iy,:);
mat_2 = material_points(ix,:); 


for k = 1:length(cells)
hold on
    i = k;
rectangle('Position',[cells(i,1), cells(i,3), cells(i,2)-cells(i,1), cells(i,4)-cells(i,3)])
end




saveas(gcf,'Displacement_cylinder','epsc')




figure

filename = strcat(path_base,'/MaterialPoints/Domains/domains_',num2str(plotFiles(1)),'.txt');
domains = csvread(filename);
filename = strcat(path_base,'/MaterialPoints/materialpoints_',num2str(plotFiles(1)),'.txt');
material_points = csvread(filename,1);



general_ellipse_drawer = @(t) draw_general_ellipse(domains(t,:),material_points(t,1),material_points(t,2));


iaa = find(disp_1(:,2) > 9.90);

plot_point = 39;

subplot(1,2,1)
hold on
plot(disp_1(:,1),disp_1(:,2),'bo');
hold on
%plot(mat_1(:,1),mat_1(:,2),'r*');
boundary_1 = disp_1(iaa,:);
hold on
plot(boundary_1(:,1),boundary_1(:,2),'b');
hold on
plot(mat_1(plot_point,1),mat_1(plot_point,2),'gx');

r =domains(ix(plot_point),:);

r = 1/sqrt(r(1));
x = mat_1(plot_point,1:2);
a = disp_1 - x;

ica = find(sqrt(sum(a.^2, 2)) <= r);

hold on


plot(disp_1(ica,1),disp_1(ica,2),'g.','markersize',10);


 ellipses_1 = general_ellipse_drawer(ix(plot_point));  

 
 
initial_connectivity = disp_1(ica,:);
initial_zoom_nodes = disp_1;
initial_zoom_ellipse = ellipses_1;
initial_boundary = boundary_1;
initial_mat_point = mat_1(plot_point,1:2);
 
 

axis equal
subplot(1,2,2)

filename = strcat(path_base,'/MaterialPoints/Domains/domains_',num2str(plotFiles(98)),'.txt');
domains = csvread(filename);

filename = strcat(path_base,'/MaterialPoints/materialpoints_',num2str(plotFiles(98)),'.txt');
material_points = csvread(filename,1);



general_ellipse_drawer = @(t) draw_general_ellipse(domains(t,:),material_points(t,1),material_points(t,2));


hold on
plot(disp_2(:,1),disp_2(:,2),'bo');
hold on
%plot(mat_2(:,1),mat_2(:,2),'r*');


hold on
plot(mat_2(plot_point,1),mat_2(plot_point,2),'gx');
axis equal
ellipses_2 = general_ellipse_drawer(ix(plot_point)) ;
boundary_2 = disp_2(iaa,:);

hold on
plot(boundary_2(:,1),boundary_2(:,2),'b');
hold on
plot(disp_2(ica,1),disp_2(ica,2),'g.','markersize',10);


final_connectivity = disp_2(ica,:);
final_zoom_nodes = disp_2;
final_zoom_ellipse = ellipses_2;
final_zoom_boundary = boundary_2;
final_mat_point = mat_2(plot_point,1:2);
 




saveas(gcf,'Displacement_cylinder','epsc')


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
