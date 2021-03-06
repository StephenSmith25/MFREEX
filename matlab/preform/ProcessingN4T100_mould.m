clear all
close all

PLOT_GRAPHS = true;
PLOT_DOMAINS_INFLUENCE = true; 


WITH_MOULD = true;
TMAX = 0.6;



%DOMAIN_TYPE = 'RECTANGULAR';
DOMAIN_TYPE = 'RADIAL'



path = './../../build/bin/preform/Displacement/';
pathSR = './../../build/bin/preform/srRod/';
boundaryNodes = csvread('./../../build/bin/preform/boundary.txt');
Domains= csvread('./../../build/bin/preform/domains.txt');

boundaryNodes = [boundaryNodes;boundaryNodes(1)];

addpath(path)
displacementdir = path ;
d = dir(displacementdir);
d1 = dir([displacementdir,'*.csv']);
numFiles = size(d,1) -3 ;

plotFiles = ceil(linspace(1,numFiles,10));

plot_point =141;
filename = strcat(path,'displacement_',num2str(plotFiles(1)),'.csv');
disp = csvread(filename,1);

% draw mould
Radius_mould = 24;
Height_mould = 30;
thickness_mould = 5;
mould_length = 180;
height = 74.2650;
Radius_out = 11.25;
top_point = 48.920;



mould_nodes = [];
mould_nodes = [mould_nodes ; [Radius_out,height+5]   ];
mould_nodes = [mould_nodes ; [Radius_out,68.92]   ];
mould_nodes = [mould_nodes ; [Radius_out+Radius_mould,top_point]   ];
mould_nodes = [mould_nodes ; [Radius_out+Radius_mould,height-mould_length]   ];
mould_nodes = [mould_nodes ; [0,height-mould_length]   ];


mould_nodes = [mould_nodes ; [0,height-mould_length-thickness_mould]   ];
mould_nodes = [mould_nodes ; [Radius_out+Radius_mould+thickness_mould,height-mould_length-thickness_mould]   ];
mould_nodes = [mould_nodes ; [Radius_out+Radius_mould+thickness_mould,top_point]   ];
mould_nodes = [mould_nodes ; [Radius_out+thickness_mould,68.92]   ];
mould_nodes = [mould_nodes ; [Radius_out+thickness_mould,height+5]   ];
mould_nodes = [mould_nodes ; [Radius_out,height+5]   ];

if ( WITH_MOULD == false)
    mould_nodes = zeros(1,2);
end


% 
% for i = 1:length(plotFiles)
%    
%    fig = figure;
% 
%    filename = strcat(path,'displacement_',num2str(plotFiles(i)),'.csv');
%     disp = csvread(filename,1);
%     subplot(1,3,2)       % add first plot in 2 x 2 grid
%     hold on 
%     fill(mould_nodes(:,1),mould_nodes(:,2),'w');
%     hold on
%     fill(-mould_nodes(:,1),mould_nodes(:,2),'w');
%     hold on
% 
% 
% 
%     plot(disp(:,1),disp(:,2),'k.','markersize',3)           % line plot
%     hold on
%     plot(-disp(:,1),disp(:,2),'k.','markersize',3) % line plot
%     hold on
%     plot(disp(plot_point,1),disp(plot_point,2),'r*')
%     axis equal
%     hold on
% 
%     hold on
%     filename = strcat(pathSR,'srRod_',num2str(plotFiles(i)),'.csv');
%     disp = csvread(filename,1);
%     hold on
%     plot(disp(:,1),disp(:,2),'r-','linewidth',1)           % line plot
%     hold on
%     plot(-disp(:,1),disp(:,2),'r-','linewidth',1)           % line plot
%     xlim([-50,50])
%     ylim([-170,78])
% 
%     filename = strcat(path,'displacement_',num2str(plotFiles(10)),'.csv');
%     disp = csvread(filename,1);
% 
%     filename = strcat('./Displacements/Displacement_',num2str(i),'.png');
%     
%     drawnow   
%     
%     saveas(fig,filename);
%    
%    
%    
%    close all
%     
%     
% end


figure




subplot(1,3,1)       % add first plot in 2 x 2 grid

hold on 
fill(mould_nodes(:,1),mould_nodes(:,2),'w');
hold on
fill(-mould_nodes(:,1),mould_nodes(:,2),'w');
hold on
plot(disp(:,1),disp(:,2),'k.','markersize',3)           % line plot
hold on
plot(-disp(:,1),disp(:,2),'k.','markersize',3)           % line plot
ymax = max(disp(:,2));
plot(disp(plot_point,1),disp(plot_point,2),'r*')



% if (strcmp(DOMAIN_TYPE,'RADIAL') == 1)
%    
%    [x,y] = circle(disp(plot_point,1),disp(plot_point,2),Domains(plot_point));
%      
% 
%    hold on
%    plot(x,y,'r.');
%    
%    hold on
%    plot(disp(plot_point,1),disp(plot_point,2),'bo')
%    
%    
% elseif( strcmp(DOMAIN_TYPE,'RECTANGULAR') == 1)
%         
%      hold on
%      rectangle('Position',[disp(plot_point,1)-Domains(plot_point,1),disp(plot_point,2)-Domains(plot_point,2),Domains(plot_point,1)*2,....
%          Domains(plot_point,2)*2]);
%        hold on
%    plot(disp(plot_point,1),disp(plot_point,2),'bo')
%    
% elseif ( strcmp(DOMAIN_TYPE, 'ELLIPTICAL') == 1)
%      invMI = [Domains(plot_point,1:2) ; Domains(plot_point,3:4)];
% 
%     MI = inv(invMI);
% 
% 
% [V,D] = eig(MI);
% 
% 
% [x_rotated,y_rotated] = ellipse(0,0,D(1,1),D(2,2));
% xyRotated = [x_rotated',y_rotated']*V';
% xy = xyRotated;
% 
% x = xy(:,1) + disp(plot_point,1);
% y = xy(:,2) + disp(plot_point,2);
% 
% 
% hold on
% plot(x,y,'r-');
% axis equal
% else
%     
% end  


axis equal 
hold on
%plot(disp(boundaryNodes,1),disp(boundaryNodes,2),'b-')
hold on
%plot(-disp(boundaryNodes,1),disp(boundaryNodes,2),'b-')


boundary_nodes_xy = disp(boundaryNodes,1:2);
Rout = max(disp(:,1));
height = max(disp(:,2));




hold on
filename = strcat(pathSR,'srRod_',num2str(plotFiles(1)),'.csv');
disp = csvread(filename,1);
hold on
plot(disp(:,1),disp(:,2),'r-','linewidth',1)           % line plot
hold on
plot(-disp(:,1),disp(:,2),'r-','linewidth',1)           % line plot
xlim([-50,50])
ylim([-170,ymax])

       % line plot




filename = strcat(path,'displacement_',num2str(plotFiles(5)),'.csv');
disp = csvread(filename,1);
subplot(1,3,2)       % add first plot in 2 x 2 grid
hold on 
fill(mould_nodes(:,1),mould_nodes(:,2),'w');
hold on
fill(-mould_nodes(:,1),mould_nodes(:,2),'w');
hold on



plot(disp(:,1),disp(:,2),'k.','markersize',3)           % line plot
hold on
plot(-disp(:,1),disp(:,2),'k.','markersize',3) % line plot
hold on
%plot(disp(plot_point,1),disp(plot_point,2),'r*')
axis equal
hold on
%plot(disp(boundaryNodes,1),disp(boundaryNodes,2),'b-')
hold on
%plot(-disp(boundaryNodes,1),disp(boundaryNodes,2),'b-')
%plot(disp(boundaryNodes,1),disp(boundaryNodes,2),'b-')
hold on
filename = strcat(pathSR,'srRod_',num2str(plotFiles(5)),'.csv');
disp = csvread(filename,1);
hold on
plot(disp(:,1),disp(:,2),'r-','linewidth',1)           % line plot
hold on
plot(-disp(:,1),disp(:,2),'r-','linewidth',1)           % line plot
xlim([-50,50])
ylim([-170,ymax])

filename = strcat(path,'displacement_',num2str(plotFiles(10)),'.csv');
disp = csvread(filename,1);


subplot(1,3,3)       % add first plot in 2 x 2 grid

hold on 
fill(mould_nodes(:,1),mould_nodes(:,2),'w');
hold on
fill(-mould_nodes(:,1),mould_nodes(:,2),'w');
hold on

plot(disp(:,1),disp(:,2),'k.','markersize',3)           % line plot
hold on
plot(-disp(:,1),disp(:,2),'k.','markersize',3)           % line plot
hold on
%plot(disp(boundaryNodes,1),disp(boundaryNodes,2),'b-')
hold on
%plot(-disp(boundaryNodes,1),disp(boundaryNodes,2),'b-')
hold on
%plot(disp(plot_point,1),disp(plot_point,2),'r*')
axis equal
hold on
filename = strcat(pathSR,'srRod_',num2str(plotFiles(10)),'.csv');
disp = csvread(filename,1);
hold on
plot(disp(:,1),disp(:,2),'r-','linewidth',1)           % line plot
hold on
plot(-disp(:,1),disp(:,2),'r-','linewidth',1)           % line plot
xlim([-50,50])
ylim([-170,ymax])



% 
% figure
% filename = strcat(path,'displacement_',num2str(plotFiles(10)),'.csv');
% 
% 
% disp = csvread(filename,1);
% 
% disp = disp(boundary_nodes,1:2);
% boundary_nodes = boundaryNodes(1:end-1);
% boundary_nodes = reshape(boundary_nodes,2,(length(boundary_nodes))/2)';
% 
% boundary_nodes = boundary_nodes - iy;
% 
% 
% % 
% % 
% %  
% %  mould_nodes_index = zeros(5,2);
% % for i = 1:5
% %     mould_nodes_index(i,:) = [i,i+1];
% %     
% %     if ( i == length(mould_nodes))
% %         mould_nodes_index(i,:) = [i,1];
% %     end
% % end
% % 
% %  num_revolves = 32;
% % [tri, xyz] = bfRevolve(mould_nodes_index, mould_nodes, num_revolves);
% % 
% % for i = 1:length(mould_nodes_index) 
% %     for j = 1:round(num_revolves/2)
% %         tri_mod((i-1)*num_revolves/2 + j,:) = tri((i-1)*num_revolves + j,:);  
% %     end
% %     
% %     end
% % s = trisurf(tri_mod, xyz(:,1),xyz(:,2),xyz(:,3),'FaceColor','y');
% % 
% % 
% % 
% [tri, xyz] = bfRevolve(boundary_nodes, disp(:,1:2), num_revolves);
% 
% % display the surface
% hold on
% trisurf(tri, xyz(:,1),xyz(:,2),xyz(:,3),'FaceColor','c');
% 
% axis equal
% axis([-40 40 -40 40 -100 77])
% 

%axis equal



print -dpng2 Displacement.png

m = dlmread('./../../build/bin/preform/pressureTime.txt',' ');
%m(ceil(length(m)/4):end,:) = smoothdata(m(ceil(length(m)/4):end,:))
figure
plot(m(:,1),m(:,2),'k-','linewidth',3);
xlabel('time');
ylabel('Pressure');
xlim([0,TMAX])
ylim([0,5])

a1 = [m(1:2:end,1),m(1:2:end,2)];

b = csvread("Experimental/N4T100_exp_pressure.csv");
hold on
plot(b(:,1),b(:,2)/1000,'b--')


c = csvread("Experimental/N4T100_fe_pressure.csv");

hold on
plot(c(:,1),c(:,2)/1000,'r--')

legend('Mesh-free','Experiemntal','FE (Abaqus)')
ylabel('Pressure (MPa)')
xlabel('Time (s)')


%set axis
set(gca, 'FontName', 'cmr12')
% set x tics and y tics
set(gca,...
'Units','normalized',...
'FontWeight','normal',...
'FontSize',14,... % size ofiguref numbers on axis
'FontName','cmr14') % font name
set(gca,'TickLabelInterpreter', 'latex');

% Y label
ylabel({'Pressure (MPa) '},...
'interpreter','latex',...
'FontSize',14,... % font size
'FontName','cmr14')
% X label
xlabel('Time (s)',...
'interpreter','latex',...
'FontWeight','normal',...
'FontSize',14,... % font size
'FontName','cmr14')
% legend

legend({'Meshfree','Experimental','FE (Abaqus )'},... % { 'legend1', 'legend2',...}
'interpreter','latex',...
'FontSize',12,...
'FontName','cmr14',...
'Location','SouthEast')


set(gcf, 'Color', 'w');

print -dpng2 pressure.png



set(gcf, 'Color', 'w');

if ( PLOT_GRAPHS)

path = './../../build/bin/preform/History/Strain';

addpath(path)
displacementdir = path ;
d = dir(displacementdir);
d1 = dir([displacementdir,'*.txt']);
numFiles = size(d,1) -3 ;

plotFiles = ceil(linspace(1,numFiles,100));
for i = 1:length(plotFiles)
    
    

    filename = strcat('./../../build/bin/preform/History/Strain/strain_',num2str(plotFiles(i)),'.txt');

    strain = csvread(filename);
    
 
    F =  [strain(1,1) strain(1,2) 0 ;
        strain(2,1) strain(2,2) 0;
        0 0 strain(3,3)]  ;
    
    
    [R U V] = poldecomp(F);
    
    
    true_strain = logm(U);

    
    hoop_strain(i) =true_strain(3,3);
    axial_strain(i) = true_strain(2,2);
    radial_strain(i) = true_strain(1,1);
    shear_strain(i) = true_strain(1,2);
    


    


    
    
    jacobian(i) = det(F);
    crit_lambda(i) = strain(2,3);
    gamma(i) = strain(3,1);
    max_lambda_n(i) = strain(1,3);
    time(i) = strain(3,2);
  
    
end


n = smoothdata(axial_strain(1:end));
axial_strain(40:end) = n(40:end);

a2 = [time(1:1:end)',axial_strain(1:1:end)'];
a3 = [time(1:1:end)',hoop_strain(1:1:end)'];

strain = figure('Position', get(0, 'Screensize'));
subplot(1,2,1)
plot(time,(axial_strain),'k','linewidth',2);
% hold on
% plot(time,log(1+shear_strain),'k');
xlim([0,TMAX])

c = csvread("Experimental/N4100_exp_strain_axial.csv");

hold on
plot(c(:,1),c(:,2),'b--')

c = csvread("Experimental/N4100_fe_strain_axial.csv");

hold on
plot(c(:,1),c(:,2),'r--')


hold on
plot(time,jacobian,'g--')

hold on
plot(time,(shear_strain),'g','linewidth',2);


%set axis
set(gca, 'FontName', 'cmr12')
% set x tics and y tics
set(gca,...
'Units','normalized',...
'FontWeight','normal',...
'FontSize',14,... % size ofiguref numbers on axis
'FontName','cmr14') % font name
set(gca,'TickLabelInterpreter', 'latex');

% Y label
ylabel({'True Strain '},...
'interpreter','latex',...
'FontSize',14,... % font size
'FontName','cmr14')
% X label
xlabel('Time',...
'interpreter','latex',...
'FontWeight','normal',...
'FontSize',14,... % font size
'FontName','cmr14')
% legend

legend({'Axial (Meshfree)', 'Axial (exp)','Axial (FE)', 'Jacobian', 'Shear' },... % { 'legend1', 'legend2',...}
'interpreter','latex',...
'FontSize',12,...
'FontName','cmr14',...
'Location','SouthEast')

grid on 




subplot(1,2,2)
plot(time,(hoop_strain),'k','linewidth',2);


c = csvread("Experimental/N4100_exp_strain_hoop.csv");

hold on
plot(c(:,1),c(:,2),'b--')


c = csvread("Experimental/N4100_fe_strain_hoop.csv");

hold on
plot(c(:,1),c(:,2),'r--')

xlim([0,TMAX])


%set axis
set(gca, 'FontName', 'cmr12')
% set x tics and y tics
set(gca,...
'Units','normalized',...
'FontWeight','normal',...
'FontSize',14,... % size ofiguref numbers on axis
'FontName','cmr14') % font name
set(gca,'TickLabelInterpreter', 'latex');

% Y label
ylabel({'True Strain '},...
'interpreter','latex',...
'FontSize',14,... % font size
'FontName','cmr14')
% X label
xlabel('Time',...
'FontUnits','points',...
'interpreter','latex',...
'FontWeight','normal',...
'FontSize',14,... % font size
'FontName','cmr14')
% legend

legend({'Hoop (Meshfree)','Hoop (exp)','Hoop (FE)' },... % { 'legend1', 'legend2',...}
'interpreter','latex',...
'FontSize',12,...
'FontName','cmr14',...
'Location','SouthEast')

grid on 



print -dpng2 strain.png




% legend
path = './../../build/bin/preform/History/Stress';




addpath(path)
displacementdir = path ;
d = dir(displacementdir);
d1 = dir([displacementdir,'*.txt']);
numFiles = size(d,1)/2 -3 ;


for i = 1:length(plotFiles)
    
    

    filename = strcat('Bond_Stress_',num2str(plotFiles(i)),'.txt');

    stress = csvread(filename);
    hoop_stress_bond(i) = stress(3,3);
    axial_stress_bond(i) = stress(2,2);
    radial_stress_bond(i) = stress(1,1);
    shear_stress_bond(i) = stress(1,2);

    filename = strcat('Conformational_Stress_',num2str(plotFiles(i)),'.txt');

    stress = csvread(filename);
    hoop_stress_conf(i) = stress(3,3);
    axial_stress_conf(i) = stress(2,2);
    radial_stress_conf(i) = stress(1,1);
     shear_stress_conf(i) = stress(1,2);

  
end



figure
subplot(1,2,1)
hold on 
plot(time,hoop_stress_conf,'r');
hold on
plot(time,hoop_stress_bond,'b');
hold on
plot(time,shear_stress_bond,'m')
hold on
plot(time,shear_stress_conf,'g')
%set axis
set(gca, 'FontName', 'cmr12')
% set x tics and y tics
set(gca,...
'Units','normalized',...
'FontWeight','normal',...
'FontSize',14,... % size ofiguref numbers on axis
'FontName','cmr14') % font name
set(gca,'TickLabelInterpreter', 'latex');

% Y label
ylabel({'Cauchy Stress '},...
'interpreter','latex',...
'FontSize',14,... % font size
'FontName','cmr14')
% X label
xlabel('Time',...
'interpreter','latex',...
'FontWeight','normal',...
'FontSize',14,... % font size
'FontName','cmr14')
% legend

legend({ 'Hoop (Conf)','Hoop (Bond)','Shear (Bond)', 'Shear (Conf)' },... % { 'legend1', 'legend2',...}
'interpreter','latex',...
'FontSize',12,...
'FontName','cmr14',...
'Location','SouthEast')

grid on 

subplot(1,2,2)
plot(time,axial_stress_conf,'r');
hold on
plot(time,axial_stress_bond,'b');

%set axis
set(gca, 'FontName', 'cmr12')
% set x tics and y tics
set(gca,...
'Units','normalized',...
'FontWeight','normal',...
'FontSize',14,... % size ofiguref numbers on axis
'FontName','cmr14') % font name
set(gca,'TickLabelInterpreter', 'latex');

% Y label
ylabel({'Cauchy Stress '},...
'interpreter','latex',...
'FontSize',14,... % font size
'FontName','cmr14')
% X label
xlabel('Time',...
'interpreter','latex',...
'FontWeight','normal',...
'FontSize',14,... % font size
'FontName','cmr14')
% legend

legend({ 'Axial (Conf)','Axial (Bond)' },... % { 'legend1', 'legend2',...}
'interpreter','latex',...
'FontSize',12,...
'FontName','cmr14',...
'Location','SouthEast')

grid on 



print -dpng2 stress.png



end
 

function [x,y] = circle(x0,y0,r)
    
    theta = linspace(0,2*pi,40);
    
    x = x0 + r*cos(theta);
    y = y0 + r*sin(theta);

end


function [x,y] = ellipse(x0,y0,a,b)
    
    t=-pi:0.01:pi;
    x=x0+a*cos(t);
    y=y0+b*sin(t);
end

