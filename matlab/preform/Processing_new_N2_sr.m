
clear all
close all

PLOT_GRAPHS = true;
PLOT_DOMAINS_INFLUENCE = true; 


WITH_MOULD = false;
TMAX = 0.7;




path = './../../build/bin/preform/Displacement/';
pathSR = './../../build/bin/preform/srRod/';
boundaryNodes = csvread('./../../build/bin/preform/boundary.txt');
path_base = './../../build/bin/preform';

boundaryNodes = [boundaryNodes;boundaryNodes(1)];

addpath(path)
displacementdir = path ;
d = dir(displacementdir);
d1 = dir([displacementdir,'*.csv']);
numFiles = size(d,1) -3 ;

plotFiles = ceil(linspace(1,numFiles,10));

plot_point =129; %%2241
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




figure



% -------------------------------------------------------------------------%
%                           PLOT 1
% -------------------------------------------------------------------------%



subplot(1,3,1)       % add first plot in 2 x 2 grid

hold on 
fill(mould_nodes(:,1),mould_nodes(:,2),'w');
hold on
fill(-mould_nodes(:,1),mould_nodes(:,2),'w');
hold on
plot(disp(:,1),disp(:,2),'k.','markersize',3)           % line plot
hold on
plot(-disp(:,1),disp(:,2),'k.','markersize',3)           % line plot
ymax = 0;
hold on
plot(-disp(plot_point,1),disp(plot_point,2),'r*')
hold on
plot(disp(plot_point,1),disp(plot_point,2),'r*')




c = [[disp(:,1);-disp(:,1)],[disp(:,2);disp(:,2)]]








save('disp_tstart.dat', 'c', '-ascii', '-double', '-tabs')



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
plot(disp(:,1),disp(:,2),'r');
hold on
plot(-disp(:,1),disp(:,2),'r');


c = [[disp(1:end-1,1);-disp(end:-1:1,1)],[disp(1:end-1,2);disp(end:-1:1,2)]]


save('disp_rod_tstart.dat', 'c', '-ascii', '-double', '-tabs')





hold on
xlim([-50,50])
ylim([-300,ymax])



% -------------------------------------------------------------------------%
%                           PLOT 2
% -------------------------------------------------------------------------%




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
hold on
plot(-disp(plot_point,1),disp(plot_point,2),'r*')
hold on
plot(disp(plot_point,1),disp(plot_point,2),'r*')


c = [[disp(:,1);-disp(:,1)],[disp(:,2);disp(:,2)]]
save('disp_tmid.dat', 'c', '-ascii', '-double', '-tabs')




%plot(disp(plot_point,1),disp(plot_point,2),'r*')
axis equal
hold on
%plot(disp(boundaryNodes,1),disp(boundaryNodes,2),'b-')
hold on
%plot(-disp(boundaryNodes,1),disp(boundaryNodes,2),'b-')
%plot(disp(boundaryNodes,1),disp(boundaryNodes,2),'b-')
hold on

hold on
filename = strcat(pathSR,'srRod_',num2str(plotFiles(5)),'.csv');
disp = csvread(filename,1);
hold on
plot(disp(:,1),disp(:,2),'r');
hold on
plot(-disp(:,1),disp(:,2),'r');


c = [[disp(1:end-1,1);-disp(end:-1:1,1)],[disp(1:end-1,2);disp(end:-1:1,2)]]


save('disp_rod_tmid.dat', 'c', '-ascii', '-double', '-tabs')





xlim([-50,50])
ylim([-300,ymax])



print -dpng2 Displacement.png

% -------------------------------------------------------------------------%
%                           PLOT 3
% -------------------------------------------------------------------------%
filename = strcat(path,'displacement_',num2str(plotFiles(10)),'.csv');
disp = csvread(filename,1);




subplot(1,3,3)       % add first plot in 2 x 2 grid

hold on 
fill(mould_nodes(:,1),mould_nodes(:,2),'w');
hold on
fill(-mould_nodes(:,1),mould_nodes(:,2),'w');
hold on

plot(disp(:,1),disp(:,2),'k.','markersize',5)           % line plot
hold on
plot(-disp(:,1),disp(:,2),'k.','markersize',5)           % line plot
hold on
%plot(disp(boundaryNodes,1),disp(boundaryNodes,2),'b-')
hold on
%plot(-disp(boundaryNodes,1),disp(boundaryNodes,2),'b-')
hold on
%plot(disp(plot_point,1),disp(plot_point,2),'r*')
hold on
plot(-disp(plot_point,1),disp(plot_point,2),'r*')
hold on
plot(disp(plot_point,1),disp(plot_point,2),'r*')


c = [[disp(:,1);-disp(:,1)],[disp(:,2);disp(:,2)]]
save('disp_tend.dat', 'c', '-ascii', '-double', '-tabs')



axis equal
hold on


hold on
filename = strcat(pathSR,'srRod_',num2str(plotFiles(10)),'.csv');
disp = csvread(filename,1);
hold on
plot(disp(:,1),disp(:,2),'r');
hold on
plot(-disp(:,1),disp(:,2),'r');


c = [[disp(1:end-1,1);-disp(end:-1:1,1)],[disp(1:end-1,2);disp(end:-1:1,2)]]


save('disp_rod_tend.dat', 'c', '-ascii', '-double', '-tabs')




xlim([-50,50])
ylim([-300,ymax])



%% PRESSURE

m = dlmread('./../../build/bin/preform/pressureTime.txt',' ');
%m(ceil(length(m)/4):end,:) = smoothdata(m(ceil(length(m)/4):end,:))
figure
plot(m(:,1),m(:,2),'k-','linewidth',3);
xlabel('time');
ylabel('Pressure');
xlim([0,TMAX])
ylim([0,1])

m = m(1:10:end,:);

save('pressure_num.dat', 'm', '-ascii', '-double', '-tabs')

b = csvread("../preform/Experimental/N2_sr_pressure.csv");
hold on
plot(b(:,1)-0.008,b(:,2)/10,'b--')
b(:,2) = b(:,2)/10;


save('pressure_exp.dat', 'b', '-ascii', '-double', '-tabs')




legend('Mesh-free','Experiemntal')
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

legend({'Meshfree','Experimental'},... % { 'legend1', 'legend2',...}
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
    
    
    true_strain = logm(V);

    
    hoop_strain(i) =true_strain(3,3);
    axial_strain(i) = true_strain(2,2);
    radial_strain(i) = true_strain(1,1);
    shear_strain(i) = true_strain(1,2);
    


    


    
    
    jacobian(i) = det(F);
    crit_lambda(i) = strain(2,3);
    maxSr(i) = strain(3,1);
    max_lambda_n(i) = strain(1,3);
    time(i) = strain(3,2);
  
    
end


strain = figure('Position', get(0, 'Screensize'));
subplot(1,2,1)
plot(time,(axial_strain),'k','linewidth',2);
% hold on
% plot(time,log(1+shear_strain),'k');
xlim([0,TMAX])


cc = [time',axial_strain'];
save('strain_axial_num.dat', 'cc', '-ascii', '-double', '-tabs')

cc = [time',hoop_strain'];
save('strain_hoop_num.dat', 'cc', '-ascii', '-double', '-tabs')



c = csvread("Experimental/N2_sr_exp_axial.csv");


save('strain_axial_exp.dat', 'c', '-ascii', '-double', '-tabs')


hold on
plot(c(:,1)+0.01,c(:,2),'b--')



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

legend({'Axial (Meshfree)', 'Axial (exp)', 'Jacobian', 'Shear' },... % { 'legend1', 'legend2',...}
'interpreter','latex',...
'FontSize',12,...
'FontName','cmr14',...
'Location','SouthEast')

grid on 




subplot(1,2,2)
plot(time,(hoop_strain),'k','linewidth',2);


c = csvread("Experimental/N2_sr_exp_hoop.csv");
save('strain_hoop_exp.dat', 'c', '-ascii', '-double', '-tabs')

hold on
plot(c(:,1)+0.01,c(:,2),'b--')



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


% 
% 
% % legend
% path = './../../build/bin/preform/History/Stress';
% 
% 
% 
% 
% addpath(path)
% displacementdir = path ;
% d = dir(displacementdir);
% d1 = dir([displacementdir,'*.txt']);
% numFiles = size(d,1)/2 -3 ;
% 
% 
% for i = 1:length(plotFiles)
%     
%     
% 
%     filename = strcat(path,'/Bond_Stress_',num2str(plotFiles(i)),'.txt');
% 
%     stress = csvread(filename);
%     hoop_stress_bond(i) = stress(3,3);
%     axial_stress_bond(i) = stress(2,2);
%     radial_stress_bond(i) = stress(1,1);
%     shear_stress_bond(i) = stress(1,2);
% 
%     filename = strcat(path,'/Conformational_Stress_',num2str(plotFiles(i)),'.txt');
% 
%     stress = csvread(filename);
%     hoop_stress_conf(i) = stress(3,3);
%     axial_stress_conf(i) = stress(2,2);
%     radial_stress_conf(i) = stress(1,1);
%      shear_stress_conf(i) = stress(1,2);
%      
%      
%    volumetric_stress(i) = 1.8e9*log(jacobian(i)) ;
% 
%   
% end
% 
% figure
% subplot(1,2,1)
% hold on 
% plot(time,hoop_stress_conf,'r');
% hold on
% plot(time,hoop_stress_bond,'b');
% hold on
% plot(time,shear_stress_bond,'m')
% hold on
% plot(time,shear_stress_conf,'g')
% hold on
% plot(time,volumetric_stress,'y');
% %set axis
% set(gca, 'FontName', 'cmr12')
% % set x tics and y tics
% set(gca,...
% 'Units','normalized',...
% 'FontWeight','normal',...
% 'FontSize',14,... % size ofiguref numbers on axis
% 'FontName','cmr14') % font name
% set(gca,'TickLabelInterpreter', 'latex');
% 
% % Y label
% ylabel({'Cauchy Stress '},...
% 'interpreter','latex',...
% 'FontSize',14,... % font size
% 'FontName','cmr14')
% % X label
% xlabel('Time',...
% 'interpreter','latex',...
% 'FontWeight','normal',...
% 'FontSize',14,... % font size
% 'FontName','cmr14')
% % legend
% 
% legend({ 'Hoop (Conf)','Hoop (Bond)','Shear (Bond)', 'Shear (Conf)' },... % { 'legend1', 'legend2',...}
% 'interpreter','latex',...
% 'FontSize',12,...
% 'FontName','cmr14',...
% 'Location','SouthEast')
% 
% grid on 
% 
% subplot(1,2,2)
% plot(time,axial_stress_conf,'r');
% hold on
% plot(time,axial_stress_bond,'b');
% 
% %set axis
% set(gca, 'FontName', 'cmr12')
% % set x tics and y tics
% set(gca,...
% 'Units','normalized',...
% 'FontWeight','normal',...
% 'FontSize',14,... % size ofiguref numbers on axis
% 'FontName','cmr14') % font name
% set(gca,'TickLabelInterpreter', 'latex');
% 
% % Y label
% ylabel({'Cauchy Stress '},...
% 'interpreter','latex',...
% 'FontSize',14,... % font size
% 'FontName','cmr14')
% % X label
% xlabel('Time',...
% 'interpreter','latex',...
% 'FontWeight','normal',...
% 'FontSize',14,... % font size
% 'FontName','cmr14')
% % legend
% 
% legend({ 'Axial (Conf)','Axial (Bond)' },... % { 'legend1', 'legend2',...}
% 'interpreter','latex',...
% 'FontSize',12,...
% 'FontName','cmr14',...
% 'Location','SouthEast')
% 
% grid on 
% 
% 
% 
% print -dpng2 stress.png
% 
% figure 
% yyaxis left
% plot(time,crit_lambda);
% hold on
% yyaxis right
% plot(time,maxSr);
% hold on
% ax = exp(axial_strain)-1;
% hoop = exp(hoop_strain)-1;
% 
% diff_ax = diff(ax)./diff(time);
% diff_hoop = diff(hoop)./diff(time);
% hold on
% plot(time(2:end),diff_ax(1:end),'g-');
% hold on
% plot(time(2:end),diff_hoop(1:end),'m-');


 end

function [] = draw_ellipse(a,b,x0,y0)

    t = linspace(0,360,50);
    
    
    a = sqrt(a);
    b = sqrt(b);
    x = x0 + a*cosd(t);
    y = y0 + b*sind(t);
    
    hold on
    plot(x,y,'r-');
    
    
    
    


end

function [] = draw_general_ellipse(M,x0,y0)

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
    
    
    
    


end



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