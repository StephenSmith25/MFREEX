
clear all
close all

TMAX = 0.3;


set(gcf, 'Color', 'w');

path = './../../build/bin/preform/History/Strain/Top';

addpath(path)
displacementdir = path ;
d = dir(displacementdir);
d1 = dir([displacementdir,'*.txt']);
numFiles = size(d,1) -3 ;

plotFiles = ceil(linspace(1,numFiles,100));
for i = 1:length(plotFiles)
    
    

    filename = strcat('./../../build/bin/preform/History/Strain/Top/strain_',num2str(plotFiles(i)),'.txt');

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
subplot(1,3,1)
plot(time,(axial_strain),'k','linewidth',2);
% hold on
% plot(time,log(1+shear_strain),'k');
xlim([0,TMAX])


cc = [time',axial_strain'];
save('N5_strain_axial_num_top.dat', 'cc', '-ascii', '-double', '-tabs')

cc = [time',hoop_strain'];
save('N5_strain_hoop_num_top.dat', 'cc', '-ascii', '-double', '-tabs')



c = csvread("Experimental/N5_sr_exp_axial_TOP.csv");


save('strain_axial_exp.dat', 'c', '-ascii', '-double', '-tabs')


hold on
plot(c(:,1),c(:,2),'b--')





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

plot(time,(hoop_strain),'k','linewidth',2);


c = csvread("Experimental/N5_sr_exp_hoop_TOP.csv");
save('strain_hoop_exp.dat', 'c', '-ascii', '-double', '-tabs')

hold on
plot(c(:,1),c(:,2),'b--')



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
title('Top')

grid on 


path = './../../build/bin/preform/History/Strain/Middle';

addpath(path)
displacementdir = path ;
d = dir(displacementdir);
d1 = dir([displacementdir,'*.txt']);
numFiles = size(d,1) -3 ;

plotFiles = ceil(linspace(1,numFiles,100));
for i = 1:length(plotFiles)
    
    

    filename = strcat('./../../build/bin/preform/History/Strain/Middle/strain_',num2str(plotFiles(i)),'.txt');

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

subplot(1,3,2)
plot(time,(axial_strain),'k','linewidth',2);
% hold on
% plot(time,log(1+shear_strain),'k');
xlim([0,TMAX])


cc = [time',axial_strain'];
save('N5_strain_axial_num_middle.dat', 'cc', '-ascii', '-double', '-tabs')

cc = [time',hoop_strain'];
save('N5_strain_hoop_num_middle.dat', 'cc', '-ascii', '-double', '-tabs')



c = csvread("Experimental/N5_sr_exp_axial.csv");


save('strain_axial_exp.dat', 'c', '-ascii', '-double', '-tabs')


hold on
plot(c(:,1),c(:,2),'b--')




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

hold on
plot(time,(hoop_strain),'k','linewidth',2);


c = csvread("Experimental/N5_sr_exp_hoop_MIDDLE.csv");
save('strain_hoop_exp.dat', 'c', '-ascii', '-double', '-tabs')

hold on
plot(c(:,1),c(:,2),'b--')



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

title('Middle')



path = './../../build/bin/preform/History/Strain/Bot';

addpath(path)
displacementdir = path ;
d = dir(displacementdir);
d1 = dir([displacementdir,'*.txt']);
numFiles = size(d,1) -3 ;

plotFiles = ceil(linspace(1,numFiles,100));
for i = 1:length(plotFiles)
    
    

    filename = strcat('./../../build/bin/preform/History/Strain/Bot/strain_',num2str(plotFiles(i)),'.txt');

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

subplot(1,3,3)
plot(time,(axial_strain),'k','linewidth',2);
% hold on
% plot(time,log(1+shear_strain),'k');
xlim([0,TMAX])


cc = [time',axial_strain'];
save('N5_strain_axial_num_bot.dat', 'cc', '-ascii', '-double', '-tabs')

cc = [time',hoop_strain'];
save('N5_strain_hoop_num_bot.dat', 'cc', '-ascii', '-double', '-tabs')



c = csvread("Experimental/N5_sr_exp_axial_BOT.csv");


save('strain_axial_exp.dat', 'c', '-ascii', '-double', '-tabs')


hold on
plot(c(:,1),c(:,2),'b--')



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

plot(time,(hoop_strain),'k','linewidth',2);


c = csvread("Experimental/N5_sr_exp_hoop_BOT.csv");
save('strain_hoop_exp.dat', 'c', '-ascii', '-double', '-tabs')

hold on
plot(c(:,1),c(:,2),'b--')



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

title('Bot')

grid on 
