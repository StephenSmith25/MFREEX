
close all
clear all 

load ./Experimental/20190515_forStephen.mat

a = process_data(5);
b = a.outer_strain;
c = a.middle_strain;
time = a.time;
pressure = a.cavity_pressure;
dt = a.duration/size(b,3);

outer_coord = a.outer_Vic3D_coord;

t = 0;
T_OFFSET = 0.1005; %0.1005;

element = 84;





path = './../../build/bin/preform/Displacement/';
pathSR = './../../build/bin/preform/srRod/';
addpath(path)


displacementdir = path ;
d = dir(displacementdir);
d1 = dir([displacementdir,'*.csv']);
numFiles = size(d,1) -3 ;

plotFiles = ceil(linspace(1,numFiles,10));

plot_point =138; %%2241
filename = strcat(path,'displacement_',num2str(plotFiles(1)),'.csv');
disp = csvread(filename,1);


boundaryNodes = csvread('./../../build/bin/preform/boundary.txt');

boundaryNodes = [boundaryNodes;boundaryNodes(1)];



% -------------------------------------------------------------------------%
%                           PLOT 1
% -------------------------------------------------------------------------%



subplot(1,3,1)       % add first plot in 2 x 2 grid


plot(disp(:,1),disp(:,2),'k.','markersize',3)           % line plot
hold on
ymax = 0;
hold on
plot(disp(plot_point,1),disp(plot_point,2),'r*')




c = [[disp(:,1);-disp(:,1)],[disp(:,2);disp(:,2)]]








save('disp_tstart.dat', 'c', '-ascii', '-double', '-tabs')



axis equal 
hold on
plot(disp(boundaryNodes,1),disp(boundaryNodes,2),'k-')




c = [[disp(1:end-1,1);-disp(end:-1:1,1)],[disp(1:end-1,2);disp(end:-1:1,2)]];


save('disp_rod_tstart.dat', 'c', '-ascii', '-double', '-tabs')





hold on
xlim([0,50])
ylim([-350,ymax])

outer_coords_i = outer_coord(:,:,1);
hold on 
plot(outer_coords_i(:,1),outer_coords_i(:,2),'b','linewidth',2);
hold on 
plot(outer_coords_i(element,1),outer_coords_i(element,2),'g*');


% -------------------------------------------------------------------------%
%                           PLOT 2
% -------------------------------------------------------------------------%




filename = strcat(path,'displacement_',num2str(plotFiles(10)),'.csv');
disp = csvread(filename,1);
subplot(1,3,2)       % add first plot in 2 x 2 grid




plot(disp(:,1),disp(:,2),'k.','markersize',3)           % line plot
hold on
plot(disp(plot_point,1),disp(plot_point,2),'r*')
hold on
plot(disp(boundaryNodes,1),disp(boundaryNodes,2),'k-')

c = [[disp(:,1);-disp(:,1)],[disp(:,2);disp(:,2)]]
save('disp_tmid.dat', 'c', '-ascii', '-double', '-tabs')




axis equal




c = [[disp(1:end-1,1);-disp(end:-1:1,1)],[disp(1:end-1,2);disp(end:-1:1,2)]];

outer_coords_i = outer_coord(:,:,floor(17*end/32));
hold on 
plot(outer_coords_i(:,1),outer_coords_i(:,2),'b','linewidth',2);
hold on 
plot(outer_coords_i(element,1),outer_coords_i(element,2),'g*');





xlim([-0,50])
ylim([-350,ymax])


% -------------------------------------------------------------------------%
%                           PLOT 3
% -------------------------------------------------------------------------%
filename = strcat(path,'displacement_',num2str(plotFiles(10)),'.csv');
disp = csvread(filename,1);




subplot(1,3,3)       % add first plot in 2 x 2 grid


plot(disp(:,1),disp(:,2),'k.','markersize',5)           % line plot
hold on
plot(disp(boundaryNodes,1),disp(boundaryNodes,2),'k-')

hold on
plot(disp(plot_point,1),disp(plot_point,2),'r*')


c = [[disp(:,1);-disp(:,1)],[disp(:,2);disp(:,2)]];
save('disp_tend.dat', 'c', '-ascii', '-double', '-tabs')



axis equal
hold on



c = [[disp(1:end-1,1);-disp(end:-1:1,1)],[disp(1:end-1,2);disp(end:-1:1,2)]];

outer_coords_i = outer_coord(:,:,floor(29*end/32));
hold on 
plot(outer_coords_i(:,1),outer_coords_i(:,2),'b','linewidth',2);
hold on 
plot(outer_coords_i(element,1),outer_coords_i(element,2),'g*');





xlim([0,60])
ylim([-350,ymax])


% ylim([-250,0])
% 
% ix = find(time > T_OFFSET);
% time = time(ix)-T_OFFSET;
% pressure = pressure(ix);
% 
% figure
% for i = 1:size(b,3)
%     strain_outer = b(:,:,i);
%     strain_middle = c(:,:,i);
%     
%     hold on
%     
%     
%     plot(t-T_OFFSET,log(strain_outer(element,1)+1),'bx');
%     hold on
%     plot(t-T_OFFSET,log(strain_outer(element,2)+1),'g.');
%     
%     y_axial(i,:) = [t-T_OFFSET,log(strain_outer(element,2)+1)];
%     y_hoop(i,:) = [t-T_OFFSET,log(strain_outer(element,1)+1)];
% 
%     t = t + dt;
% end
% xlim([0,0.35])
% 
% 
% csvwrite("./Excel/N5_sr_exp_axial.csv",y_axial);
% csvwrite("./Excel/N5_sr_exp_hoop.csv",y_hoop);
% figure
% 
% yyaxis left
% plot(time,pressure*100);
% xlim([0,0.35]);
% 
% sr = a.internal_force(ix);
% 
% hold on
% yyaxis right
% plot(time,sr);
% 
% ab = diff(pressure*100)./diff(time);
% iy = find(time < 0.01);
% 
% pressure_grad = mean(ab(iy));
% 
% y = [time,pressure];
% csvwrite("./Excel/N5_sr_pressure.csv",y);





% figure
% num_loops = size(outer_coord,3);
% F(num_loops) = struct('cdata',[],'colormap',[]);
% 
% for i = 1:num_loops
%     
%     outer_coords_i = outer_coord(:,:,i);
% plot(outer_coords_i(:,1),outer_coords_i(:,2),'b.');
%     axis equal
% 
%    ylim([-250,0])
%    xlim([0,60])
%     drawnow
%     F(i) = getframe(gcf);
% 
%     
% end
% 
% v = VideoWriter('newfile.avi','Motion JPEG AVI');
% v.Quality = 100;
% open(v)
% writeVideo(v,F)
% close(v)
