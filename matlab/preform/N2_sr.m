
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
element_top = 18;
element_middle = 58;
element_bot = 95;


figure

outer_coords_i = outer_coord(:,:,floor(3*end/4));

% 
plot(outer_coords_i(:,1),outer_coords_i(:,2),'b.');
hold on
outer_coords_i = outer_coord(:,:,1);
plot(outer_coords_i(:,1),outer_coords_i(:,2),'b.');
hold on
plot(outer_coords_i(element_top,1),outer_coords_i(element_top,2),'r*');
hold on
plot(outer_coords_i(element_middle,1),outer_coords_i(element_middle,2),'r*');
hold on
plot(outer_coords_i(element_bot,1),outer_coords_i(element_bot,2),'r*');





hold on
outer_coords_i = outer_coord(:,:,floor(end/3));
plot(outer_coords_i(:,1),outer_coords_i(:,2),'b.');
hold on
plot(outer_coords_i(element_top,1),outer_coords_i(element_top,2),'r*');
hold on
plot(outer_coords_i(element_middle,1),outer_coords_i(element_middle,2),'r*');
hold on
plot(outer_coords_i(element_bot,1),outer_coords_i(element_bot,2),'r*');



axis equal
ylim([-250,0])

ix = find(time > T_OFFSET);
time = time(ix)-T_OFFSET;
pressure = pressure(ix);


t=0;
figure

y_axial = [];
y_hoop = [];


for i = 1:size(b,3)
    strain_outer = b(:,:,i);
    strain_middle = c(:,:,i);
    
    hold on
    
    
    plot(t-T_OFFSET,log(strain_outer(element_top,1)+1),'bx');
    hold on
    plot(t-T_OFFSET,log(strain_outer(element_top,2)+1),'g.');
    
    y_axial(i,:) = [t-T_OFFSET,log(strain_outer(element_top,2)+1)];
    y_hoop(i,:) = [t-T_OFFSET,log(strain_outer(element_top,1)+1)];

    t = t + dt;
end




csvwrite("./Experimental/N2_sr_exp_axial_TOP.csv",y_axial);
csvwrite("./Experimental/N2_sr_exp_hoop_TOP.csv",y_hoop);
y_axial = [];
y_hoop = [];



t=0;
for i = 1:size(b,3)
    strain_outer = b(:,:,i);
    strain_middle = c(:,:,i);
    
    hold on
    
    
    plot(t-T_OFFSET,log(strain_outer(element_middle,1)+1),'bx');
    hold on
    plot(t-T_OFFSET,log(strain_outer(element_middle,2)+1),'g.');
    
    y_axial(i,:) = [t-T_OFFSET,log(strain_outer(element_middle,2)+1)];
    y_hoop(i,:) = [t-T_OFFSET,log(strain_outer(element_middle,1)+1)];

    t = t + dt;
end
xlim([0,0.35])
csvwrite("./Experimental/N2_sr_exp_axial_MIDDLE.csv",y_axial);
csvwrite("./Experimental/N2_sr_exp_hoop_MIDDLE.csv",y_hoop);



t=0;

y_axial = [];
y_hoop = [];


for i = 1:size(b,3)
    strain_outer = b(:,:,i);
    strain_middle = c(:,:,i);
    
    hold on
    
    
    plot(t-T_OFFSET,log(strain_outer(element_bot,1)+1),'bx');
    hold on
    plot(t-T_OFFSET,log(strain_outer(element_bot,2)+1),'g.');
    
    y_axial(i,:) = [t-T_OFFSET,log(strain_outer(element_bot,2)+1)];
    y_hoop(i,:) = [t-T_OFFSET,log(strain_outer(element_bot,1)+1)];

    t = t + dt;
end
xlim([0,0.35])


csvwrite("./Experimental/N2_sr_exp_axial_BOT.csv",y_axial);
csvwrite("./Experimental/N2_sr_exp_hoop_BOT.csv",y_hoop);



figure

yyaxis left
plot(time,pressure*100);
xlim([0,0.35]);

sr = a.internal_force(ix);

hold on
yyaxis right
plot(time,sr);

ab = diff(pressure*100)./diff(time);
iy = find(time < 0.01);

pressure_grad = mean(ab(iy));

y = [time,pressure];
csvwrite("./Experimental/N2_sr_pressure.csv",y);





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
