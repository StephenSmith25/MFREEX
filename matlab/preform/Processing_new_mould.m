
clear all
close all

PLOT_GRAPHS = true;
PLOT_DOMAINS_INFLUENCE = true; 


WITH_MOULD = false;
TMAX = 0.30;




path = './../../build/bin/preform/Displacement/';
pathSR = './../../build/bin/preform/srRod/';
boundaryNodes = csvread('./../../build/bin/preform/boundary.txt');
path_base = './../../build/bin/preform/';

boundaryNodes = [boundaryNodes;boundaryNodes(1)];



segments = [];
for i = 1:length(boundaryNodes)
    segments(i,:) = [i,i+1];
    
    if ( i == length(boundaryNodes))
        
    segments(i,:) = [i,1];
    end
end
        


addpath(path)
displacementdir = path ;
d = dir(displacementdir);
d1 = dir([displacementdir,'*.csv']);
numFiles = size(d,1) -3 ;

plotFiles = ceil(linspace(1,numFiles,10));

plot_point =155; %%2241
filename = strcat(path,'displacement_',num2str(plotFiles(1)),'.csv');
disp = csvread(filename,1);





 [mould_nodes] = draw_mould_ribs();
ymax = 0;

figure



% -------------------------------------------------------------------------%
%                           PLOT 1
% -------------------------------------------------------------------------%



subplot(1,3,1)       % add first plot in 2 x 2 grid

hold on 
hold on 
plot(mould_nodes(:,1),mould_nodes(:,2),'k');
hold on
plot(-mould_nodes(:,1),mould_nodes(:,2),'k');
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

boundarydisp = disp(boundaryNodes,:);

hold on
filename = strcat(pathSR,'srRod_',num2str(plotFiles(1)),'.csv');
disp = csvread(filename,1);
hold on
plot(disp(:,1),disp(:,2),'r');
hold on
plot(-disp(:,1),disp(:,2),'r');


c = [[disp(1:end-1,1);-disp(end:-1:1,1)],[disp(1:end-1,2);disp(end:-1:1,2)]]


save('disp_rod_tstart.dat', 'c', '-ascii', '-double', '-tabs')



write_triangle_file(boundarydisp,segments,'bottle.poly');
system('triangle -pDq30 bottle.poly')
convert_tri_file_to_vtk('bottle.1','bottle_tstart.vtk',0)



hold on
xlim([-50,50])
ylim([-350,ymax])



% -------------------------------------------------------------------------%
%                           PLOT 2
% -------------------------------------------------------------------------%




filename = strcat(path,'displacement_',num2str(plotFiles(3)),'.csv');
disp = csvread(filename,1);
subplot(1,3,2)       % add first plot in 2 x 2 grid
hold on 
plot(mould_nodes(:,1),mould_nodes(:,2),'k');
hold on
plot(-mould_nodes(:,1),mould_nodes(:,2),'k');
hold on



plot(disp(:,1),disp(:,2),'k.','markersize',3)           % line plot
hold on
plot(-disp(:,1),disp(:,2),'k.','markersize',3) % line plot
hold on
hold on
plot(-disp(plot_point,1),disp(plot_point,2),'r*')
hold on
plot(disp(plot_point,1),disp(plot_point,2),'r*')


c = [[disp(:,1);-disp(:,1)],[disp(:,2);disp(:,2)]];
save('disp_tmid.dat', 'c', '-ascii', '-double', '-tabs')


boundarydisp = disp(boundaryNodes,:);


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


c = [[disp(1:end-1,1);-disp(end:-1:1,1)],[disp(1:end-1,2);disp(end:-1:1,2)]];



write_triangle_file(boundarydisp,segments,'bottle.poly');
system('triangle -pDq30 bottle.poly')
convert_tri_file_to_vtk('bottle.1','bottle_tmid.vtk',0)

save('disp_rod_tmid.dat', 'c', '-ascii', '-double', '-tabs')





xlim([-50,50])
ylim([-350,ymax])



print -dpng2 Displacement.png

% -------------------------------------------------------------------------%
%                           PLOT 3
% -------------------------------------------------------------------------%
filename = strcat(path,'displacement_',num2str(plotFiles(10)),'.csv');
disp = csvread(filename,1);




subplot(1,3,3)       % add first plot in 2 x 2 grid

hold on 
plot(mould_nodes(:,1),mould_nodes(:,2),'k');
hold on
plot(-mould_nodes(:,1),mould_nodes(:,2),'k');
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





    

% write the triangle file

boundarydisp = disp(boundaryNodes,:);
write_triangle_file(boundarydisp,segments,'bottle.poly');
system('triangle -pDq30 bottle.poly')
convert_tri_file_to_vtk('bottle.1','bottle_tend.vtk',0)




hold on
plot(-disp(plot_point,1),disp(plot_point,2),'r*')
hold on
plot(disp(plot_point,1),disp(plot_point,2),'r*')


c = [[disp(:,1);-disp(:,1)],[disp(:,2);disp(:,2)]];
save('disp_tend.dat', 'c', '-ascii', '-double', '-tabs')

csvwrite('paraview_test.csv',disp(boundaryNodes,:)) ;
csvwrite('mould_test.csv',mould_nodes) ;

axis equal
hold on




hold on
filename = strcat(pathSR,'srRod_',num2str(plotFiles(10)),'.csv');
disp = csvread(filename,1);
hold on
plot(disp(:,1),disp(:,2),'r');
hold on
plot(-disp(:,1),disp(:,2),'r');


c = [[disp(1:end-1,1);-disp(end:-1:1,1)],[disp(1:end-1,2);disp(end:-1:1,2)]];


save('disp_rod_tend.dat', 'c', '-ascii', '-double', '-tabs')




xlim([-60,60])
ylim([-350,ymax])





function [] = write_triangle_file(X,segments,filename)


fileID = fopen(filename,'w');

fprintf(fileID,'%d %d %d \n',length(X),2,0);

for i = 1:length(X)
   fprintf(fileID,'%d %f %f\n',i,X(i,1),X(i,2));
    
end
fprintf(fileID,'%d \n',length(segments));

for i = 1:length(segments)
   fprintf(fileID,'%d %d %d\n',i,segments(i,1),segments(i,2));
    
end
fprintf(fileID,'%d \n',0);

fclose(fileID);



end



