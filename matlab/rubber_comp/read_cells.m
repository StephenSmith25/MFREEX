clear all
close all 


nodes = csvread('./../../build/bin/rubber_comp/nodes.csv');


tri = csvread('./../../build/bin/rubber_comp/triangles.csv');


figure
triplot(tri,nodes(:,1),nodes(:,2));

boundary_nodes = csvread('./../../build/bin/rubber_comp/boundary.txt');
boundary_nodes = [boundary_nodes;boundary_nodes(1)];





fileID = fopen("./../../build/bin/rubber_comp/cells.txt");


tline = fgetl(fileID);
C = strsplit(tline);
num_verticies = str2double(C{1});
dim = str2double(C{2});

verticies = zeros(num_verticies,dim);

for i = 1:num_verticies


    tline = fgetl(fileID);
    C = strsplit(tline);
    
    for k = 1:dim
        verticies(i,k) = str2double(C{k});
    end
    
        
    
end


figure
plot(verticies(:,1),verticies(:,2),'b*')

axis equal

tline = fgetl(fileID);
C = strsplit(tline);
num_cells = str2double(C{1});

figure

for i = 1:num_cells

    tline = fgetl(fileID);
    C = strsplit(tline);
    
    num_cell_verticies = str2double(C{1});
    
    poly = zeros(num_cell_verticies,2);
    for k = 1:num_cell_verticies
        poly(k,1) = verticies(str2double(C{1+k})+1,1);
        poly(k,2) = verticies(str2double(C{1+k})+1,2);

    end
    area(i) = polyarea(poly(:,1),poly(:,2));
    
    hold on 
    fill(poly(:,1),poly(:,2),rand(1,3));
end
axis equal
% path_base = './../../build/bin/cylinder/';


%material_points = csvread(filename,1);

% path = './../../build/bin/cylinder/Displacement';
% 
% filename = strcat(path,'/displacement_','1','.txt');
% disp = nodes;
% plot(disp(:,1),disp(:,2),'kd','markersize',5,...
%     'MarkerEdgeColor','black',...
%     'MarkerFaceColor',[1 1 1])    
% xlim([0,22])
% ylim([0,22])
% 
% 
% 
% % Y label
% ylabel({'y'},...
% 'interpreter','latex',...
% 'FontSize',14,... % font size
% 'FontName','cmr14')
% % X label
% xlabel('x',...
% 'interpreter','latex',...
% 'FontSize',14,... % font size
% 'FontName','cmr14')
% 
% 
% 
% set(gcf,'color','white');
% 
% 
% 
% export_fig voronoi_cylinder.png -m8
% 
% 
% 
% 
% figure
% triangles = csvread('./../../build/bin/cylinder/triangles.csv');

% 
% for i = 1:length(triangles)
%     
%    poly = disp(triangles(i,:),:) ;
%    
%        
%     hold on 
%     fill(poly(:,1),poly(:,2),rand(1,3));
%    
% end
% 
% %triplot(triangles,disp(:,1),disp(:,2));
% hold on
% plot(disp(:,1),disp(:,2),'kd','markersize',5,...
%     'MarkerEdgeColor','black',...
%     'MarkerFaceColor',[1 1 1])    
% 
% 
% axis equal
% 
% 
% 
% xlim([0,22])
% ylim([0,22])
% 
% % Y label
% ylabel({'y'},...
% 'interpreter','latex',...
% 'FontSize',14,... % font size
% 'FontName','cmr14')
% % X label
% xlabel('x',...
% 'interpreter','latex',...
% 'FontSize',14,... % font size
% 'FontName','cmr14')
% 
% 
% 
% set(gcf,'color','white');
% 
% 
% 
% export_fig triangles_cylinder.png -m8
% 
% while ischar(tline)
%     disp(tline)
%     [token,remain] = strtok(tline)
%     while ~isempty(remain)
%         
%         [token,remain] = strtok(remain)
%         
%     end
%     tline = fgetl(fileID);
% end