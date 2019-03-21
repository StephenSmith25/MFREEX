
close all
clear all





material_points = csvread('./../../build/bin/cylinder/materialpoints.csv');



plot(material_points(:,1),material_points(:,2),'r*');

%axis equal 


nodes = csvread('./../../build/bin/cylinder/nodes_2.csv');



%tri = csvread('./../../build/bin/cylinder/triangles.csv');

%triplot(tri,nodes(:,1),nodes(:,2));


cells = csvread('./../../build/bin/cylinder/search_cells.csv');


active_cells =  csvread('./../../build/bin/cylinder/active_cells.csv');
nodes_origin = csvread('./../../build/bin/cylinder/nodes.csv');


for k = 1:length(active_cells)
hold on
    i = active_cells(k)+1;
rectangle('Position',[cells(i,1), cells(i,3), cells(i,2)-cells(i,1), cells(i,4)-cells(i,3)])
end

plot(nodes_origin(:,1),nodes_origin(:,2),'b.');

% 
% v = VideoWriter('cells.avi');
% v.FrameRate = 5;
% open(v);
% for i = 1:49
% 
%     filename = strcat('./../../build/bin/cylinder/nodes_',num2str(i),'.csv');
%    
% 
%     
%     
%     nodes = csvread(filename);
%     hold off
%     plot(nodes(:,1),nodes(:,2),'b.');
%     
%     hold on
%     plot(nodes_origin(:,1),nodes_origin(:,2),'r.');
% 
%     filename = strcat('./../../build/bin/cylinder/search_cells_',num2str(i),'.csv');
%     
%     active_cells = csvread(filename);
%     
%     
%     for k = 1:length(active_cells)
%         hold on
%         l = active_cells(k)+1;
%         rectangle('Position',[cells(l,1), cells(l,3), cells(l,2)-cells(l,1), cells(l,4)-cells(l,3)])
%     end
% 
%     
%    
% ylim([0,35]) 
%     xlim([0,35])
% 
%     
%     drawnow
%    
%     F(i) = getframe(gcf);
%     writeVideo(v,F(i));
% 
% end
% 
% 
% close(v);
% 
% 
% active_nodes = csvread("./../../build/bin/cylinder/active_nodes.csv");
% 
% nodes_active = [];
% for i = 1:length(active_nodes)
%     
%     
%    for j = 1:size(active_nodes,2)
%         ix = active_nodes(i,j);
%         
%         if ( ix > 0 )
%            nodes_active = [nodes_active ; ix]; 
%         end
%        
%    end
%    
% 
% %   figure
% %  ylim([0,35]) 
% %     xlim([0,35])
% %     ax = gca;
% %    
% % movie(gcf,F)
% 
% end
% nodes_active = sort(nodes_active);