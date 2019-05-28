% Draw Planar-Straight-Line-Graph of a preform
clc;
clear all;
close all;

%% RADII






[nodes,segments] = get_node_temperature(105,18);

plot(nodes(:,1),nodes(:,2),'b.')


axis off
axis equal

%boundary nodes


tri = csvread('./../../build/bin/preform/triangles.csv');

figure
p = patch('Faces',tri,'Vertices',nodes(:,1:2),'FaceVertexCData',nodes(:,3));
p.FaceColor = 'flat';
colorbar
axis equal
figure 

plot(nodes(:,1),nodes(:,2),'k.');
hold on
%plot(nodes(boundaryNodes(:,1),1),nodes(boundaryNodes(:,1),2),'b-')
axis off
axis equal
hold on
C = [];



% segments(2:(NUM_NODES_SPHERICAL_CAP + NUM_NODES_SIDEWALL + + NUM_NODES_RADIUS_TAPER_IN + NUM_NODES_TAPER -1 +NUM_NODES_TOP_FIXTURE ),3) = 2;
% count = NUM_NODES_SPHERICAL_CAP + NUM_NODES_SIDEWALL +NUM_NODES_TAPER + NUM_NODES_RADIUS_TAPER_IN - 2 + NUM_NODES_TOP_FIXTURE;
% segments(count:count+NUM_NODES_TOP_FIXTURE*1-1,3) = 5;
% count = count+NUM_NODES_TOP*1-1;
% segments(count:count + NUM_NODES_SIDEWALL+NUM_NODES_SPHERICAL_CAP + NUM_NODES_TAPER + 4 , 3) = 6;
% segments(end:-1:end-(NUM_NODES_BOT_FIXTURE-2),3) = 4;
% 
%    

%% plot segments
for i = 1:length(segments)
    if ( segments(i,3) == 2)
        color = 'r';
        hold on
        plot(nodes(segments(i,1:2),1),nodes(segments(i,1:2),2),color);
        
    elseif ( segments(i,3) == 5)
        color = 'g';
        hold on
        plot(nodes(segments(i,1:2),1),nodes(segments(i,1:2),2),color);
        
    elseif ( segments(i,3) == 0)
        color = 'r';
        hold on
        plot(nodes(segments(i,1:2),1),nodes(segments(i,1:2),2),color);
        
    elseif ( segments(i,3) == 4)
        color = 'g';
        hold on
        plot(nodes(segments(i,1:2),1),nodes(segments(i,1:2),2),color);
        
    else
        hold on
        plot(nodes(segments(i,1:2),1),nodes(segments(i,1:2),2),'m');
        
        
    end
    
    
end
% 


%% write files
dlmwrite('../../problems/preform/preform.nodes',[length(nodes),1],'delimiter',' ')
dlmwrite('../../problems/preform/preform.nodes',nodes,'-append','delimiter',' ')


%segments
dlmwrite('../../problems/preform/preform.segs',length(segments),'delimiter',' ')
dlmwrite('../../problems/preform/preform.segs',segments,'-append',.........
    'delimiter',' ')







