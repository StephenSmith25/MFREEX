close all
clear all
clc

W = 2.5;
L = 20;
rin =60;

% nodes
nodes = [];

count =1;
nodes(count,:) = [0,0];
count = count + 1;
nodes(count,:) = [0,W];
count = count + 1;
nodes(count,:) = [L,W];
count = count + 1;
nodes(count,:) = [L,0];
count = count + 1;

nodes(:,1) = nodes(:,1) + rin;
segments = [];
% Write the segments 
for i = 1:length(nodes)

    if ( i == length(nodes))
        segments(i,:) = [i,1,i];
    else
        segments(i,:) = [i,i+1,i]; 
    end
        
    
end
segments(1,3) = 2;
segments(2,3) = 5;
segments(3,3) = 0;
segments(4,3) = 4;


% Plot it 
figure
plot(nodes(:,1),nodes(:,2),'k-');
axis off
axis equal 
    rand1 = rand(1,3);
hold on 
fill(nodes(:,1),nodes(:,2),rand1);
hold on 
plot(nodes(:,1),nodes(:,2),'k-');


%% write files
dlmwrite('square.nodes',length(nodes),'delimiter',' ')
dlmwrite('square.nodes',nodes,'-append','delimiter',' ')

%segments
dlmwrite('square.segs',length(segments),'delimiter',' ')
dlmwrite('square.segs',segments,'-append',.........
    'delimiter',' ')

