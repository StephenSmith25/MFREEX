% Draw Planar-Straight-Line-Graph of a preform
clc;
clear all;
close all;

%% RADII


% BOTTOM SPHERICAL RADIUS
ROUT_BOT = 11.2522;
RIN_BOT = 6.986;


% LENGTHS
LENGTH_PREFORM = 77.6986;
LENGTH_NECK = 3;
LENGTH_TAPER = 17.8308;
LENGTH_SIDEWALL = LENGTH_PREFORM - LENGTH_NECK - LENGTH_TAPER - ROUT_BOT;

% BOTTOM CYLINDRICAL RADIUS

RIN_BOT_SIDEWALL = 13.97/2;
ROUT_BOT_SIDEWALL = 22.5044/2;


%TOP CYLIDER RADIUS ( BELOW NECK)
RIN_TOP_SIDEWALL = 14.8844/2;
ROUT_TOP_SIDEWALL = 23.1902/2;


% TOP NECK
RIN_NECK = 24.3078/2;
ROUT_NECK = 29.4894/2;


THICKNESS_BOT_CAP = 3.2;

CIRCLE_BOT_OFFSET = (ROUT_BOT - RIN_BOT) - THICKNESS_BOT_CAP;

% number of nodes through the thickness
NUM_NODES_THICKNESS = 4;



%% NUMBER OF NODES
NUM_NODES_SIDEWALL =56;
NUM_NODES_TAPER = 18;
NUM_NODES_SPHERICAL_CAP = 20;
NUM_NODES_TOP_FIXTURE = 5;
NUM_NODES_BOT_FIXTURE = 6;




% starting from 0,0
nodes = [];
count = 1;

theta = linspace(-90,0,NUM_NODES_SPHERICAL_CAP);

ymax = -CIRCLE_BOT_OFFSET;



for i = 1:length(theta)
    nodes(count,:) = [RIN_BOT*cosd(theta(i)), -CIRCLE_BOT_OFFSET +  RIN_BOT*sind(theta(i))];
    count = count +1;
    
end




% Left wall
%% traction nodes
nodes(count:1:count+NUM_NODES_SIDEWALL-1,:) = [linspace(RIN_BOT_SIDEWALL,RIN_TOP_SIDEWALL,NUM_NODES_SIDEWALL)',linspace(0,LENGTH_SIDEWALL,NUM_NODES_SIDEWALL)'];
count = count + NUM_NODES_SIDEWALL;
nodes(count:1:count+NUM_NODES_TAPER-1,:) = [linspace(RIN_TOP_SIDEWALL,RIN_NECK,NUM_NODES_TAPER)',linspace(LENGTH_SIDEWALL,LENGTH_SIDEWALL+LENGTH_TAPER,NUM_NODES_TAPER)'];
count = count + NUM_NODES_TAPER;
%% Essential boundary nodes
nodes(count:1:count+NUM_NODES_TOP_FIXTURE-1,:) = [linspace(RIN_NECK,RIN_NECK,NUM_NODES_TOP_FIXTURE)',linspace(LENGTH_SIDEWALL+LENGTH_TAPER,LENGTH_SIDEWALL + LENGTH_TAPER + LENGTH_NECK,NUM_NODES_TOP_FIXTURE)'];
count = count + NUM_NODES_TOP_FIXTURE;
height_top = max(nodes(:,2));

nodes(count:1:count+NUM_NODES_TOP_FIXTURE-1,:) = [linspace(RIN_NECK,ROUT_NECK,NUM_NODES_TOP_FIXTURE)',linspace(height_top,height_top,NUM_NODES_TOP_FIXTURE)'];
count = count + NUM_NODES_TOP_FIXTURE;

nodes(count:1:count+NUM_NODES_TOP_FIXTURE-1,:) = [linspace(ROUT_NECK,ROUT_NECK,NUM_NODES_TOP_FIXTURE)',linspace(height_top,LENGTH_SIDEWALL + LENGTH_TAPER,NUM_NODES_TOP_FIXTURE)'];
count = count + NUM_NODES_TOP_FIXTURE;
%% Outside nodes
nodes(count:1:count+NUM_NODES_TAPER-1,:) = [linspace(ROUT_NECK,ROUT_TOP_SIDEWALL,NUM_NODES_TAPER)',linspace(LENGTH_SIDEWALL+LENGTH_TAPER,LENGTH_SIDEWALL,NUM_NODES_TAPER)'];
count = count + NUM_NODES_TAPER;
nodes(count:1:count+NUM_NODES_SIDEWALL-1,:) = [linspace(ROUT_TOP_SIDEWALL,ROUT_BOT_SIDEWALL,NUM_NODES_SIDEWALL)',linspace(LENGTH_SIDEWALL,0,NUM_NODES_SIDEWALL)'];
count = count + NUM_NODES_SIDEWALL;

theta = linspace(0,-90,NUM_NODES_SPHERICAL_CAP);

for i = 1:length(theta)
    nodes(count,:) = [ROUT_BOT*cosd(theta(i)), ROUT_BOT*sind(theta(i))];
    count = count +1;
    
end
% Bottom edge
nodes(count:1:count+NUM_NODES_BOT_FIXTURE-1,:) = [linspace(0,0,NUM_NODES_BOT_FIXTURE)',linspace(-ROUT_BOT,-RIN_BOT-CIRCLE_BOT_OFFSET,NUM_NODES_BOT_FIXTURE)'];

[~,ib] = unique(nodes,'rows');

ib = sort(ib);
nodes = nodes(ib,:);




%boundary nodes

boundaryNodes = linspace(1,length(nodes),length(nodes))';
boundaryNodes(1:(NUM_NODES_SPHERICAL_CAP + NUM_NODES_SIDEWALL +NUM_NODES_TAPER -1  ),2) = 2;
count = NUM_NODES_SPHERICAL_CAP + NUM_NODES_SIDEWALL +NUM_NODES_TAPER -1 ;
boundaryNodes(count:count+NUM_NODES_TOP_FIXTURE*3-3,2) = 5;
count = count+NUM_NODES_TOP_FIXTURE*3-2;
boundaryNodes(count:count + NUM_NODES_SIDEWALL+NUM_NODES_SPHERICAL_CAP + NUM_NODES_TAPER , 2) = 6;

boundaryNodes(end:-1:end-(NUM_NODES_BOT_FIXTURE-2),2) = 4;
% 


%  Through thickness nodes
nodes1 = [];
count = 1;
% 
RADII_BOT = linspace(RIN_BOT,ROUT_BOT,NUM_NODES_THICKNESS+2);
RADII_BOT_SIDEWALL = linspace(RIN_BOT_SIDEWALL,ROUT_BOT_SIDEWALL,NUM_NODES_THICKNESS+2);
CIRCLE_BOT_OFFSETS = linspace(CIRCLE_BOT_OFFSET, 0, NUM_NODES_THICKNESS+2);
RADII_TOP_SIDEWALL = linspace(RIN_TOP_SIDEWALL,ROUT_TOP_SIDEWALL,NUM_NODES_THICKNESS+2);
RADII_NECK = linspace(RIN_NECK,ROUT_NECK,NUM_NODES_THICKNESS+2);

nodes1 = [];
count = 1;



 for i = 1:NUM_NODES_THICKNESS

     

theta_max = asind((ymax + CIRCLE_BOT_OFFSETS(i+1))/RADII_BOT(i+1));
theta = linspace(-90,theta_max,NUM_NODES_SPHERICAL_CAP);



for j = 1:length(theta)-1
    nodes1(count,:) = [RADII_BOT(i+1)*cosd(theta(j+1)),  RADII_BOT(i+1)*sind(theta(j+1)) - CIRCLE_BOT_OFFSETS(i+1)   ];
    count = count +1;
    
end




%Left wall
%%traction nodes
nodes1(count:1:count+NUM_NODES_SIDEWALL-1,:) = [linspace(RADII_BOT_SIDEWALL(i+1),RADII_TOP_SIDEWALL(i+1),NUM_NODES_SIDEWALL)',linspace(0,LENGTH_SIDEWALL,NUM_NODES_SIDEWALL)'];
 count = count + NUM_NODES_SIDEWALL;
 nodes1(count:1:count+NUM_NODES_TAPER-1,:) = [linspace(RADII_TOP_SIDEWALL(i+1),RADII_NECK(i+1),NUM_NODES_TAPER)',linspace(LENGTH_SIDEWALL,LENGTH_SIDEWALL+LENGTH_TAPER,NUM_NODES_TAPER)'];
 count = count + NUM_NODES_TAPER;
% %% Essential boundary nodes

heights = linspace(LENGTH_SIDEWALL+LENGTH_TAPER,LENGTH_SIDEWALL + LENGTH_TAPER + LENGTH_NECK  ,NUM_NODES_TOP_FIXTURE)';

nodes1(count:1:count+NUM_NODES_TOP_FIXTURE-2,:) = [linspace(RADII_NECK(i+1),RADII_NECK(i+1),NUM_NODES_TOP_FIXTURE-1)', heights(1:end-1)];
count = count + NUM_NODES_TOP_FIXTURE-1;
height_top = max(nodes(:,2));





 end
  
[~,ib] = unique(nodes1,'rows');
ib = sort(ib);
nodes1 = nodes1(ib,:);

nodes_total = [nodes;nodes1];

nodes = nodes_total;
nodes(:,2) = nodes(:,2) - min(nodes(:,2));
% Find temperature of each node
tempProfile = xlsread('IRprofile_P_AND_G.xlsx');
tempProfile(:,2) = flipud(tempProfile(:,2));
tempProfile(:,3) = flipud(tempProfile(:,3));

for i = 1:(length(nodes)+1)/2
nodes(i,3) = interp1q(tempProfile(:,1),tempProfile(:,3),nodes(i,2));
    
end

for i = round ((length(nodes)+1)/2):length(nodes)
   nodes(i,3) = interp1q(tempProfile(:,1),tempProfile(:,2),nodes(i,2));
end
figure
plot(nodes(:,1),nodes(:,2),'k.');
hold on
plot(nodes(boundaryNodes(:,1),1),nodes(boundaryNodes(:,1),2),'b-')
axis off
axis equal
hold on
C = [];
for i = 1:length(boundaryNodes)
    
    if ( boundaryNodes(i,2) == 2)
        color = 'r.';
        hold on
        plot(nodes(boundaryNodes(i,1),1),nodes(boundaryNodes(i,1),2),color);
        
    elseif (boundaryNodes(i,2)  == 5)
        color = 'go';
        hold on
        plot(nodes(boundaryNodes(i,1),1),nodes(boundaryNodes(i,1),2),color);
        
    elseif ( boundaryNodes(i,2)  == 0)
        color = 'k.';
        hold on
        plot(nodes(boundaryNodes(i,1),1),nodes(boundaryNodes(i,1),2),color);
        
    elseif ( boundaryNodes(i,2)  == 4)
        color = 'go';
        hold on
        plot(nodes(boundaryNodes(i,1),1),nodes(boundaryNodes(i,1),2),color);
                
    elseif ( boundaryNodes(i,2)  == 6)
        color = 'mo';
        hold on
        plot(nodes(boundaryNodes(i,1),1),nodes(boundaryNodes(i,1),2),color);   
        
        
    else
      
         
    end
    
    
end
% 
% 
% 
% write files
dlmwrite('../../problems/preform/preform.nodes',[length(nodes),1],'delimiter',' ')
dlmwrite('../../problems/preform/preform.nodes',nodes,'-append','delimiter',' ')

%segments
dlmwrite('../../problems/preform/preform.boundary',length(boundaryNodes),'delimiter',' ')
dlmwrite('../../problems/preform/preform.boundary',boundaryNodes,'-append',.........
    'delimiter',' ')

