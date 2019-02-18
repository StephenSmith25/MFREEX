% Draw Planar-Straight-Line-Graph of a preform
clc;
clear all;
close all;

%% GOEMETRY

Radius_in = 13.97/2;
Radius_out = 22.50/2;

length_preform = 70;
length_wall = length_preform - Radius_in;

NUM_NODES_THICKNESS=4;
%% NUMBER OF NODES
N1 = 80;
N2 = 5;
N3 = 5;
ntheta = 20;

% starting from 0,0
nodes = [];
count = 1;
theta = linspace(-90,0,ntheta);

for i = 1:length(theta)
    nodes(count,:) = [Radius_in*cosd(theta(i)), Radius_in*sind(theta(i))];
    count = count +1;
    
end
nodes;

% Left wall
%% traction nodes
nodes(count:1:count+N1-1,:) = [linspace(Radius_in,Radius_in,N1)',linspace(0,length_wall,N1)'];
count = count + N1;
height = nodes(end,2);
%% Essential boundary nodes
nodes(count:1:count+N2-1,:) = [linspace(Radius_in,Radius_out,N2)',linspace(height,height,N2)'];
count = count + N2;
%% Outside nodes
% right wall
nodes(count:1:count+N1-1,:) = [linspace(Radius_out,Radius_out,N1)',linspace(height,0,N1)'];
count = count + N1;

theta = linspace(0,-90,ntheta);

for i = 1:length(theta)
    nodes(count,:) = [Radius_out*cosd(theta(i)), Radius_out*sind(theta(i))];
    count = count +1;
    
end
% Bottom edge
nodes(count:1:count+N3-1,:) = [linspace(0,0,N3)',linspace(-Radius_out,-Radius_in,N3)'];

nodes
[~,ib] = unique(nodes,'rows');

ib = sort(ib);
nodes = nodes(ib,:);

%nodes(:,2) = nodes(:,2) - min(nodes(:,2));




% boundary nodes
boundaryNodes = linspace(1,length(nodes),length(nodes))';
boundaryNodes(1:(ntheta + N1  ),2) = 2;
boundaryNodes((ntheta+(N1)-6):((ntheta+(N1+ N2-2)+5)),2) = 5;
boundaryNodes(((ntheta+(N1+ N2-2)+6)):end-(N3-2),2) = 6;

boundaryNodes(end:-1:end-(N3-2),2) = 4;

%%  Geoemtric measures ( Change these for different geometries) 
% bottom
nodes1 = [];
count = 1;

RADIUS = linspace(Radius_in,Radius_out,NUM_NODES_THICKNESS+2);
for i = 1:NUM_NODES_THICKNESS
    
    

Rin_bot = RADIUS(i+1) ;


theta = linspace(-85,0,ntheta);

for i = 1:length(theta)
    nodes1(count,:) = [Rin_bot*cosd(theta(i)),  Rin_bot*sind(theta(i))];
    count = count +1;
    
end

% Left wall
length_wall_i = linspace(0,length_wall,N1);
length_wall_i  = length_wall_i(1:end-1);
%% traction nodes
nodes1(count:1:count+N1-2,:) = [linspace(Rin_bot,Rin_bot,N1-1)',length_wall_i'];
count = count + N1-1;

end


[~,ib] = unique(nodes1,'rows');
ib = sort(ib);
nodes1 = nodes1(ib,:);

nodes_total = [nodes;nodes1];

nodes = nodes_total;
nodes(:,2) = nodes(:,2) - min(nodes(:,2));
axis equal
%% Find temperature of each node
tempProfile = xlsread('IRprofile.xlsx');
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


%% write files
dlmwrite('../../problems/preform/preform.nodes',[length(nodes),1],'delimiter',' ')
dlmwrite('../../problems/preform/preform.nodes',nodes,'-append','delimiter',' ')

%segments
dlmwrite('../../problems/preform/preform.boundary',length(boundaryNodes),'delimiter',' ')
dlmwrite('../../problems/preform/preform.boundary',boundaryNodes,'-append',.........
    'delimiter',' ')



% draw mould
Radius_mould = 11.25;
Height_mould = 30;
radius_bot_mould = 40;
thickness_mould = 5;
mould_length = 160;
height = 74.2650;


mould_nodes = [];
mould_nodes = [mould_nodes ; [Radius_out,height+5]   ];
mould_nodes = [mould_nodes ; [Radius_out,68.92]   ];
mould_nodes = [mould_nodes ; [Radius_out+Radius_mould,68.92-Radius_mould]   ];
mould_nodes = [mould_nodes ; [Radius_out+Radius_mould,height-mould_length]   ];
mould_nodes = [mould_nodes ; [0,height-mould_length]   ];


mould_nodes = [mould_nodes ; [0,height-mould_length-thickness_mould]   ];
mould_nodes = [mould_nodes ; [Radius_out+Radius_mould+thickness_mould,height-mould_length-thickness_mould]   ];
mould_nodes = [mould_nodes ; [Radius_out+Radius_mould+thickness_mould,68.92-Radius_mould]   ];
mould_nodes = [mould_nodes ; [Radius_out+thickness_mould,68.92]   ];
mould_nodes = [mould_nodes ; [Radius_out+thickness_mould,height+5]   ];
mould_nodes = [mould_nodes ; [Radius_out,height+5]   ];


hold on 
plot(mould_nodes(:,1),mould_nodes(:,2),'b.');

C = rand(3,1);
hold on 
fill(mould_nodes(:,1),mould_nodes(:,2),'b');