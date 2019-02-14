% Draw Planar-Straight-Line-Graph of a preform
clc;
clear all;
close all;

%% RADII

% bot diameter
DIN_BOT = 13.97;
DOUT_BOT = 22.50;


% mid diameter
DOUT_MID = 23.01;
DIN_MID = 14.61;


% top diamter
DOUT_TOP = 29.49;
DIN_TOP = 24.31;


% thickness of preform
THICKNESS_BOT = (DOUT_BOT - DIN_BOT )/2;
THICKNESS_MID = (DOUT_MID - DIN_MID )/2;
THICKNESS_TOP = (DOUT_TOP - DIN_TOP)/2;
THICKNESS_BOT_BOT = 3.2;

% total preform length
TOTAL_LENGTH = 97.18;

% number of nodes through the thickness
NUM_NODES_THICKNESS = 3;

R_REDUCED_BOT = (THICKNESS_BOT - THICKNESS_BOT_BOT);


%% LENGTHS
% middle
L_m = TOTAL_LENGTH -37.31-11.25;
L_t1 = 3;
L_t = 37.31-19.48-L_t1;



%% OUTSIDE SURFACE
% bottom
Din_b = DIN_BOT;
Dout_b = DOUT_BOT;

Din_m = DIN_MID;
Dout_m = DOUT_MID;


Rin_bot = Din_b/2 + (3/4)*R_REDUCED_BOT;
Rout_bot = Dout_b/2 - (1/4)*R_REDUCED_BOT; 
Rin_bot_s = Rin_bot;
Rout_bot_s = Rout_bot;


Din_t = DIN_TOP;
Dout_t = DOUT_TOP;


%% NUMBER OF NODES
N1_a =50;
N1_b = 12;
N1 = N1_a + N1_b;
N2 = 22;
ntheta = 15;
N3 = 5;
N4 = 4;
N5 = 5;
N6 = 5;
N7 = 6;
N1 = N1+N7;




% starting from 0,0
nodes = [];
count = 1;

theta_max = acosd((Din_b/2)/(Rin_bot));
theta = linspace(-90,-theta_max,ntheta);




for i = 1:length(theta)
    nodes(count,:) = [Rin_bot*cosd(theta(i)), Rin_bot*sind(theta(i))];
    count = count +1;
    
end

height_bot = max(nodes(:,2));
nodes;





% Left wall
%% traction nodes
nodes(count:1:count+N7-1,:) = [linspace(Din_b/2,Din_b/2,N7)',linspace(height_bot,0,N7)'];
count = count + N7;
nodes(count:1:count+N1_a-1,:) = [linspace(Din_b/2,Din_m/2,N1_a)',linspace(0,(5/6)*L_m,N1_a)'];
count = count + N1_a;
nodes(count:1:count+N1_b-1,:) = [linspace(Din_m/2,Din_m/2,N1_b)',linspace((5/6)*L_m,L_m,N1_b)'];
count = count + N1_b;
nodes(count:1:count+N2-1,:) = [linspace(Din_m/2,Din_t/2,N2)',linspace(L_m,L_t+L_m,N2)'];
count = count + N2;
nodes(count:1:count+N3-1,:) = [linspace(Din_t/2,Din_t/2,N3)',linspace(L_m+L_t,L_t+L_m+L_t1,N3)'];
count = count + N3;

height = nodes(end,2);
%% Essential boundary nodes
nodes(count:1:count+N4-1,:) = [linspace(Din_t/2,Dout_t/2,N4)',linspace(height,height,N4)'];
count = count + N4;
%% Outside nodes
% right wall
nodes(count:1:count+N3-1,:) = [linspace(Dout_t/2,Dout_t/2,N3)',linspace(height,L_t+L_m,N3)'];
count = count + N3;
nodes(count:1:count+N2-1,:) = [linspace(Dout_t/2,Dout_m/2,N2)',linspace(L_t+L_m,L_m,N2)'];
count = count + N2;
% right wall
nodes(count:1:count+N1_b-1,:) = [linspace(Dout_m/2,Dout_m/2,N1_b)',linspace(L_m,(5/6)*L_m,N1_b)'];
count = count + N1_b;
nodes(count:1:count+N1_a-1,:) = [linspace(Dout_m/2,Dout_b/2,N1_a)',linspace((5/6)*L_m,0,N1_a)'];
count = count + N1_a;

% bottom joining portion
theta_min = asind(height_bot/Rout_bot);
R_bot_stop = Rout_bot*cosd(theta_min);
nodes(count:1:count+N6-1,:) = [linspace(Dout_b/2,R_bot_stop,N6)',linspace(0,height_bot,N6)'];
count = count + N6;



theta = linspace(theta_min,-90,ntheta);

for i = 1:length(theta)
    nodes(count,:) = [Rout_bot*cosd(theta(i)), Rout_bot*sind(theta(i))];
    count = count +1;
    
end
% Bottom edge
nodes(count:1:count+N5-1,:) = [linspace(0,0,N5)',linspace(-Rout_bot,-Rin_bot,N5)'];

[~,ib] = unique(nodes,'rows');

ib = sort(ib);
nodes = nodes(ib,:);

%nodes(:,2) = nodes(:,2) - min(nodes(:,2));




% boundary nodes
boundaryNodes = linspace(1,length(nodes),length(nodes))';
boundaryNodes(1:(ntheta + N1+N2 -5 ),2) = 2;
boundaryNodes((ntheta+(N1+N2)-4):((ntheta+(N1+N2+(N3-4)-1+N3))+ N4-2),2) = 5;
boundaryNodes(end:-1:end-(N5-2),2) = 4;

%%  Geoemtric measures ( Change these for different geometries) 
% bottom
nodes1 = [];
count = 1;

BOT_DIAMETER = linspace(DIN_BOT,DOUT_BOT,NUM_NODES_THICKNESS+2);
BOT_BOT_DIAMETER = linspace(Rin_bot*2,Rout_bot*2,NUM_NODES_THICKNESS+2);
MID_DIAMETER = linspace(DIN_MID,DOUT_MID,NUM_NODES_THICKNESS+2);
TOP_DIAMETER = linspace(DIN_TOP,DOUT_TOP,NUM_NODES_THICKNESS+2);
for i = 1:NUM_NODES_THICKNESS
    
    
Din_b = BOT_DIAMETER(i+1);

Din_m = MID_DIAMETER(i+1);
Din_t = TOP_DIAMETER(i+1);
Rin_bot = BOT_BOT_DIAMETER(i+1)/2;

theta_max = asind((height_bot)/(Rin_bot));

theta = linspace(-85,theta_max,ntheta);


for i = 1:length(theta)
    nodes1(count,:) = [Rin_bot*cosd(theta(i)), Rin_bot*sind(theta(i))];
    count = count +1;
    
end
height = Rin_bot*sind(theta(end));
R_stop = Rin_bot*cosd(theta(end));
% Left wall
%% traction nodes
nodes1(count:1:count+N7-1,:) = [linspace(R_stop,Din_b/2,N7)',linspace(height,0,N7)'];
count = count +N7;
nodes1(count:1:count+N1_a-1,:) = [linspace(Din_b/2,Din_m/2,N1_a)',linspace(0,(5/6)*L_m,N1_a)'];
count = count + N1_a;
nodes1(count:1:count+N1_b-1,:) = [linspace(Din_m/2,Din_m/2,N1_b)',linspace((5/6)*L_m,L_m,N1_b)'];
count = count+N1_b;
nodes1(count:1:count+N2-1,:) = [linspace(Din_m/2,Din_t/2,N2)',linspace(L_m,L_t+L_m,N2)'];
count = count + N2;
nodes1(count:1:count+N3-1,:) = [linspace(Din_t/2,Din_t/2,N3)',linspace(L_m+L_t,L_t+L_m+L_t1-0.5,N3)'];
count = count + N3;

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

