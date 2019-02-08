% Draw Planar-Straight-Line-Graph of a preform
clc;
clear all;
close all;


%%  Geoemtric measures ( Change these for different geometries) 
% bottom
Din_b = (13.97+14.43)/2;
Dout_b = (23.19+22.50)/2;


Rin_bot = Din_b/2;
Rout_bot = Dout_b/2;


% middle
L_m = 97.18-37.31-11.25;
% top
Din_t = 24.31;
Dout_t = 29.49;
L_t1 = 3;
L_t = 37.31-19.48-L_t1;


N1 = 45;
N2 = 25;



% Draw it
ntheta = 15;
% starting from 0,0
nodes = [];
count = 1;


theta = linspace(-90,0,ntheta);

for i = 1:length(theta)
    nodes(count,:) = [Rin_bot*cosd(theta(i)), Rin_bot*sind(theta(i))];
    count = count +1;
    
end
nodes;
N3 = 5;
N4 = 4;
N5 = 5;
% Left wall
%% traction nodes
nodes(count:1:count+N1-1,:) = [linspace(Din_b/2,Din_b/2,N1)',linspace(0,L_m-1,N1)'];
count = count + N1 ;
nodes(count:1:count+N2-1,:) = [linspace(Din_b/2,Din_t/2,N2)',linspace(L_m,L_t+L_m,N2)'];
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
nodes(count:1:count+N2-1,:) = [linspace(Dout_t/2,Dout_b/2,N2)',linspace(L_t+L_m,L_m,N2)'];
count = count + N2;
% right wall
nodes(count:1:count+N1-1,:) = [linspace(Dout_b/2,Dout_b/2,N1)',linspace(L_m-1,0,N1)'];





count = count + N1;

theta = linspace(0,-90,ntheta);

for i = 1:length(theta)
    nodes(count,:) = [Rout_bot*cosd(theta(i)), Rout_bot*sind(theta(i))];
    count = count +1;
    
end
% Bottom edge
nodes(count:1:count+N5-1,:) = [linspace(0,0,N5)',linspace(-Rout_bot,-Rin_bot,N5)'];

nodes
[~,ib] = unique(nodes,'rows');

ib = sort(ib);
nodes = nodes(ib,:);

%nodes(:,2) = nodes(:,2) - min(nodes(:,2));




% boundary nodes
boundaryNodes = linspace(1,length(nodes),length(nodes))';
boundaryNodes(1:(ntheta + N1+N2)-1,2) = 2;
boundaryNodes((ntheta+(N1+N2)):((ntheta+(N1+N2+(N3-1)-1+N3))+ N4-2),2) = 5;
boundaryNodes(end:-1:end-(N5-2),2) = 4;

%%  Geoemtric measures ( Change these for different geometries) 
% bottom

Din_b = (13.97+14.43)/2 +1.5;
Dout_b = (23.19+22.50)/2  -1.5;



Rin_bot = (Din_b/2) ;
Rout_bot = (Dout_b/2);


% middle
L_m = 97.18-37.31-11.25;
% top
Din_t = 24.31   +1;
Dout_t = 29.49  -1  ;
L_t1 = 3;
L_t = 37.31-19.48-L_t1 ;


N1 = 45;
N2 = 25;



% Draw it
ntheta = 15;
% starting from 0,0
nodes1 = [];
count = 1;


theta = linspace(-85,0,ntheta);

for i = 1:length(theta)
    nodes1(count,:) = [Rin_bot*cosd(theta(i)),  Rin_bot*sind(theta(i))];
    count = count +1;
    
end
N3 = 5;
N4 = 4;
N5 = 3;
% Left wall
%% traction nodes
nodes1(count:1:count+N1-1,:) = [linspace(Din_b/2,Din_b/2,N1)',linspace(0,L_m-1,N1)'];
count = count + N1 ;
nodes1(count:1:count+N2-1,:) = [linspace(Din_b/2,Din_t/2,N2)',linspace(L_m,L_t+L_m,N2)'];
count = count + N2;
nodes1(count:1:count+N3-1,:) = [linspace(Din_t/2,Din_t/2,N3)',linspace(L_m+L_t,L_t+L_m+L_t1-0.5,N3)'];
count = count + N3;
% 
height = nodes1(end,2);
%% Outside nodes
% right wall
nodes1(count:1:count+N3-1,:) = [linspace(Dout_t/2,Dout_t/2,N3)',linspace(height,L_t+L_m,N3)'];
count = count + N3;
nodes1(count:1:count+N2-1,:) = [linspace(Dout_t/2,Dout_b/2,N2)',linspace(L_t+L_m,L_m,N2)'];
count = count + N2;
% right wall
nodes1(count:1:count+N1-1,:) = [linspace(Dout_b/2,Dout_b/2,N1)',linspace(L_m-1,0,N1)'];
% 
% 
% 
%
% 
count = count + N1;

theta = linspace(0,-85,ntheta);

for i = 1:length(theta)
    nodes1(count,:) = [Rout_bot*cosd(theta(i)), Rout_bot*sind(theta(i))];
    count = count +1;
end

Din_b = (13.97+14.43)/2 +3.5;
Dout_b = (23.19+22.50)/2  -3.5;



Rin_bot = (Din_b/2) ;
Rout_bot = (Dout_b/2);


% middle
L_m = 97.18-37.31-11.25;
% top
Din_t = 24.31   +2;
Dout_t = 29.49  -2  ;
L_t1 = 3;
L_t = 37.31-19.48-L_t1 ;


N1 = 45;
N2 = 25;



% Draw it
ntheta = 15;
% starting from 0,0


theta = linspace(-85,0,ntheta);

for i = 1:length(theta)
    nodes1(count,:) = [Rin_bot*cosd(theta(i)),  Rin_bot*sind(theta(i))];
    count = count +1;
    
end
N3 = 5;
N4 = 4;
N5 = 3;
% Left wall
%% traction nodes
nodes1(count:1:count+N1-1,:) = [linspace(Din_b/2,Din_b/2,N1)',linspace(0,L_m-1,N1)'];
count = count + N1 ;
nodes1(count:1:count+N2-1,:) = [linspace(Din_b/2,Din_t/2,N2)',linspace(L_m,L_t+L_m,N2)'];
count = count + N2;
nodes1(count:1:count+N3-1,:) = [linspace(Din_t/2,Din_t/2,N3)',linspace(L_m+L_t,L_t+L_m+L_t1-0.5,N3)'];
count = count + N3;
% 
height = nodes1(end,2);
%% Outside nodes
% right wall
nodes1(count:1:count+N3-1,:) = [linspace(Dout_t/2,Dout_t/2,N3)',linspace(height,L_t+L_m,N3)'];
count = count + N3;
nodes1(count:1:count+N2-1,:) = [linspace(Dout_t/2,Dout_b/2,N2)',linspace(L_t+L_m,L_m,N2)'];
count = count + N2;
% right wall
nodes1(count:1:count+N1-1,:) = [linspace(Dout_b/2,Dout_b/2,N1)',linspace(L_m-1,0,N1)'];
% 
% 
% 
%
% 
count = count + N1;

theta = linspace(0,-85,ntheta);

for i = 1:length(theta)
    nodes1(count,:) = [Rout_bot*cosd(theta(i)), Rout_bot*sind(theta(i))];
    count = count +1;
end
[~,ib] = unique(nodes1,'rows');

ib = sort(ib);
nodes1 = nodes1(ib,:);
% 
% end
% % Bottom edge
% nodes1(count:1:count+N5-1,:) = [linspace(0,0,N5)',linspace(-Rout_bot,-Rout_bot,N5)'];
% 

% nodes1(:,2) = nodes1(:,2) - min(nodes1(:,2));
% %% Find temperature of each node
% tempProfile = xlsread('IRprofile.xlsx');
% tempProfile(:,2) = flipud(tempProfile(:,2));
% tempProfile(:,3) = flipud(tempProfile(:,3));
% 
% for i = 1:(length(nodes1)+1)/2
% nodes1(i,3) = interp1q(tempProfile(:,1),tempProfile(:,3),nodes(i,2));
%     
% end
% 
% for i = round ((length(nodes1)+1)/2):length(nodes1)
%    nodes1(i,3) = interp1q(tempProfile(:,1),tempProfile(:,2),nodes(i,2));
% end
% 
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

