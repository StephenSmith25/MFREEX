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
thickness_bot = 3.2;
vDim = Dout_b/2 -Rin_bot -thickness_bot;
% middle
Din_b = (13.97+14.43)/2;
Dout_b = (23.19+22.50)/2;
L_m = 97.18-37.31-11.25;
% top
Din_t = 24.31;
Dout_t = 29.49;
L_t1 = 3;
L_t = 37.31-19.48-L_t1;


N1 = 35;
N2 = 15;



% Draw it
ntheta = 10;
% starting from 0,0
nodes = [];
count = 1;


theta = linspace(-90,-5,ntheta);

for i = 1:length(theta)
    nodes(count,:) = [Rin_bot*cosd(theta(i)), -vDim + Rin_bot*sind(theta(i))];
    count = count +1;
    
end
nodes;
N3 = 4;
N4 = 3;
N5 = 3;
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

theta = linspace(-10,-90,ntheta);

for i = 1:length(theta)
    nodes(count,:) = [Rout_bot*cosd(theta(i)), Rout_bot*sind(theta(i))];
    count = count +1;
    
end
% Bottom edge
nodes(count:1:count+N5-1,:) = [linspace(0,0,N5)',linspace(-Rout_bot,-Rout_bot+thickness_bot,N5)'];

nodes
[~,ib] = unique(nodes,'rows');

ib = sort(ib);
nodes = nodes(ib,:);

nodes(:,2) = nodes(:,2) - min(nodes(:,2));
%% Find temperature of each node
tempProfile = xlsread('IRprofile.xlsx');
tempProfile(:,2) = flipud(tempProfile(:,2));
tempProfile(:,3) = flipud(tempProfile(:,3));

for i = 1:(length(nodes)+1)/2
nodes(i,3) = interp1q(tempProfile(:,1),tempProfile(:,3),nodes(i,2));
    
end

for i = (length(nodes))/2+1:length(nodes)
   nodes(i,3) = interp1q(tempProfile(:,1),tempProfile(:,2),nodes(i,2));
end

% overwrite temperature profile
%%nodes(:,3) = 97.45;


segments = [];
% Write the segments 
for i = 1:length(nodes)

        
    if ( i == length(nodes))
        segments(i,:) = [i,1];
    else
        segments(i,:) = [i,i+1]; 
    end
  
    
end

segments(1:(ntheta + N1+N2+(N3-1)),3) = 2;
segments((ntheta+(N1+N2+(N3-1))):((ntheta+(N1+N2+(N3-1)))+ N4-2),3) = 5;
segments(end:-1:end-(N5-2),3) = 4;


% Plot it 
% figure
% plot(nodes(:,1),nodes(:,2),'k-');
% axis equal 
%     rand1 = rand(1,3);
% hold on 
% fill(nodes(:,1),nodes(:,2),rand1);
% hold on 
% plot(nodes(:,1),nodes(:,2),'k.','linewidth',3);


%% plot segments
figure
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
        color = 'k';
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
axis equal
hold on
plot(nodes(:,1),nodes(:,2),'ko','markersize',3);

% Save files 
% segments and nodes 

% first line nubmer of nodes / number of segments 



%% write files
dlmwrite('../../problems/preform/preform.nodes',[length(nodes),1],'delimiter',' ')
dlmwrite('../../problems/preform/preform.nodes',nodes,'-append','delimiter',' ')

%segments
dlmwrite('../../problems/preform/preform.segs',length(segments),'delimiter',' ')
dlmwrite('../../problems/preform/preform.segs',segments,'-append',.........
    'delimiter',' ')


area = 2*3.2^2*pi + 2*pi*7.1*(60.04-10.3) + pi*(7.1+3.2)*sqrt((7.1-3.2)^2 + (10.3-3.2)^2) + pi*(12.15+7.1)*sqrt((12.15-7.1)^2 + (74.87-60.04)^2) + 2*pi*12.15*(77.87-74.87)
volume_in = (2.00/3)*7.1^2*pi + pi*7.1^2*(60.043-10.3) + (1.00/3)*pi*(7.1^2 + 7.1*12.155+12.155^2)*(74.873-60.043) + pi*12.155^2*(77.87-74.87)



%%volume_in = (2.00/3)*11.42^2*pi + (1.00/3)*pi*(3.2^2 + 3.2*7.1+7.1^2)*(10.3-3.2) + pi*7.1^2*(60.043-10.3) + (1.00/3)*pi*(7.1^2 + 7.1*12.15+12.15^2)*(74.87-60.04) + pi*12.15^2*(77.87-74.87)

%%volume_out = volume_in - volume_out;

