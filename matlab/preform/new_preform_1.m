clear all 
close all




%% Diameters

% top
internal_diameter_top = 21.4;
external_diameter_top = 24.5;
length_top = 2; 


% taper
external_diameter_taper_top = 24.5;
external_diameter_taper_bot = 24.5;
internal_diameter_taper_top = 21.4;
internal_diameter_taper_bot = 18.63;
length_taper = 7.5;



% 
radius_neck = 20;
radius_top = 10;


%sidewall
length_of_side_wall = 96.7-length_top - length_taper - 12.6;


%bot
internal_diameter_bot = 17.56;
external_diameter_bot = 23.59;
wall_bot_thickness = 2.39;


bot_outside_diameter = 23.59;


% NUM NODES
N_JOIN = 3;
N_THETA = 15;
N_wall_in = 60;

N_taper_in = 5;
N_taper_angle = 6;
N_wall_out = N_wall_in;
N_wall_bot = 3;
N_THETA_OUT = 15;
N_neck_top = 3;
N_neck_in = 3;
N_neck_out = 3;
theta_taper_max = 15;


%% DRAW THE PREFORM OUTSIDE SURFACE

% draw inner surface
r_bot = (bot_outside_diameter/2-wall_bot_thickness);

theta_max = asind(internal_diameter_bot/(2*r_bot));
% draw bottom

theta = linspace(0,theta_max,N_THETA);

for i = 1:length(theta)
   nodes(i,1) = r_bot * sind(theta(i));
   nodes(i,2) = -r_bot * cosd(theta(i));
end


% side wall to zero
x_max = max(nodes(:,1));
y_max = max(nodes(:,2));

y = linspace(y_max,0,N_JOIN+1);

for i = 1:N_JOIN
    nodes = [nodes ; [x_max,y(i+1)]];
    
end

% internal side wall
x_max = max(nodes(:,1));
y_max = max(nodes(:,2));





y_taper_start = length_of_side_wall-20*sind(theta_taper_max /2);
x_taper_max = 18.44/2 -(20-20*cosd(theta_taper_max /2));
x = linspace(x_max,x_taper_max,N_wall_in+1);
y = linspace(y_max,y_taper_start,N_wall_in+1);

for i = 1:N_wall_in
    nodes = [nodes ; [x(i+1),y(i+1)]];
    
end



% Inside taper 
x_max = max(nodes(:,1));
y_max = max(nodes(:,2));

x = linspace(x_max,internal_diameter_taper_top/2,N_taper_in+1);
y = linspace(y_max,length_of_side_wall+length_taper,N_taper_in+1);




theta = linspace(0,theta_taper_max,N_taper_angle);
y_taper = [];
for i = 1:length(theta)
    
    nodes = [nodes ; [(x_max+20) - 20*cosd(theta(i)),y_max + 20*sind(theta(i))]];
    y_taper(i) = y_max + 20*sind(theta(i));
    
end


% Inside taper 
x_max = max(nodes(:,1));
y_max = max(nodes(:,2));

x = linspace(x_max,internal_diameter_taper_top/2,N_taper_in+1);
y = linspace(y_max,length_of_side_wall+length_taper,N_taper_in+1);

count  = length(y_taper);

for i = 1:N_taper_in
    nodes = [nodes ; [x(i+1),y(i+1)]];
    y_taper(count + i) = y(i+1);

end

y_taper = flip(y_taper);



%neck inner wall
x_max = max(nodes(:,1));
y_max = max(nodes(:,2));

x = linspace(x_max,x_max,N_neck_in+1);
y = linspace(y_max,y_max+length_top,N_neck_in+1);

for i = 1:N_neck_in
    nodes = [nodes ; [x(i+1),y(i+1)]];
    
end


% neck top
x_max = max(nodes(:,1));
y_max = max(nodes(:,2));

x = linspace(x_max,external_diameter_taper_top/2,N_neck_top+1);
y = linspace(y_max,y_max,N_neck_top+1);

for i = 1:N_neck_top
    nodes = [nodes ; [x(i+1),y(i+1)]];
    
end



% neck outer wall
x_max = max(nodes(:,1));
y_max = max(nodes(:,2));

x = linspace(x_max,x_max,N_neck_out+1);
y = linspace(y_max,y_max-length_top,N_neck_out+1);

for i = 1:N_neck_out
    nodes = [nodes ; [x(i+1),y(i+1)]];
    
end

%outside taper
x_max = max(nodes(:,1));
y_max = max(nodes(:,2));

x = linspace(x_max,external_diameter_taper_bot/2,length(y_taper));

for i = 1:length(y_taper)
    nodes = [nodes ; [x(i),y_taper(i)]];
    
end
% outside wall
x_max = max(nodes(:,1));
y_max = nodes(end,2);

x = linspace(x_max,external_diameter_bot/2,N_wall_out+1);
y = linspace(y_max,0,N_wall_out+1);

for i = 1:N_wall_out
    nodes = [nodes ; [x(i+1),y(i+1)]];
    
end
% outside cap

theta = linspace(90,0,N_THETA_OUT);

count = length(nodes);

for i = 1:length(theta)
   nodes(count + i,1) =  external_diameter_bot/2 * sind(theta(i));
   nodes(count +i,2) = -external_diameter_bot/2 *cosd(theta(i));
end


%outside wall

x = linspace(0,0,N_wall_bot+1);
y = linspace(-external_diameter_bot/2,-r_bot,N_wall_bot+1);

for i = 1:N_wall_bot
    nodes = [nodes ; [x(i+1),y(i+1)]];
    
end
[ix,iy] = uniquetol(nodes,'ByRows',true)
nodes = nodes(sort(iy),:);
boundary_nodes = nodes;





%% THICKNESS NODES

NUM_THICKNESS = 2;

BOT_WALL_DIAMETER = linspace(internal_diameter_bot,external_diameter_bot,NUM_THICKNESS+2);

BOT_TAPER_DIAMETER = linspace(internal_diameter_bot,external_diameter_taper_bot,NUM_THICKNESS+2);
TOP_TAPER_DIAMETER = linspace(internal_diameter_taper_top,external_diameter_taper_bot,NUM_THICKNESS+2);
R_BOT = linspace(r_bot,external_diameter_bot/2,NUM_THICKNESS+2);


for k = 1:NUM_THICKNESS

count = length(nodes);
% draw inner surface
r_bot = R_BOT(k+1);
N_THETA = 15;

theta = linspace(0,theta_max,N_THETA+1);

for i = 1:length(theta)-1
   nodes(count+i,1) = r_bot * sind(theta(i+1));
   nodes(count+i,2) = -r_bot * cosd(theta(i+1));
end


% side wall to zero
x_max = nodes(end,1);
y_max = nodes(end,2);

y = linspace(y_max,0,N_JOIN+1);
x_max = BOT_WALL_DIAMETER(k+1)/2;
for i = 1:N_JOIN
    nodes = [nodes ; [x_max,y(i+1)]];
    
end

% internal side wall
x_max = nodes(end,1);
y_max = nodes(end,2);

x = linspace(x_max,BOT_TAPER_DIAMETER(k+1)/2,N_wall_in+1);
y = linspace(y_max,y_taper_start,N_wall_in+1);

for i = 1:N_wall_in
    nodes = [nodes ; [x(i+1),y(i+1)]];
    
end


% Inside taper 
x_max = nodes(end,1);
y_max = nodes(end,2);

y_taper = linspace(length_of_side_wall+length_taper,length_of_side_wall,N_taper_in+1);

theta = linspace(0,theta_taper_max,N_taper_angle);

for i = 1:length(theta)
    
    nodes = [nodes ; [(x_max+20) - 20*cosd(theta(i)),y_max + 20*sind(theta(i))]];
    
end

% Inside taper 
x_max = nodes(end,1);
y_max = nodes(end,2);

x = linspace(x_max,TOP_TAPER_DIAMETER(k+1)/2,N_taper_in+1);
y = linspace(y_max,length_of_side_wall+length_taper,N_taper_in+1);
for i = 1:N_taper_in
    nodes = [nodes ; [x(i+1),y(i+1)]];
    
end



%neck inner wall
x_max = nodes(end,1);
y_max = nodes(end,2);


x = linspace(x_max,x_max,N_neck_in+1);
y = linspace(y_max,y_max+length_top,N_neck_in+1);

for i = 1:N_neck_in-1
    nodes = [nodes ; [x(i+1),y(i+1)]];
    
end




end

[ix,iy] = uniquetol(nodes,'ByRows',true)
nodes = nodes(sort(iy),:);

%% boundary segments

% SEGEMNTS% Write the segments 
segments = [];
for i = 1:length(boundary_nodes)

        
    if ( i == length(boundary_nodes))
        segments(i,:) = [i,1];
    else
        segments(i,:) = [i,i+1]; 
    end
  
    
end
N_taper_in = N_taper_in + N_taper_angle - 1 ;
% traction
segments(1:N_THETA+N_taper_in+N_wall_in+N_JOIN+2,3) = 2;
count = N_THETA+N_taper_in+N_wall_in+N_JOIN+2;
% h = h_max, encastre
segments(count+1:count+N_neck_top,3) = 5;
count = count + N_neck_top;
% contact
segments(count+1:count+N_taper_in + N_wall_out+N_THETA_OUT-2,3) = 6;
% r = 0 , symmetry
segments(end:-1:end-(N_wall_bot)+1,3) = 4;



%% TEMPERATURE
nodes(:,3) = 101.9;

nodes(:,2) = nodes(:,2) - max(nodes(:,2));

for i = 1:length(nodes)
    
   if ( nodes(i,2) < -85)
       nodes(i,3) = 101.9-1.2*abs(nodes(i,2)+85);
   end

%     if ( nodes(i,2) < -90)
%        nodes(i,3) = 93;
%    end
   
       
%     if ( nodes(i,2) > -6)
%        nodes(i,3) = 98;
%    end
   
   
end




%% PLOT


figure

hold on
plot(nodes(:,1),nodes(:,2),'b.');




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
        color = 'k';
        hold on
        plot(nodes(segments(i,1:2),1),nodes(segments(i,1:2),2),color);
        
    elseif ( segments(i,3) == 4)
        color = 'g';
        hold on
        plot(nodes(segments(i,1:2),1),nodes(segments(i,1:2),2),color);
        
    elseif ( segments(i,3) == 7)
        color = 'y';
        hold on
        plot(nodes(segments(i,1:2),1),nodes(segments(i,1:2),2),color);
    else
        hold on
        plot(nodes(segments(i,1:2),1),nodes(segments(i,1:2),2),'m');
        
        
        
    end

end

axis equal


%% write files
dlmwrite('../../problems/preform/preform.nodes',[length(nodes),1],'delimiter',' ')
dlmwrite('../../problems/preform/preform.nodes',nodes,'-append','delimiter',' ')


%segments
dlmwrite('../../problems/preform/preform.segs',length(segments),'delimiter',' ')
dlmwrite('../../problems/preform/preform.segs',segments,'-append',.........
    'delimiter',' ')

% 

for i = 1:length(segments)
    
   x1 = nodes(segments(i,1),1:2);
   x2 = nodes(segments(i,2),1:2);
   
   a(i) = norm(x1-x2,2)
end

% figure
% plot(boundary_nodes(:,1),boundary_nodes(:,2),'g-')
% axis off
% axis equal

