clear all
close all 
clc;

r0 = 20;
r1 = 21.57;




numCircle = 32;
num_thickness = 1;


% draw boundary first
theta = linspace(0,90,numCircle);
nodes = [];
countNodes = 1;
for i = 1:length(theta)
    nodes(countNodes,1) = r0 * cosd(theta(i));
    nodes(countNodes,2) = r0*sind(theta(i));
    countNodes = countNodes + 1;
end

for i = 1:num_thickness+2
    thickness_coords = linspace(r0,r1,num_thickness+2)';
end


nodes(countNodes+1:countNodes+num_thickness,:) = [linspace(0,0,num_thickness)',thickness_coords(2:end-1)];

countNodes = countNodes + num_thickness;

theta = linspace(90,0,numCircle);

for i = 1:length(theta)
    nodes(countNodes,1) = r1 * cosd(theta(i));
    nodes(countNodes,2) = r1*sind(theta(i));
    countNodes = countNodes + 1;
end

nodes(countNodes+1:countNodes+num_thickness,:) = [linspace(r0,r1,num_thickness)',linspace(0,0,num_thickness)'];




segments = [];
for i = 1:length(nodes)

    if ( i < (numCircle ))   
        if ( i == length(nodes))
            segments(i,:) = [i,1,2];
        else
            segments(i,:) = [i,i+1,2];
        end
        
    elseif ( (i < (numCircle+1)) && ( i > (numCircle -1))  )
        if ( i == length(nodes))
            segments(i,:) = [i,1,3];
        else
            segments(i,:) = [i,i+1,3];
        end
    elseif ( ( i < (numCircle + numCircle )) && ( i > numCircle ))
        if ( i == length(nodes))
            segments(i,:) = [i,1,0];
        else
            segments(i,:) = [i,i+1,0];
        end
    else
        if ( i == length(nodes))
            segments(i,:) = [i,1,8];
        else
            segments(i,:) = [i,i+1,8];
        end
    end
    
end

% Plot it 
figure
plot(nodes(:,1),nodes(:,2),'k.');
axis off
axis equal 
    rand1 = rand(1,3);
hold on 
fill(nodes(:,1),nodes(:,2),rand1);
hold on 
plot(nodes(:,1),nodes(:,2),'k-');




%% write files
dlmwrite('cylinder.nodes',[length(nodes),0],'delimiter',' ')
dlmwrite('cylinder.nodes',nodes,'-append','delimiter',' ')

%segments
dlmwrite('cylinder.segs',length(segments),'delimiter',' ')
dlmwrite('cylinder.segs',segments,'-append',.........
    'delimiter',' ')


