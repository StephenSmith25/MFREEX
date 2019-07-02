clear all
close all 
clc;

r0 = 20.00;
r1 = 21.57;




numCircle = 20;
num_thickness = 3;


% draw boundary first
theta = linspace(0,90,numCircle);
nodes = [];
countNodes = 1;
for i = 1:length(theta)
    nodes(countNodes,1) = r0 * cosd(theta(i));
    nodes(countNodes,2) = r0*sind(theta(i));
    countNodes = countNodes + 1;
end

y = linspace(r0,r1,4);

for i = 1:2
   nodes(countNodes,1) = 0;
   nodes(countNodes,2) = y(i+1);
   countNodes = countNodes + 1;

end


theta = linspace(90,0,numCircle);

for i = 1:length(theta)
    nodes(countNodes,1) = r1 * cosd(theta(i));
    nodes(countNodes,2) = r1*sind(theta(i));
    countNodes = countNodes + 1;
end


x = linspace(r1,r0,4);

for i = 1:2
   nodes(countNodes,1) = x(i+1);
   nodes(countNodes,2) = 0;
   countNodes = countNodes + 1;

end


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
            segments(i,:) = [i,1,0];
        else
            segments(i,:) = [i,i+1,0];
        end
    elseif ( ( i < (numCircle + numCircle )) && ( i > numCircle ))
        if ( i == length(nodes))
            segments(i,:) = [i,1,0];
        else
            segments(i,:) = [i,i+1,0];
        end
    else
        if ( i == length(nodes))
            segments(i,:) = [i,1,0];
        else
            segments(i,:) = [i,i+1,0];
        end
    end
    
end




% write the polyfile

% 

% Plot it 
figure
plot(nodes(:,1),nodes(:,2),'k.');
axis off
axis equal 
    rand1 = rand(1,3);
% hold on 
% fill(nodes(:,1),nodes(:,2),rand1);
% hold on 
% plot(nodes(:,1),nodes(:,2),'k-');


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
write_triangle_file(nodes,segments,'cylinder.poly');


%% write files
dlmwrite('cylinder.nodes',[length(nodes),0],'delimiter',' ')
dlmwrite('cylinder.nodes',nodes,'-append','delimiter',' ')

%segments
dlmwrite('cylinder.segs',length(segments),'delimiter',' ')
dlmwrite('cylinder.segs',segments,'-append',.........
    'delimiter',' ')




function [] = write_triangle_file(X,segments,filename)


fileID = fopen(filename,'w');

fprintf(fileID,'%d %d %d \n',length(X),2,0);

for i = 1:length(X)
   fprintf(fileID,'%d %f %f\n',i,X(i,1),X(i,2));
    
end
fprintf(fileID,'%d \n',length(segments));

for i = 1:length(segments)
   fprintf(fileID,'%d %d %d\n',i,segments(i,1),segments(i,2));
    
end
fprintf(fileID,'%d \n',0);

fclose(fileID);



end

