function [nodes,segments] = get_node_temperature(OIL_TEMPERATURE,COOLING_TIME)


% if less lose the through thickness temperature profile



% if in neck region we have to modify this problem 
[X,Y,TEMPERATURE] = temperature_profile(OIL_TEMPERATURE,COOLING_TIME);

X;
Y;
Y = Y + 77.44;

Z = TEMPERATURE(:);


%% 

nodes = [];


% r = 0 marker = 4
count = length(nodes);
for i = 1:size(Y,2)
    
    nodes(count+i,1:2) = [X(1,size(Y,2)-(i-1)),Y(1,size(Y,2)-(i-1))];
    nodes(count+i,3) = TEMPERATURE(1,size(Y,2)-(i-1));
    
    
    if ( i < size(Y,2))
    segments(count+i,1:3) = [count+i,count+i+1,4];   
    end
end



% inside traction boundary marker  = 2
count = length(nodes);
count_1 = length(segments);

for i = 2:length(Y)-1
    nodes(count+i-1,1:2) = [X(i,1),Y(i,1)];
    nodes(count+i-1,3) = TEMPERATURE(i,1);

    
    if ( i < length(Y) )
     segments(count_1+(i-1),1:3) = [count_1+(i-1),count_1+(i-1)+1,2];   
    end
     
end
count_1 = length(segments);

segments(count_1+1,1:3) = [count_1,count_1+1,2];
segments(count_1+2,1:3) = [count_1+1,count_1+2,2];


count = length(nodes)
count_1 = length(segments);

% top boundary marker = 5
for i = 1:size(Y,2)
    
    nodes(count+i,1:2) = [X(end,i),Y(end,i)];
    nodes(count+i,3) = TEMPERATURE(end,i);
    
    if ( i < size(Y,2))
    
    segments(count_1+i,1:3) = [count_1+i-1,count_1+i,5];   
    end

end

count = length(nodes);
count_1 = length(segments);

%outside marker = 6
for i = 1:length(Y)-2
    nodes(count+i,1:2) = [X(end-i,end),Y(end-i,end)];
    
     nodes(count+i,3) = TEMPERATURE(end-i,end);

    
    
    if ( i < length(Y) - 2)
     segments(count_1+i,1:3) = [count_1+i-1,count_1+i,6];   
       
    end

end
count_1 = length(segments);

segments(count_1+1,1:3) = [count_1,count_1+1,6];
segments(count_1+2,1:3) = [count_1+1,1,6];

count = length(nodes);


num_boundary = length(nodes);

    


for i = 2:length(Y)-1
    
    for j = 2:size(Y,2) -1
     nodes(count+j-1,1:2) = [X(i,j),Y(i,j)]  ; 
      nodes(count+j-1,3) = TEMPERATURE(i,j); 
    end
    
    count = length(nodes);
     
end

length(nodes)



x = nodes(:,1);
y = nodes(:,2);
length(nodes);

