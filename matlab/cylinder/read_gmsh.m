close all
clear all

% SET UP PROBLEM DIM 
DIM = 2;

% OPEN THE FILE
fid = fopen('cylinder.msh');

% GET FIRST LINE
tline = fgetl(fid);

% LOOP OVER EACH LINE
while ischar(tline)
    
    
    % CHECK FOR PHYSICAL GROUPS
    if ( strcmp(tline,"$PhysicalNames") == 1)  
         [physcial_groups,num_physcial_grups] = read_physical_names(fid,DIM);
    end
    
    % CHECK IF NODES
    if ( strcmp(tline,"$Nodes") == 1)  
          [nodes,numnodes] = read_nodes(fid,DIM);
    end
    
    % CHECK IF ELEMENTS
    if ( strcmp(tline,"$Elements") == 1)  
          [elements,num_elements] = read_elements(fid,DIM);
    end
    
 
    % MOVE TO NEW LINE
    tline = fgetl(fid);
end
fclose(fid);


figure
plot(nodes(:,1),nodes(:,2),'r.');


function [physical_groups,num_physical_groups] = read_physical_names(fid,DIM)
      
   tline = fgetl(fid);
   num_physical_groups = str2num(tline);
   
   physical_groups = [];
    
   count = 1;

    while (strcmp(tline,'$EndPhysicalNames') ~= 1)
           newStr = split(tline);
           
           
           % CASE PRESSURE
            

           
           tline = fgetl(fid);
    end




end

function [nodes,num_nodes] = read_nodes(fid,DIM)

           % get next line to read number of nodes
           tline = fgetl(fid);
           num_nodes = str2num(tline);
           
 
          
           nodes = zeros(num_nodes,3);
           count = 1;
           tline = fgetl(fid);
           % loop over each line to form nodes

       while (strcmp(tline,'$EndNodes') ~= 1)
           newStr = split(tline);
           % read into nodes array 
           for i = 1:DIM
               nodes(count,i) = str2num(newStr{i+1});
           end
           tline = fgetl(fid);
            count = count+1;

       end

end

function [elements,num_elements] = read_elements(fid,DIM)

     % get next line to read number of elements
      % elementType 
               % 1 -  2D line
               % 2 -  3-Node triangle
               % 3 -  4-node quadrangle 
               % 4 -  4 node tet
               % 5 -  8 node hex
               % 6  - 6 node prism
               
       BLOCKSETS = {};
       % EACH ENTRY CONTAINS
            % entityDim(int) elementType(int; see below) numElementsInBlock(size_t) %elementTag(size_t) nodeTag(size_t) ...
           % get which type of element it is
           
      %      
      % IF DIM OF ENTITY == DIM OF PROBLEM THEN THESE ARE DOMAIN_ELEMENTS
      % IF DIM OF ENTITY = DIM -1 THEN BOUNDARY CONDITION ELEMENTS, READ
      % FROM PHYSICAL GROUPS
           
           
      tline = fgetl(fid);
      num_elements = str2num(tline);
      elements = [];
      
      % BLOCK SET
       % CONTAINS DOMAIN SET - TRIANGLES WITH MAT PROPERTIETIES
       % TRIANGLES
       % NODES
       % MAT_TYPE
       % MAT CONSTANTS 
      BLOCK_SET = [];
      
      % NODESET
      % - CONTAINS DISPLACEMENT AND VELOCITY CONSTRAINTS
      NODE_SET_1 = [];
      NODESETS = {1,NODE_SET_1};

      % SIDESET
       % - CONTAINS PRESSURE
       % - PMAG 
       % - RAMP?
      SIDESET_1 = [];
    
                 
       while (strcmp(tline,'$EndElements') ~= 1)
           tline = fgetl(fid);
       end      
           
           
           

end

