function [] = convert_tri_file_to_vtk(tri_file,vtk_file, plot_flag)



% write a vtk file 

nodesfile = strcat(tri_file,'.node');
fileID = fopen(nodesfile,'r');

line= fgetl(fileID);
line = split(line,' ');
numnodes = str2num(line{1});

for i = 1:numnodes
    line = fgetl(fileID);
    line = split(line,' ');

    
    
    count = 1;
    
        
    for k = 1:length(line)
        
       if ( ~isempty(line{k}))
           
           if ( count > 1)
               
              nodes(i,count-1) = str2num(line{k}); 
           end
           
           
           
           count = count+1;
           
       end
       % else move on to next line
    end
    
      
end

fclose(fileID);

elefile = strcat(tri_file,'.ele');
fileID = fopen(elefile,'r');


tri = [];
line= fgetl(fileID);
line = split(line,' ');
numele = str2num(line{1});

tri = zeros(numele,3);

for i = 1:numele
    line = fgetl(fileID);
    line = split(line,' ');

    count = 1;
    
        
    for k = 1:length(line)
        
       if ( ~isempty(line{k}))
           
           if ( count > 1)
               
              tri(i,count-1) = str2num(line{k}); 
           end
           
           
           
           count = count+1;
           
       end
       % else move on to next line
    end
    
end

fclose(fileID);

if ( plot_flag)
figure
plot(nodes(:,1),nodes(:,2),'r.');

triplot(tri,nodes(:,1),nodes(:,2));

end
z = zeros(length(nodes),1);

vtkwrite(vtk_file,'polydata','triangle',nodes(:,1),nodes(:,2),z,tri);    









