clear all
close all 

fileID = fopen("./../../build/bin/tube/cells.txt");


tline = fgetl(fileID);
C = strsplit(tline);
num_verticies = str2double(C{1});
dim = str2double(C{2});

verticies = zeros(num_verticies,dim);

for i = 1:num_verticies


    tline = fgetl(fileID);
    C = strsplit(tline);
    
    for k = 1:dim
        verticies(i,k) = str2double(C{k});
    end
    
        
    
end


figure
plot(verticies(:,1),verticies(:,2),'b*')

axis equal

tline = fgetl(fileID);
C = strsplit(tline);
num_cells = str2double(C{1});

figure

for i = 1:num_cells

    tline = fgetl(fileID);
    C = strsplit(tline);
    
    num_cell_verticies = str2double(C{1});
    
    poly = zeros(num_cell_verticies,2);
    for k = 1:num_cell_verticies
        poly(k,1) = verticies(str2double(C{1+k})+1,1);
        poly(k,2) = verticies(str2double(C{1+k})+1,2);

    end
    area(i) = polyarea(poly(:,1),poly(:,2));
    
    hold on 
    fill(poly(:,1),poly(:,2),rand(1,3));
end
axis equal

% while ischar(tline)
%     disp(tline)
%     [token,remain] = strtok(tline)
%     while ~isempty(remain)
%         
%         [token,remain] = strtok(remain)
%         
%     end
%     tline = fgetl(fileID);
% end