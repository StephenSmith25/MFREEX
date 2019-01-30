function [vorCells,nodes] = readCells(fileName)


nodes = csvread(strcat(fileName,".sites"));
fileID = fopen(strcat(fileName,".edges"));
tline = fgetl(fileID);
countEdges = 1;
vorCells = cell(length(nodes),1);
countSite = 1;
cellEdges = [];

while (tline ~= -1)
    
    tline = fgetl(fileID);
    % seperate line into tokens
    
    % first reading x coord
    coord = 1;
    [token,remain] = strtok(tline,' ');
    while ( token ~=-1)
        if ( strcmp(token,'\') == 1)
            countEdges = 0;
            vorCells{countSite} = cellEdges;
            countSite = countSite + 1;
            cellEdges = [];
        else
            cellEdges(countEdges,coord) = str2num(token);
            % now set to read y coordinate
            coord = 2;
        end
        [token,remain] = strtok(remain,' ');
    end
    countEdges = countEdges +1;   
end
end

