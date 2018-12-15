function [poly,num_contours,num_verticies] = read_polygon(filename)

fid = fopen(filename,'rt');

count = 1;

while ~feof(fid)
    thisline = fgetl(fid);
    if ( count == 1)
        num_contours = str2num(thisline);
    end
    if ( count == 2)
        num_verticies = str2num(thisline);
        i = 1;
    end
    
    if ( count > 2)  
        [token,remain] = strtok(thisline);
        poly(i,1) = str2double(token);
        poly(i,2) = str2double(remain);
        i = i + 1;
    end
    
    count = count +1;
end

fclose(fid);

end

