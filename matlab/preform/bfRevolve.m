function [tri, pts] = bfRevolve(bf, bp, nSpokes)
% Revolves the output of boundaryFacets() by a number of times and combines
% the points into faces that can be displayed with trisurf()
%
% bf is the list of boundary facet pairs from boundaryFacets()
% bp is the list of boundary points from boundaryFacets()
% nSpokes is the number of rovolved copies of the profile you want
%
% bp/bf must be a 2D profile
% nSpokes is the TOTAL number of profiles (including the original)

% turn the 2D profile points into 3D aligned on the x/z plane (y=0)
pts = [bp(:,1) zeros(size(bp, 1), 1) bp(:,2)];


% the angles of each spoke
w = 2*pi.*(1:nSpokes-1)./nSpokes;

% revolve each point about z-axis (x=y=0)
for p=1:size(pts, 1)
    x = pts(p,1);
    y = pts(p,2);
    z = pts(p,3);
    
    x2 = x.*cos(w);
    y2 = x.*sin(w);
    z2 = z.*ones(nSpokes-1, 1);
    
    if p>1
        pts2 = [pts2; [x2' y2' z2];];
    else
        pts2 = [x2' y2' z2];
    end
end

% order the points.  For every point in the original profile, add it and
% all of its revolved copies to the list (in CCW order), then move on to
% the next profile point and its copies.  This works because bp is already
% ordered to form a continuous polyline (pt 1 is connected to pt 2, pt 2
% is connected to pt 3, ..., pt n-1 is connected to pt n)
pts3 = [];
for p=1:size(pts,1)
    a = (nSpokes-1)*(p-1)+1;
    b = a + (nSpokes-1)-1;
    
    pts3 = [pts3; pts(p,:); pts2(a:b, :)];
end
pts = pts3;


% Build the faces.  For each face, start at the bottom left point and work
% CCW
%
% eg:
% D --- C
%       |
%       |
% A --- B

% march each bf and make a face to the one rotated.  Basically B is index
% A+1, C is above B (index i+nSpokes), D is above A.  Use Modulus so points
% wrap around (if A is pt n, B is pt 1)
%
% TODO:
% this still needs some work on the wrapping, since it forms some weird
% crossing pattern
nPts = size(pts, 1);
tri = [];
for p = 1:size(bf,1)
    for s=1:nSpokes
        A = s+(p-1)*nSpokes;
        B = (mod(A, nSpokes)+1)+(p-1)*nSpokes;
        C = mod(A + nSpokes, nPts) + 1;
        D = mod(C - 2, nPts) + 1;
        
        % list of the pt indices that are connected together in order.
        % This forms a closed polygon.  Using 4 pts makes a tetragon.
        tri = [tri; [A B C D]];
    end
end


