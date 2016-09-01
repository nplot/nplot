function tri = delaunayfromvoronoi(faces, coords, dims, nv)

% Deduce a triangulation from a known Voronoi diagram
% Use the fact that each Voronoi vertex is the center of a circle (sphere)
% that contains the corners of one delaunay triangle. Each Voronoi cell
% corresponds to one cornerpoint 
%
% "faces" does not have to define a real Voronoi diagram, i.e. cells may be
% cut off etc., as this is usually the case, either for too distant data
% points, and (practically always) at the borders of the point set.
% In this case, "tri" is no complete Delaunay triangulation - roughly speaking, 
% it does not contain triangles connecting points that are too far away from 
% each other, which is convenient for later use (interpolation etc.).
% Apart that, "tri" is an exact solution.
%
% Works currently in two dimensions, but can be generalized
%
% P. Steffens, 06/2010



if nargin<3, dims = 2; end  % ** !
if nargin<4, nv = max(faces(:)); end  % Take this as the number of vertices 
%(if faces is given as cell array like output of voronoin, then need to provide nv!)
facecount = zeros(nv,1);  % for each vertex, count to how many faces it belongs
tri_raw = nan(nv,dims+1); % a line in tri_raw contains the numbers of the faces that contain it


if ~iscell(faces)   % standard (matrix) format of faces
    for i = 1:size(faces,1)
        for j = 1:size(faces,2)
            fij = faces(i,j);
            if ~isnan(fij)
                facecount(fij) = facecount(fij) + 1;
                tri_raw(fij, facecount(fij)) = i;
            end
        end
    end
    
else

    for i = 1:length(faces)
        faci = faces{i};
        for j = 1:numel(faci)
            fij = faci(j)-1;
            if fij > 0
                facecount(fij) = facecount(fij) + 1;
                tri_raw(fij,facecount(fij)) = i;
            end
        end
    end
end


% lines in tri_raw with less than dims+1 entries do not correspond to triangles
% lines in tri_raw with more than -"- need to be split in triangles (more
% than dims+1 on the same circle, so there is some freedom in the choice)


extracount = sum(facecount(facecount>dims+1)-dims-1); % to estimate number of triangles
tri = zeros(nv-dims+extracount,dims+1); % initialize with estimate for number
tcount = 1;
for i = 1:size(tri_raw,1)
    if facecount(i) == dims+1
        tri(tcount, :) = tri_raw(i, 1:dims+1);
        tcount = tcount+1;
    elseif facecount(i) > dims+1

        % In 2D, do the following:
        % sort points by x-coord., and add a triangle for each, one after
        % the other. Each time adding a point creates two new sides. To exactly one of
        % these, the following triangle will be connected.
        
        [sx, xorder] = sort(coords(tri_raw(i,1:facecount(i)),1));
        
        tri(tcount,:) = tri_raw(i, xorder(1:3));
        tcount = tcount + 1;
        
        for k= (dims+2):facecount(i)
            
            facet1 = coords(tri(tcount-1,[2,3]),:);
            % facet2 = coords(tri(tcount-1,[1,3]),:);
            
            % facet normal = [- (facet_y2-facet_y1) ; facet_x2-facet_x1];
            % Sign of Distance of next point to facet: (nextpoint - facetpoint) * normal
            % So test if this has a different sign than the distance of the opposite(old) triangle point
            % If yes, connect new triagle to this facet. If not, connect to the other.
            
            % Remark: do not need to consider the "ratios" here!
            
            normal1 = [ - (facet1(2,2)-facet1(1,2)) ; facet1(2,1)-facet1(1,1) ];
            if ((coords(tri_raw(i,xorder(k)),:) - coords(tri(tcount-1,3),:)) * normal1 ) * ((coords(tri(tcount-1,1),:) - coords(tri(tcount-1,3),:)) * normal1 ) < 0
                % coords(tri_raw(i,xorder(k)),:) is the new point
                % coords(tri(tcount-1,3),:) is one of the facet's points
                % coords(tri(tcount-1,1),:) is the third point of the last triangle (not on the facet)
                % So, if the new point lies on the correct(i.e. other) side: (product<0)
                tri(tcount,:) = [tri(tcount-1, [2,3]), tri_raw(i, xorder(k)) ];
            else
                tri(tcount,:) = [tri(tcount-1, [1,3]), tri_raw(i, xorder(k)) ];
            end
            tcount = tcount +1;
        end
        
        % In higher dimensions, the simplification by sorting will not
        % work, even the number of (hyper-)tetrahedra is not constant.
        % Generalize such that:
        % - start with a tetrahedron
        % - maintain list of "active" facets, i.e. facets to which no other tetrahedron has been connected
        % - take an active facet and look if there is a point that can be connected to it (i.e. lying in the "outside" half-plane)
        % (look at http://www.iue.tuwien.ac.at/phd/fleischmann/node61.html)

    end
end

tri = tri(1:tcount-1, :);  % delete empty lines

        


