function contourindices = findcontour(faces)

% Find border of polygon set (= edges that belong to only one polygon)
% in 2 dimensions (coords in arbitrary dimensions)
% contourindices: ordered list of vertices (1d-array of indices) that are on the borderline
% note: several borders may exist (for holes or unconnected sets), cell array is returned 

% P. Steffens, 8/2014


finind = logical(isfinite(faces));
edgecount = sum(finind,2); % (#edges of each polygon)

faces=faces';
faces = [faces(size(faces,1)*(0:size(faces,2)-1)+edgecount'); faces];
% now faces is transposed (easier indexing); first row contains last vertex of each face (again)

edgelist = [faces([false(1,size(faces,2));finind']), faces([finind';false(1,size(faces,2))])];
% linear list, each line is one edge

% now, eliminate doubles (by pairs!)
edgelist = sortrows(sort(edgelist,2));
% % for eind=1:sum(edgecount)-1
% %     if all(edgelist(eind,:)==edgelist(eind+1,:))
% %         edgelist(eind:eind+1,:)=nan;
% %     end
% % end

% %(Replacement of loop:)
rmind = false(size(edgelist,1),1);
eqind = (edgelist(1:end-1,1)==edgelist(2:end,1)) & (edgelist(1:end-1,2)==edgelist(2:end,2)); % lines equalling the following line
while any(eqind)
    lpind = eqind(1:end-1) & ~eqind(2:end);   % indicates last pair of sequence
    rmind([lpind; 0] | [0; lpind]) = true;    % mark this and the following line for removal
    eqind(lpind | [lpind(2:end); 0]) = false; % remove these and the preceding line from eqind
end
edgelist(rmind,:)=nan;

edgelist = edgelist(isfinite(edgelist(:,1)),:);
edgelist(edgelist(:,1)==edgelist(:,2),:)=nan; % remove zero length edges (if improper polygon structure) 

% now, bring them in the right order
nedg = size(edgelist,1);
% contourindices = zeros(nedg,1);
ctind = 0;

while any(isfinite(edgelist(:,1)))
    % start a new contourline (there may be several (if holes etc.))
    nextedge = find(isfinite(edgelist(:,1)),1,'first'); % first remaining edge in list
    ctind = ctind+1;
    contourindices{ctind}(1:2) = edgelist(nextedge,:); %#ok<*AGROW>
    edgelist(nextedge,:) = nan;  % mark "used" edges as nan
    eind = 3;
    % now, continue contourline as long as further connected edges exist
    while any(any(edgelist == contourindices{ctind}(eind-1)))
        nextedge = mod(find(edgelist == contourindices{ctind}(eind-1), 1, 'first')-1,nedg)+1; % row index of next edge
        contourindices{ctind}(eind) = edgelist(nextedge,  edgelist(nextedge,:)~= contourindices{ctind}(eind-1)); % in this row, take the one which is different
        edgelist(nextedge,:)=nan;
        eind=eind+1;
    end
end




    