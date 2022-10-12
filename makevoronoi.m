function [vertexlist, faces, triang, vol] = makevoronoi(points,maxsize,ratio,switches)

% [vertexlist, faces, triang, vol] = makevoronoi(points,maxsize,ratio,switches)
%
% Calculates the voronoi cells for N-dim data points. The maximum size of the
% cells is maxsize in each dimension
% ratio scales the coordinate axes (metrics)
% switches: may contain 'calcvol' and/or 'calctri' to compute volume / delaunay triangulation
%
% vertexlist:    rows are the coordinates of all vertices
% faces:         rows are lists of the vertices (in vertexlist) defining each cell
% triang:        Delaunay triangulation
% vol:           volume / surface area of the cells

% replaces calcpatch.m as generalization to higher dimensions
%
% P. Steffens, 08/2014
%

ndims = size(points,2);     % Number of dimensions
npoints = size(points,1);   % Number of (real) points

calctriang  = (nargin>3) && any(strcmpi(switches,'calctri'));    % Calculate Delaunay triangulation?
% calcvol     = (nargin>3) && any(strcmpi(switches,'calcvol'));    % Calculate Volume of cells?     % ** so far unused


% for each dimension, add two artificial points to make all cells finite
points2 = points;
meanpoint = mean(points,1);
for d = 1:ndims
    s = max(points(:,d)) - min(points(:,d));
    points2 = [points2; meanpoint(1:(d-1)), max(points(:,d)) + s + 3*maxsize(d), meanpoint((d+1):end)]; %#ok<*AGROW>
    points2 = [points2; meanpoint(1:(d-1)), min(points(:,d)) - s - 3*maxsize(d), meanpoint((d+1):end)];
end

% Scale axes
if nargin<3, ratio = ones(1,ndims); end
for d = 1:ndims
    points2(:,d) = points2(:,d) * ratio(d);
end

%compute Voronoi vertices V and cells C
[V,C] = voronoin(points2,{'Qz'}); % ,{'Qbb','Qz'});

% Scale back
for d = 1:ndims
    V(:,d) = V(:,d) / ratio(d);
end


%% new
% Alternative:  % no Vol calculation
vertexlist = V;
faces = zeros(npoints,4); 
facecount = zeros(npoints,1);
for nf=1:npoints 
    thisface = C{nf};
    facecount(nf)=numel(thisface);
    faces(nf,1:facecount(nf)) = thisface;
end
helpfaces = faces; % for coordinates, create a helper matrix with finite entries only
for c = 2:size(faces,2), helpfaces(helpfaces(:,c)==0,c) = helpfaces(helpfaces(:,c)==0,1); end
% needcut = false(npoints,ndims); % index for which faces need cutting
needcut = false(npoints,1); % index for which faces need cutting
for d=1:ndims
    allcoord = vertexlist(:,d);
    thisdimcoord = allcoord(helpfaces);
%     needcut(:,d) = (max(thisdimcoord,[],2) > points(:,d) + maxsize(d)/2) | (min(thisdimcoord,[],2) < points(:,d) - maxsize(d)/2);
    needcut = needcut | (max(thisdimcoord,[],2) > points(:,d) + maxsize(d)/2) | (min(thisdimcoord,[],2) < points(:,d) - maxsize(d)/2);
end


% % **
 faces(faces==0)=NaN;
triang = delaunayfromvoronoi(faces, points, 2);
% % **

nvert = size(vertexlist,1)+1;
vertexlist = [vertexlist; nan( sum(any(needcut,2))* size(faces,2), ndims)]; % some preallocation by guessing

% for nc = find(any(needcut,2))'
for nc = find(needcut)'
    % Cut the cells that need to be cut
    newcell = V(faces(nc,1:facecount(nc)),:);
    for d=1:ndims
        if max(newcell(:,d)) > points(nc,d) + maxsize(d)/2
            if ndims==2, K = [(1:size(newcell,1))', [(2:size(newcell,1))';1]]; else K=[]; end
            [~, newcell] = polyhypercut( newcell,  [zeros(1,d-1), 1, zeros(1,ndims-d)],   points(nc,d) + maxsize(d)/2,    K );
%             cutpoints = [cutpoints; addcutpoints];
%             cut= true;
        end
        if min(newcell(:,d)) < points(nc,d) - maxsize(d)/2
            if ndims==2, K = [(1:size(newcell,1))', [(2:size(newcell,1))';1]]; else K=[]; end
            [~, newcell] = polyhypercut( newcell, -[zeros(1,d-1), 1, zeros(1,ndims-d)], -(points(nc,d) - maxsize(d)/2),   K );
%             cutpoints = [cutpoints; addcutpoints];
%             cut = true;
        end
    end
    
%     newcell = unique(round(newcell*1E10)/1E10,'rows','stable'); % delete evtl. doubles %** do this at the end!!
    
%     if  ndims ~= 2 %optimize and get Volume
%         %(treat ndims=2 below)
%         [K,vol(nc)] = convhulln(newcell);
%         ind = unique(K);
%         newcell = newcell(ind,:);
%     end
    
    
    % Append vertices of this cell to vertexlist
    vertexlist( nvert : (nvert+size(newcell,1)-1), 1:ndims)  = newcell;
    nvert = nvert + size(newcell,1);
    
    % Append indices of these vertices to cell list (undefined order, not safe for 2D!)
%     oldfaces = faces(i,:);
    faces(nc,1:size(newcell,1))= (nvert-size(newcell,1)) : (nvert-1);
    faces(nc,size(newcell,1)+1:end) = nan;
    
end


vertexlist = vertexlist(1:nvert-1,:);
faces(faces==0)=NaN;


%% old

% % % count vertices to preallocate memory for vertexlist and faces
% % % (however, not really predictable due to cutting below)
% % % Note also that this largely overcounted, because it corresponds to a
% % % non-optimized list (multiple occurences), before "optimizepatch"
% % nvert=0;
% % for i=1:length(C)
% %     nvert = nvert + numel(C{i});
% % end
% % % nvert = size(V,1) - 1;  % ("-1" because first point is "inf")
% % vertexlist=zeros(nvert,ndims); 
% % faces = zeros(npoints,4);
% % vol = zeros(npoints,1);
% % 
% % 
% % % % Obtain Triangulation, if desired   % --> later below!
% % % % if calctriang 
% % %     if ndims == 2
% % %         triang = delaunayfromvoronoi({C{1:npoints}}, points, 2, nvert);
% % % %     else
% % %     elseif calctriang
% % %         fprintf('Cannot do triangulation so far for dimension higher than 2!\n'); return;
% % %     end
% % % % end
% % 
% % 
% % nvert = 1;
% % for i=1:npoints    
% %     cut = false;
% % %     cutpoints = [];
% %     newcell = V(C{i},:);
% %     % Find out if cell needs to be cut - do this for each dimension
% %     for d=1:ndims
% %         if max(newcell(:,d)) > points(i,d) + maxsize(d)/2
% %             if ndims==2, K = [(1:size(newcell,1))', [(2:size(newcell,1))';1]]; else K=[]; end
% %             [~, newcell] = polyhypercut( newcell,  [zeros(1,d-1), 1, zeros(1,ndims-d)],   points(i,d) + maxsize(d)/2,    K );
% % %             cutpoints = [cutpoints; addcutpoints];
% %             cut= true;
% %         end
% %         if min(newcell(:,d)) < points(i,d) - maxsize(d)/2
% %             if ndims==2, K = [(1:size(newcell,1))', [(2:size(newcell,1))';1]]; else K=[]; end
% %             [~, newcell] = polyhypercut( newcell, -[zeros(1,d-1), 1, zeros(1,ndims-d)], -(points(i,d) - maxsize(d)/2),   K );
% % %             cutpoints = [cutpoints; addcutpoints];
% %             cut = true;
% %         end
% %     end
% %     
% %     if cut && (ndims ~= 2) %optimize and get Volume
% %         %(treat ndims=2 below)
% %         [K,vol(i)] = convhulln(newcell);
% %         ind = unique(K);
% %         newcell = newcell(ind,:);
% %     end
% %     
% %     
% %     % Append vertices of this cell to vertexlist
% %     vertexlist( nvert : (nvert+size(newcell,1)-1), 1:ndims)  = newcell;
% %     nvert = nvert + size(newcell,1);
% %     
% %     % Append indices of these vertices to cell list (undefined order, not safe for 2D!)
% % %     oldfaces = faces(i,:);
% %     faces(i,1:size(newcell,1))= (nvert-size(newcell,1)) : (nvert-1);
% %     
% %     
% % end
% % 
% % faces(faces==0)=NaN;

%%
%**
% beforeok = checkfaces(faces);
[vertexlist, faces] = optimizepatch(vertexlist, faces);

% Remove eventual doubles in faces that have been cut (replaces "unique", but detects only those following each other)
cutlines = find(needcut);
for fc=size(faces,2):-1:2
    rmind =  cutlines(faces(needcut,fc) == faces(needcut,fc-1)); % lines in which a double is in column fc
    faces(rmind,fc-1:end-1) = faces(rmind,fc:end);
    faces(rmind,end) = nan;
end

%**
% afterok = checkfaces(faces);
% if beforeok && ~afterok, 
%     fprintf('!!!!!****\n'); 
% end

%% end here for dimensions > 2
if ndims>2, return; end

%%
% triang = delaunayfromvoronoi(faces, points, 2);

%%
% Now, check the following: it may happen that a facet(edge) which was shared with
% one of the neighbors has been cut, and a new vertex is now lying on this
% neighbor's facet. This vertex should now be incorporated in that facet.
% (For now, do this in 2D only.) **
% (One could in principle also simplify this by keeping more information from the cutting...**)




% as candidates, look for vertices that appear in only one polygon
onceverts = find(histc(faces(:),.5:max(faces(:))+.5) == 1)';
for vertind = onceverts
    owner = mod(find(faces==vertind)-1,size(faces,1))+1; % this is the face to which it belongs
    % find the neighbors of it
    nbs = triang(any(triang == owner,2),:);
    nbs = unique(nbs(:)); nbs = nbs(nbs~=owner);
    % among them only those having the same distance to vert(vertind) as owner are candidates
    nbs = nbs(abs( (points(nbs,1)-vertexlist(vertind,1)).^2   + (points(nbs,2)-vertexlist(vertind,2)).^2 ...
                  -(points(owner,1)-vertexlist(vertind,1)).^2 - (points(owner,2)-vertexlist(vertind,2)).^2 ) < 1e-8);
    for nbind = nbs'
        % check the edges of this neighbor if vert(vertind) lies on them
        point1 = vertexlist(faces(nbind,1:sum(isfinite(faces(nbind,:)))),:);
        point2 = point1([end,1:end-1],:); % start and end point of each edge
        
        vec1x = point1(:,1)-vertexlist(vertind,1); vec1y = point1(:,2)-vertexlist(vertind,2); norm1sq = vec1x.*vec1x + vec1y.*vec1y;
        vec2x = point2(:,1)-vertexlist(vertind,1); vec2y = point2(:,2)-vertexlist(vertind,2); norm2sq = vec2x.*vec2x + vec2y.*vec2y;
        % points are between, if scalar prod of vec1, vec2 = -|v1|*|v2|, and v1,v2 not too short (not too near p1,p2)
        sprod = vec1x.*vec2x + vec1y.*vec2y;
        edgind = sprod < 0 & (norm1sq.*norm2sq - sprod.^2) < 1e-15 & norm1sq > 1e-15 & norm2sq > 1e-15;
        if sum(edgind)==1
            faces(nbind,find(edgind):sum(isfinite(faces(nbind,:)))+1) = [vertind, faces(nbind,find(edgind):sum(isfinite(faces(nbind,:))))];
            faces(faces==0)=nan; %(if matrix has been extended)
        end
    end
end
        


%%
if ndims==2 && nargout>3
    vol = abs(polygonarea( vertexlist(:,1), vertexlist(:,2), faces));
end





%%


% % Check the 
%     if cut && ndims==2
%         % identify edges shared with any other polygon (possibly affected)
%         oldfaces = [oldfaces(oldfaces>0),oldfaces(1)];
%         for nf=1:numel(oldfaces-1)
%             % is this a shared edge?
%             ind = any(faces==oldfaces(nf),2) & any(faces==oldfaces(nf+1),2); % (do not need to check order because of convexity)
%             if ~any(ind), continue; end
%             point1 = vertexlist(oldfaces(nf),:); point2 = vertexlist(oldfaces(nf+1),:); % end points of this edge
%             % which cutpoints lie on this edge
%             vec1x = point1(1)-cutpoints(:,1); vec1y = point1(2)-cutpoints(:,2); norm1sq = vec1x.*vec1x + vec1y.*vec1y;
%             vec2x = point2(1)-cutpoints(:,1); vec2y = point2(2)-cutpoints(:,2); norm2sq = vec1x.*vec1x + vec1y.*vec1y;
%             % scalar prod of vec1, vec2 must be -|v1|*|v2| (on the connection line), and v1,v2 not too short (not too near p1,p2)
%             sprod = vec1x.*vec2x + vec1y.*vec2y;
%             cpind = sprod < 0 & (norm1sq.*norm2sq - sprod.^2) < 1e-15 & norm1sq > 1e-15 & norm2sq > 1e-15;
%             if ~any(cpind), continue; 
%             elseif sum(cpind)==1    % put this point in between
%                 for insertind = find(ind') %(should be only one!)
%                     pos = min([find(faces(insertind,:)==oldfaces(nf)), find(faces(insertind,:)==oldfaces(nf+1))]);
%                     faces(insertind, (pos:sum(faces(insertind,:)>0))+1) = [nvert+1, faces(insertind, pos:sum(faces(insertind,:)>0))];
%                     vertexlist = [vertexlist; cutpoints(cpind,:)];
%                     nvert = nvert + 1;
%                 end
%             else    % insert more than one point (in right order)
%                 for insertind = find(ind')
%                     ch = convhull([vertexlist(oldfaces(1:end-1),:);cutpoints(cpind,:)]); %make it easy...
%                     faces(insertind, 1:numel(ch)) = nvert + ch(:)';
%                     vertexlist = [vertexlist; vertexlist(oldfaces(1:end-1),:);cutpoints(cpind,:)];
%                     nvert = nvert + numel(ch);
%                 end
%             end
%         end
%     end




