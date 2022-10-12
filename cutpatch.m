function [newfaces, newvertices, cellassignment, partition] = cutpatch(faces, vertices, n, c, d2)

% Cut the patch defined by (faces, vertices) along the hyperplane n*x=c
% (where n is a vector (the normal vector) in n-dim space and c a constant)
% i.e. cells that are intersected by the hyperplane are divided along this
% line
% cellassignment contains the assignment of the new faces to the original
% ones (use for instance for color assignment)
% partition (logical) is true for faces above the hyperplane (n*x>c), false otherwise 
%
% If d2=true given, use algorithm for 2D (flat patches in higher dimensions)
%
% P. Steffens 08/2008



n = n(:);
ndim = numel(n);
nvert = size(vertices,1);
if size(vertices,2)~=ndim, fprintf('Dimensional mismatch in cutpatch.\n'); return; end
K = [];

vertexpos = zeros(nvert,1);                       % Calculate for each vertex its position w.r.to hyperplane
for nd = 1:ndim
    vertexpos = vertexpos + vertices(:,nd) * n(nd);
end
vside = sign(vertexpos -c);
vsidepos = [vside>0; true];
vsideneg = [vside<0; true];
testfaces = faces;
testfaces(isnan(testfaces)) = nvert + 1;

if size(faces,1) == 1  %strange, here "all" behaves in a bizarre way...
    allpos = all(vsidepos(testfaces));
    allneg = all(vsideneg(testfaces));
else
    allpos = all(vsidepos(testfaces),2);
    allneg = all(vsideneg(testfaces),2);
end


indexoneside = find(allpos | allneg);  % These faces are completely on one side of the cut 
indextwoside = find(~allpos & ~allneg); % These are the faces to be cut
ncuts = numel(indextwoside);


% Estimation: each face that is really cut produces 2 new vertices 
% (this could be made more intelligently by considering the shared vertices),
% and one additional new face
newvertices = [vertices ; zeros( 2*(size(faces,1)-numel(indexoneside)), size(vertices,2))] ;
%newfaces = zeros( size(faces,1) + ncuts, size(faces,2));
cellassignment = zeros(size(faces,1) + ncuts, 1);
partition      = zeros(size(faces,1) + ncuts, 1);

% Assign those that are unchanged
newfacecount = numel(indexoneside);
%newfaces(1:newfacecount,:)     = faces(indexoneside,:);
cellassignment(1:newfacecount) = indexoneside;
partition(1:newfacecount)      = allpos(indexoneside);
newfaces = [faces(indexoneside,:); NaN(2*ncuts, size(faces,2))];
faccol = size(newfaces,2); %Number of columns in newfaces


if ndim>2 && (nargin<5 || ~d2)
    % use polyhypercut to treat higher dimensions

    for nface = indextwoside'

        fa = faces(nface, isfinite(faces(nface,:)));    % vertices defining the current face

        %if ndim==2, K = [(1:numel(fa))', [(2:numel(fa))';1]]; end
                                                        % in 2D, use the correct order for convex hull
        [cutpoints, poly1, poly2, np1, np2] = polyhypercut( vertices(fa,:), n, c);

        ncp = size(cutpoints,1);                        % Number of cutpoints
        newvertices((nvert+1):(nvert+ncp),:) = cutpoints;%Add cutpoints

        faind = [fa, (nvert+1):(nvert+ncp)];            % Vertices defining the current face including cut points

        if ~isempty(poly1)                              % Add new cell below plane
            newfacecount = newfacecount+1;
            if numel(np1)>faccol, newfaces(:,(faccol+1):numel(np1))=NaN; faccol = numel(np1); end % Ensure Nans in empty entries
            newfaces(newfacecount, 1:numel(np1)) = faind(np1);
            cellassignment(newfacecount) = nface;
            partition(newfacecount) = false;
        end

        if ~isempty(poly2)                              % Add new cell above plane
            newfacecount = newfacecount+1;
            if numel(np2)>faccol, newfaces(:,(faccol+1):numel(np2))=NaN; faccol = numel(np2); end
            newfaces(newfacecount, 1:numel(np2)) = faind(np2);
            cellassignment(newfacecount) = nface;
            partition(newfacecount) = true;
        end

        nvert = nvert+ncp;

    end
    
    
else   % in 2D, use more efficient (old) code
    
    faneu{1}=[0,0,0,0,0];  % initialize two vertex lists for the eventually two polygons (for use within loop)
    faneu{2}=[0,0,0,0,0];

    % Loop over all faces that are to be cut
    for nface = indextwoside'

        newfacecount = newfacecount+1;                  % index in newfaces
        fa = faces(nface, isfinite(faces(nface,:)));    % vertices defining the current face
        fa = fa([1:end,1]);
        facount = [0,0];                                % actual length of faneu{1} and {2}
        whichfa = 1;                                    % on which of the two currently working?
        sidefound = false;                              % side on which polygon lies already determined? (for partition)
        lastside = vside(fa(1));                        % lastside, newside: on which side are vertex nv and nv+1?
        for nv=1:(numel(fa)-1)
            newside = vside(fa(nv+1)); if newside == 0, newside=lastside; end
            if lastside == 0, lastside = newside; end
            if ~sidefound && whichfa==1 && vsidepos(fa(nv)) && ~vsideneg(fa(nv)), partition(newfacecount)=1; sidefound=true; end
            if ~sidefound && whichfa==1 && ~vsidepos(fa(nv)) && vsideneg(fa(nv)), partition(newfacecount)=0; sidefound=true; end
            facount(whichfa) = facount(whichfa) + 1;
            faneu{whichfa}(facount(whichfa)) = fa(nv);  % add point to current list                                            

            if lastside ~= newside  % There is an intersection between vertex No. nv and nv+1
                v1 = vertices(fa(nv),:);     dv = vertices(fa(nv+1),:) - v1;
                lambda = (c - v1 * n) / (dv * n); % calculate intersection point  
                cutcoord = v1 + lambda * dv ;
                nvert = nvert+1;
                newvertices(nvert,:) = cutcoord;        % Create new entry in vertex list
                                                        % Add intersection point to both polygons...
                facount = facount + 1; % (increment both)
                faneu{whichfa}(facount(whichfa)) = nvert; 
                whichfa = 3 - whichfa;                  % Switch to other polygon (1->2 or 2->1)
                faneu{whichfa}(facount(whichfa)) = nvert;
            end  
            lastside = newside;
        end
                        %Update all entries for the one or two obtained faces
        if facount(1)>faccol, newfaces(:,(faccol+1):facount(1))=NaN; faccol = facount(1); end
        newfaces(newfacecount,1:facount(1)) = faneu{1}(1:facount(1));
        cellassignment(newfacecount) = nface;
        if ~isempty(faneu{2})
            newfacecount=newfacecount+1; 
            if facount(2)>faccol, newfaces(:,(faccol+1):facount(2))=NaN; faccol = facount(2); end
            newfaces(newfacecount,1:facount(2)) = faneu{2}(1:facount(2)); 
            cellassignment(newfacecount) = nface;
            partition(newfacecount) = ~partition(newfacecount-1);
        end
    end
    
    
end
    

newvertices    = newvertices(1:nvert,:);
newfaces       = newfaces(1:newfacecount,:);
cellassignment = cellassignment(1:newfacecount);
partition      = partition(1:newfacecount);

% newfaces(newfaces==0) = nan;

return;




% ******************************************************
% Code below is old version; will never run
% ******************************************************
% ******************************************************



n = n(:); %#ok<UNRCH>
ndim = numel(n);
nvert = size(vertices,1);
if size(vertices,2)~=ndim, fprintf('Dimensional mismatch in cutpatch.\n'); return; end


vside = zeros(size(vertices,1),1);                       % Calculate for each vertex its position w.r.to hyperplane
for nd = 1:ndim
    vside = vside + vertices(:,nd) * n(nd);
end
vside = sign(vside -c);
vsidepos = [vside>0; true];
vsideneg = [vside<0; true];
testfaces = faces;
testfaces(isnan(testfaces))=size(vertices,1)+1;

if size(faces,1) == 1  %strange, here "all" behaves in a bizarre way...
    allpos = all(vsidepos(testfaces));
    allneg = all(vsideneg(testfaces));
else
    allpos = all(vsidepos(testfaces),2);
    allneg = all(vsideneg(testfaces),2);
end

indexoneside = find(allpos | allneg);  % These faces are completely on one side of the cut 
indextwoside = find(~allpos & ~allneg); % These are the faces to be cut
ncuts = numel(indextwoside);


% Estimation: each face that is really cut produces 2 new vertices 
% (this could be made more intelligently by considering the shared vertices),
% and one additional new face
newvertices = [vertices ; zeros( 2*(size(faces,1)-numel(indexoneside)), size(vertices,2))] ;
%newfaces = zeros( size(faces,1) + ncuts, size(faces,2));
cellassignment = zeros(size(faces,1) + ncuts, 1);
partition      = zeros(size(faces,1) + ncuts, 1);

% Assign those that are unchanged
newfacecount = numel(indexoneside);
%newfaces(1:newfacecount,:)     = faces(indexoneside,:);
cellassignment(1:newfacecount) = indexoneside;
partition(1:newfacecount)      = allpos(indexoneside);
newfaces = [faces(indexoneside,:); NaN(2*ncuts, size(faces,2))];

faneu{1}=[0,0,0,0,0];  % initialize two vertex lists for the eventually two polygons (for use within loop)
faneu{2}=[0,0,0,0,0];

% Loop over all faces that are to be cut
for nface = indextwoside'
    
    newfacecount = newfacecount+1;                  % index in newfaces
    fa = faces(nface, isfinite(faces(nface,:)));    % vertices defining the current face
%     if all(vside(fa)<=0) || all(vside(fa)>=0)       % test if all vertices on same side
%         newfaces(newfacecount,1:(length(fa))) = fa([1:end]);   % ** removed a '+1' here on both sides...
%         cellassignment(newfacecount) = nface;
%         partition(newfacecount) = any(vside(fa)>0);
%         continue;                                   % if so, do not need to do the following, polygon stays the same
%     end
    fa = fa([1:end,1]);
    facount = [0,0];                                % actual length of faneu{1} and {2}
    whichfa = 1;                                    % on which of the two currently working?
    sidefound = false;                              % side on which polygon lies already determined? (for partition)
    lastside = vside(fa(1));                        % lastside, newside: on which side are vertex nv and nv+1?
    for nv=1:(numel(fa)-1)
        newside = vside(fa(nv+1)); if newside == 0, newside=lastside; end
        if lastside == 0, lastside = newside; end
        if ~sidefound && whichfa==1 && vsidepos(fa(nv)) && ~vsideneg(fa(nv)), partition(newfacecount)=1; sidefound=true; end
        if ~sidefound && whichfa==1 && ~vsidepos(fa(nv)) && vsideneg(fa(nv)), partition(newfacecount)=0; sidefound=true; end
        facount(whichfa) = facount(whichfa) + 1;
        faneu{whichfa}(facount(whichfa)) = fa(nv);  % add point to current list                                            
        
        if lastside ~= newside  % There is an intersection between vertex No. nv and nv+1
            v1 = vertices(fa(nv),:);     dv = vertices(fa(nv+1),:) - v1;
            lambda = (c - v1 * n) / (dv * n); % calculate intersection point  
            cutcoord = v1 + lambda * dv ;
            nvert = nvert+1;
            newvertices(nvert,:) = cutcoord;        % Create new entry in vertex list
                                                    % Add intersection point to both polygons...
            facount = facount + 1; % (increment both)
            faneu{whichfa}(facount(whichfa)) = nvert; 
            whichfa = 3 - whichfa;                  % Switch to other polygon (1->2 or 2->1)
            faneu{whichfa}(facount(whichfa)) = nvert;
        end  
        lastside = newside;
    end
                    %Update all entries for the one or two obtained faces
    newfaces(newfacecount,1:facount(1)) = faneu{1}(1:facount(1));
    cellassignment(newfacecount) = nface;
    if ~isempty(faneu{2})
        newfacecount=newfacecount+1; 
        newfaces(newfacecount,1:facount(2)) = faneu{2}(1:facount(2)); 
        cellassignment(newfacecount) = nface;
        partition(newfacecount) = ~partition(newfacecount-1);
    end
end

newvertices    = newvertices(1:nvert,:);
newfaces       = newfaces(1:newfacecount,:);
cellassignment = cellassignment(1:newfacecount);
partition      = partition(1:newfacecount);

%newfaces(newfaces==0) = nan;



 