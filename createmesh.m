function [vertices, faces, numbers] = createmesh (vertexlist, cells, vectors, origin)

% Obtain an n-dim subspace (typically n=2) of an N-dim cell set
% (defined by vertexlist, cells)
% Plane is defined by vectors (rows are vectors) and origin
% Returned coordinates of vertices refer to the given vectors (coefficients for linear combination):
% (N-dim) coordinates in original system may be obtained via origin+coords*vectors
% numbers are the indices of the used cells
%
% P. Steffens, 08/2008


ndims = size(vectors,2);
ndimneu = size(vectors,1);
origin = origin(:)';


% Create N-n hyperplanes whose intersection defines the subspace

% Obtain N-n vectors that are orthogonal to the given 'vectors'
% "Guessing" by taking N vectors and regarding their orthogonal parts

testvecs = eye(ndims);
for d = 1 : ndims
    for dn = 1:ndimneu
        testvecs(d,:) = testvecs(d,:) - testvecs(d,:) * vectors(dn,:)' / sum(vectors(dn,:).^2) * vectors(dn,:);
    end
end
[s, order] = sort(sum(testvecs.^2),2,'descend');

normals = testvecs(order(1:(ndims-ndimneu)),:); 
nhyp = ndims-ndimneu;
sortedges = false;

numbers = 1:(size(cells,1));
faces   = cells;
vertices = vertexlist;

for h = 1:nhyp
    if ndimneu == 2 && h==nhyp, sortedges = true; end 
    const = normals(1,:)*origin';
    % normals(h,:) and const define the hyperplane h
    cellnumbers = numbers;
    cells = faces;
    vertexlist = vertices;
    numbers = [];
    vertices = [];
    faces = [];
    vol = [];
    nfac = 0;
    nvert = 0;
    for c = 1:numel(cellnumbers)
        thiscell = cells(c,:);
        thiscell = thiscell(isfinite(thiscell));
        thiscellvertex = vertexlist(thiscell,:);
        % Now cut:
        cutpoints = polyhypercut(thiscellvertex, normals(1,:), const);

        
        if ~isempty(cutpoints)
            numbers = [numbers; cellnumbers(c)];
            % transform to projected coordinates, reduce dimensions by one
            str.coordlist = cutpoints;
            str = coordtransform(str, 'projection', [vectors; normals(2:end,:)], origin); 
            cutpoints = str.coordlist;
            [K,vol(c)] = convhulln(cutpoints);
            if sortedges  % in 2D, sort vertices in continuous order
                pointorder = [K(1,:), zeros(1,size(K,1)-2)];
                lastpoint = K(1,2);
                K = K(2:end,:);
                for i = 2:(size(K,1))
                    ind = find(K'==lastpoint); 
                    % ind = ind(1);
                    row = ceil(ind/2);
                    lastpoint = K(row,mod(ind,2)+1);
                    pointorder(i+1) = lastpoint;
                    K = K([1:(row-1),(row+1):end],:);
                end
                cutpoints = cutpoints(pointorder,:);
                nnewpoints = i+1;
            else % order does not matter
                ind = unique(K);
                cutpoints = cutpoints(ind,:);
                nnewpoints = numel(ind);
            end
            vertices = [vertices; cutpoints];
            nfac = nfac +1;
            faces(nfac,1:nnewpoints) = (nvert + 1) : (nvert + nnewpoints);
            nvert = nvert + nnewpoints;
        end
    end
    
    faces(faces==0) = NaN;
    
    % As all other coords have been transformed, need to transform also
    % remaining testvecs and origin
    str.coordlist = [vectors + repmat(origin,size(vectors,1),1); normals(2:end,:) + repmat(origin,size(normals,1)-1,1); origin];
    str = coordtransform(str, 'projection', [vectors; normals(2:end,:)], origin);
    vectors = str.coordlist(1:ndimneu,:);
    normals = str.coordlist((ndimneu+1):(end-1),:);
    origin = str.coordlist(end,:);
    
    
end


[vertices, faces] = optimizepatch(vertices, faces);


    