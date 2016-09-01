function neighbors = findneighbors(nface, coordlist, deltri)

% Find neighbors of face 'nface' (index to 'coordlist')
% patch implicitly defined by 'coordlist'
%
% Idea: Delaunay-triang. is dual to voronoi-diagram.
% The neighbors are thus the Voronoi cells that share a facet with the ones
% in nface.

% P. Steffens, 05/2008 - 08/2014



% Note: Problem: cells may partially share a facet, i.e. touch, but not have the
% vertices in common. This is because tht too big cells are artificially
% cut off.
% ** This is now treated in Makevoronoi.


if nargin<3 %Compute Delaunay triangulation, if not provided
    deltri = delaunay(coordlist(:,1), coordlist(:,2));
end

% if nargin>2, s2=size(faces,2); else s2=10; end

% adjacentlist = zeros(size(coordlist,1), s2);
% index = ones(size(coordlist,1), 1);


% neighbors = [];

% for nf = 1:numel(nface)
% %     s = sum(deltri==nface(nf), 2) > 0;    % lines of deltri that contain nf (logical index)
%     s = any(deltri==nface(nf), 2);
%     allv = deltri(s,:);                   % content of these lines
%     neighbors = [neighbors; setdiff( unique(allv(:)), nface(nf) ) ]; % -"-, unique and without nf
% end   
% 
% N = false(size(deltri,1),1);

% for nf = 1:numel(nface)
%     N = N | any(deltri==nface(nf), 2);
% end

N = any(ismember(deltri,nface),2);    % lines of deltri that contain a member of nface (logical index)

allv = deltri(N,:);                   % content of these lines
if numel(nface)~=1
    neighbors = setdiff( unique(allv(:)), nface ); 
else % do same without setdiff (faster)
    allv = unique(allv(:));
    neighbors = allv(allv~=nface);
end

% neighbors = unique(neighbors);

