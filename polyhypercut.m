function [cutpoints, poly1, poly2, npoly1, npoly2] = polyhypercut(pointset, normal, const, K, pos)

% Cut an N-dim convex polyeder in two parts
% normal * x = const  defines an (N-1 dim.) hyperplane in N-dim space
% pointset: a set of points in N-dim space that define a convex polyeder
% K :       if given, save the (costly) call to convhulln
% pos :     if given, do not calculate normal*x for each point, but use pos instead. (normal is ignored) 
%           This allows to use other surfaces than hyperplanes to cut (i.e. nonlinear).
%           In this case however, the instersection point is obtained via linear interpolation. 
% 
% cutpoints :    vertices of the intersection with the hyperplane
% poly1,2 :      point sets describing the two parts  (poly1 are those with x*n < c)
% npoly1,npoly2: indices into [pointset;cutpoints]
%
% The resulting pointlists are not ordered. In 2D, one can obtain ordered
% poly1 and poly 2 by giving an ordered pointset and K=[1,2; 2,3; 3,4; ...; n,1].
% In this special case it also works for non-convex polygons. 

% Attention with flat facets (e.g. flat 2D polygons in 3d): result contains
% more points than necessary (but is not false)
% Attention 2: point list is not necessarily unique!

% P. Steffens, 04/2009


normal = normal(:);
npoints = size(pointset,1);
spacedims = size(pointset,2);

% For all points in pointlist, determine position with respect to
% hyperplane, if argument pos is not given
if nargin < 5 || isempty(pos)
    pos = zeros(npoints,1);
    for d = 1:spacedims
        pos = pos + pointset(:,d) * normal(d);   % ("normal" is used only here!)
    end
end

% If Polygon entirely lies on one side (within precision), do not do the rest and
% return no cutpoints
if all(pos > const-1E-10) % if min(pos) >= const - 1E-10
    poly1 = [];
    poly2 = pointset;
    cutpoints = [];
    return;
elseif all(pos < const + 1E-10) % max(pos) <= const + 1E-10
    poly1 = pointset;
    poly2 = [];
    cutpoints = [];
    return;
end

side1 = pos <= const;% + 1E-10;
side2 = pos >= const;% - 1E-10;

if nargin<4 || isempty(K)
    K = convhulln(pointset);
    % Man könnte versuchen, ohne auszukommen, indem man ALLE Paare
    % betrachtet (weiterer Aufruf von convhulln ist für das Ergebnis
    % sowieso nötig)
end

% find all edges
% edges = zeros(size(K,1) * size(K,2)*(size(K,2)-1)/2,  2);
% linind = zeros(1, size(edges,1));
% nfacs = zeros(1, size(edges,1));
% ncuts=0;

% Find all edges (pairs of points) that cut the plane

K1 = side1(K);

e1 = zeros(0,1);
e2 = zeros(0,1);

celldims = size(K,2);
for i=1:(celldims-1)
    for j=(i+1):celldims
        crossedges = (K1(:,i) ~= K1(:,j)); % points on different sides
        e1 = [e1; K(crossedges,i)]; % indices of these points
        e2 = [e2; K(crossedges,j)];
    end
end
% create a number edgeid that represents each pair (symmetrically)
a = e1 + e2;
b = abs(e1 - e2);
edgeid = (a+b) .* (a+b-1)/2 + a;

%[b,m,n] = unique(edgeid);

% Find identical edgeid to avoid calculating the same edge several times   ** (can avoid this in low Dim??)
[sortededgeid,order] = sort(edgeid);
goodind = [true; sortededgeid(2:end)~=sortededgeid(1:(end-1))];

goodind(order) = goodind;
e1 = e1(goodind);
e2 = e2(goodind);
ncuts = sum(goodind);
% now e1 and e2 are the endpoints of the ncuts different edges     

% Find intersection of edges with hyperplane
% Calculate lambda = (c-p1*n)/((p2-p1)*n) and cutpoints = p1 + lambda*(p2-p1) 

% p1n = zeros(ncuts,1);
% denom=zeros(ncuts,1);
% for d = 1:ndims
%     p1 = pointset(e1,d);
%     p2 = pointset(e2,d);
%     p1n   = p1n + p1 * normal(d);
%     denom = denom + (p2 - p1) * normal(d);
% end
% lambda = (const - p1n) ./ denom;

lambda = (const-pos(e1))./(pos(e2)-pos(e1)); % This is better
cutpoints = pointset(e1,:);
p2mp1 = pointset(e2,:) - pointset(e1,:);
for d = 1:spacedims
    cutpoints(:,d) = cutpoints(:,d) + lambda .* p2mp1(:,d);
end
   

% 
% %for fac = 1:size(K,1)
% facets = find((sum(K1,2)<3) & (sum(K2,2)<3));
% for fac = 1:numel(facets)
%     % enumerate vertices of facet pairwise
%     for i=1:(ndims-1)
%         for j=(i+1):ndims
%             e1 = K(facets(fac),i);
%             e2 = K(facets(fac),j);
%             if side1(e1)==side1(e2)
%                 continue;
%                 % Take only those edges that connect points on different sides of the plane
%             end
%             ncuts=ncuts+1;
%             edges(ncuts,:) = [e1 , e2];
%             nfacs(ncuts)=facets(fac);
%             % Construct a single number linind that represents the pair (e1,e2):
%             a = e1+e2;
%             b = abs(e1-e2);
%             linind(ncuts) = (a+b)*(a+b-1)/2 + a;
%         end
%     end
% end 
% edgenumber = nfacs(1:ncuts);
% [linind,order] = sort(linind(1:ncuts));
% edges = edges(order,:);
% cutpoints = zeros(ncuts,ndims);
% cutused=false(ncuts,1);
% for e=1:ncuts
%     if (e>1 && linind(e)==linind(e-1)), continue; end
%     cutused(e)=true;
%     p1 = pointset(edges(e,1),:);
%     p2 = pointset(edges(e,2),:);
%     % Find intersection of edge with hyperplane
%     lambda = (const - p1*normal) / ((p2-p1)*normal);
%     cutpoints(e,:) =  p1 + lambda*(p2-p1);
% end
% cutused(order)=cutused;
% cutpoints(order,:) = cutpoints;
% cutpoints = cutpoints(cutused,:);


if nargout < 2, return; end

if celldims==2 && nargin > 3 && ~isempty(K) % in 2D, try to preserve order
    edgenumber = find(crossedges) + (1:(ncuts))';  % cutpoints at these positions
    % Create long pointlist allpoints with original and cutpoints and
    % adjust logical indices from above...
    orig = true(npoints+ncuts,1);
    side1all = orig;
    side2all = orig;
    
    orig(edgenumber) = false;
    side1all(orig) = side1;
    side2all(orig) = side2;
    
    allpoints=zeros(npoints+ncuts,spacedims);
    allpoints(orig,:) = pointset;
    allpoints(edgenumber,:) = cutpoints;
    
    poly1 = allpoints(side1all,:);
    poly2 = allpoints(side2all,:);
    
    if nargout>3
        pointnumberlist(orig) = 1:npoints;
        pointnumberlist(edgenumber) = (npoints+1):(npoints+ncuts);
        npoly1 = pointnumberlist(side1all);
        npoly2 = pointnumberlist(side2all);
    end
    
else  % otherwise, just return the points in simplest order
    poly1 = [pointset(side1, :); cutpoints]; %if size(poly1,1) <= ndims, poly1 = zeros(0,ndims); end
    poly2 = [pointset(side2, :); cutpoints]; %if size(poly2,1) <= ndims, poly2 = zeros(0,ndims); end   % if <=ndims points, they do not define an n-dim. Polyeder (just cut on the edge)
    if nargout>3
        npoly1 = [find(side1); ((npoints+1):(npoints+ncuts))'];
        npoly2 = [find(side2); ((npoints+1):(npoints+ncuts))'];    
    end
end

 
    
    
    