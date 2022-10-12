function [cutvertices, cutfaces, numbers] = createintersection (vertexlist, cells, pos, celldim)

% Obtain the intersection of the patch defined (vertexlist,cells) with the
% hypersurface f(x)=0, where f(x) for each point is given in pos.
% As it is not evident from spacedim what effective dimensionality have the
% cells, the argument celldim is necessary (=2 for "plottable" surfaces in N-dim space).
% (For celldim=2 the rows in cells must have the right order.)
%
% Relation to cutpatch: give intersection instead of two halfs
% Relation to createmesh: result not necessarily planar, dimensionality is always reduced by one. createmesh gives ordered output of 2d cuts.

% P. Steffens, 04/2009


% %%
% function orderpoints
% % bring the point set cutpoints in the right order
% % cutpoints lie on a plane (or nearly in more general case). Find 2D coordinates by projection on plane
% [normal,c] = fithyperplane ( cutpoints );
% 
% end
% 
% %%

spacedim = size(vertexlist,2);
nvert    = size(vertexlist,1);

% Determine which cells cross the border (code like in cutpatch)
pos = pos(:);
vside = sign(pos);
vsidepos = [vside>0; true];
vsideneg = [vside<0; true];
testfaces = cells;
testfaces(isnan(testfaces)) = nvert + 1;

if size(cells,1) == 1  %strange, here "all" behaves in a bizarre way...
    allpos = all(vsidepos(testfaces));
    allneg = all(vsideneg(testfaces));
else
    allpos = all(vsidepos(testfaces),2);
    allneg = all(vsideneg(testfaces),2);
end

indextwoside = find(~allpos & ~allneg); % These are the faces to be cut


cutvertices = zeros(0,spacedim);
cutfaces = [];
cutvertexcount = 0;
cutfacecount = 1;

for ncell = indextwoside'
    
    fa = cells(ncell, isfinite(cells(ncell,:)));    % vertices defining the current cell
    if celldim == 2 
        K = [ (1:numel(fa))', [2:numel(fa),1]' ];   % in 2D there is an order, otherwise find K in polyhypercut
    else K=[]; 
    end
    
    cutpoints = polyhypercut( vertexlist(fa,:),[], 0, K, pos(fa'));
    
    ncp = size(cutpoints,1);
%     if celldim==3 && nargout>1, orderpoints; end  % The cuts are (2d)polygons in space. Find the right order of the points
    cutvertices = [cutvertices; cutpoints];
    cutfaces(cutfacecount,1:ncp) = cutvertexcount + (1:ncp);
    numbers(cutfacecount) = ncell;
    cutvertexcount = cutvertexcount + ncp;
    cutfacecount = cutfacecount + 1;
    
end

cutfaces(cutfaces==0) = NaN;


[cutvertices, cutfaces] = optimizepatch(cutvertices, cutfaces);
    
    
end
    
    
    
    