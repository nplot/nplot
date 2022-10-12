function [grid,assign]=makegrid(coords, centers, steps)

% Determine a set of grid points that is best suited to bin the data to
% (for instance contains only points to which points in (x,y) will likely
% be assigned).
%
% Works in n dimensions.

% centers may be empty (--> calculate it) or multiple

% To be extended and improved...


% P. Steffens 08/2008


npoints = size(coords,1);
ndims = size(coords,2);

% if centers not given, try to find a center such that the binning leaves as
% many points as possible unshifted:
if isempty(centers)
   testcenter = median(coords, 1);
   for d=1:ndims
       mods = mod(coords(:,d) - testcenter(d), steps(d));
       mods = round(mods*10/steps(d)) *steps(d)/10;
       cc=0;
       for i=unique(mods)'
           if sum(mods==i)>cc, cc=sum(mods==i); sh=i; end
       end
       centers(1,d) = testcenter(d) + sh;
   end       
end

ncent = size(centers,1);

if size(steps,1)==1
    steps = repmat(steps,[ncent,1]);
end

dist = zeros(npoints, size(centers,1));
for nc = 1:ncent
    for d = 1:ndims
        dist(:,nc) = dist(:,nc) + ( coords(:,d) - centers(nc,d) ).^2;
    end
end
[m,ic] = min(dist,[],2);
stepdist = zeros(npoints,ndims);
for d = 1:ndims
    stepdist(:,d) = round( (coords(:,d) - centers(ic,d)) ./ steps(ic,d) );
end
list = [ic, stepdist];


[shortlist, n, assign] = unique(list, 'rows');

grid = zeros(size(shortlist,1),ndims);

for i=1:size(shortlist,1)
    grid(i,:) = centers(shortlist(i,1),:) + steps(shortlist(i,1),:) .* shortlist(i,2:end);
end