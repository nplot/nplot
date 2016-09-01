function erg = bintogrid(datapoints,gridpoints,maxdist)

% Find grid points that are closest to given data points ("binning")
%
% Works in arbitrary number of dimensions (like e.g. qx,qy,qz,en,temp,...)
% and with arbitrary grids (equidistant or not...).
% Points having larger than given maximum distance are ignored.
% No averaging is done in this routine.
%
% maxdist:  vector of maximum acceptable distances to nearest grid point
%           along the respective dimension
% datapoints, gridpoints:
%           coordinate list of data and grid points; 
%           rows are the points, columns the coordinates, i.e.
%           size(..)=[#points,#dims]
%
% output:   list of data points (first column) with the number of respective
%           grid point (second column)
%
% P. Steffens 06/2012


erg = [];
np=size(datapoints,1); %Number of data points
ng=size(gridpoints,1); %Number of grid points
dims = size(datapoints,2);
if size(gridpoints,2) ~= dims, fprintf('Error: Dimension mismatch between data and grid in "bintogrid.m"\n'); return; end
if numel(maxdist)==1, maxdist = repmat(maxdist,1,dims); end

% try     % version without loop, but very memory-consuming...
%     
%     %create matrices in convenient format...
%     G = repmat(gridpoints,[1,1,np]);
%     P = repmat(datapoints,[1,1,ng]);
%     P = permute(P,[3,2,1]);
% 
%     %Calculate the distances...
%     %Differences in all coordinates
%     dif  = G-P;
%     % Test if any is larger than limit (--> set to NaN)
%     dif(logical(abs(dif)>repmat(maxdist(:)',[ng,1,np])))=NaN;
%     %Distances by pairs
%     dist = sum(dif.^2,2);
%     dist = squeeze(dist);
% 
%     %For each datapoint: nearest distance (m) to grid point (ind)
%     [mindist,ind] = min(dist,[],1);
% 
% catch   % do it in a loop
%     
%     clear G P dif dist
    mindist = zeros(1,np);
    ind  = zeros(1,np);
    for i=1:np
       dgrpt = zeros(ng,1);
       for d=1:dims
           dgrpt = dgrpt + (gridpoints(:,d) - datapoints(i,d)).^2;
       end
       [mindist(i), ind(i)] = min(dgrpt);
       if any( abs(gridpoints(ind(i),:)-datapoints(i,:)) > maxdist ), mindist(i)=nan; end
    end
    
% end
    
%Discard too distant data points (mindist==NaN)
pointlist = find(isfinite(mindist));
ind = ind(pointlist);

erg = [pointlist', ind'];
