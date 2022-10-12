function newdat = rebindata(data, stepx, stepy)

% function newdat = rebindata(data, stepx, stepy)
%
% 2-dim. rebinning of data points. Give ALL x and y values in stepx,stepy.
% (New grid is meshgrid(stepx, stepy))

[xx,yy] = meshgrid(stepx, stepy);
bins = bintogrid(data.coordlist, [xx(:),yy(:)],[.1,.1]); % ** vernünftige maxdist!!!!
newdat = sublist(data, bins(:,1));
newdat.coordlist = [xx(bins(:,2)), yy(bins(:,2))];
newdat = cmbavg(newdat);
