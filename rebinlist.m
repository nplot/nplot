function [newx,newy,newdy] = rebinlist(x,y,dy,gridx)

% function [newx,newy,newdy] = rebinlist(x,y,dy,gridx)


gridx = sort(gridx(:));
x=x(:); y=y(:); dy=dy(:);
binlist = bintogrid(x, gridx, max(abs(gridx(2:end)-gridx(1:end-1))));

[ubl2, m , ubl2num] = unique(binlist(:,2));
% newx = zeros(numel(ubl2),size(coordlist,2));
% newy = zeros(numel(ubl2), 2);

ind = 1;
for i=ubl2'         %i is the "bin"-number
    lines = binlist(logical(binlist(:,2)==i),1);  %"lines" are all datapoints that belong to this bin
    [newy(ind), newdy(ind)] = weightedmean( y(lines), dy(lines) ); %(standard weights 1/sig^2)
    newx (ind) = gridx(i);
    ind = ind + 1;
end
