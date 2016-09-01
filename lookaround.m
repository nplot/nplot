function [zi,errzi] = lookaround(datalist,coord,range)

% Obtain values at coord by averaging over a given range of datalist.
% (Weighted average of points within range)
% Considers the n Dims of coord, evtl. additional dims of datalist are ignored
% (Corresponds to smoothdata.)

% P. Steffens, 06/2008

zi = []; errzi = [];
ndims = size(coord,2);
if numel(range)~=ndims, fprintf('Error: (lookaraound): numel(range) does not correspond to number of dimensions.\n');  return; end
normval = getoption('normval');
zi    = NaN(size(coord,1),1);
errzi = NaN(size(coord,1),1);

for i=1:size(coord,1)
    ind = ( abs( datalist.coordlist(:,1) - coord(i,1) ) <= range(1) );
    for nd = 2:ndims
        ind = ind & (abs( datalist.coordlist(:,nd) - coord(i,nd) ) <= range(nd) );
    end
    if datalist.raw ==1 % For counting data: just add, normalize and obtain error as sqrt
        counts = datalist.valuelist(ind,1) .* datalist.monitorlist(ind,2) / normval;  % These are the real original counts!
        newlist.monitorlist(i,:) = sum( datalist.monitorlist(ind,:), 1 );
        zi(i)    =       sum(counts)  / newlist.monitorlist(i,2) * normval;
        errzi(i) =  sqrt(sum(counts)) / newlist.monitorlist(i,2) * normval;
    else % Otherwise weighted average
        [zi(i), errzi(i)] = weightedmean ( datalist.valuelist(ind,1), datalist.valuelist(ind,2) ); 
    end
end       
        
newlist.raw = false;