function newlist = smoothdata(datalist,range)

% function newlist = smoothdata(datalist,range)
% Smooth data by averaging over a given range (in n dim.)
% (Weighted average of points within range)

% P. Steffens, 06/2008

newlist = [];
ndims = size(datalist.coordlist,2);


if numel(range)~=ndims, fprintf('Error: (smoothlist): numel(range) does not correspond to number of dimensions.\n');  return; end
newlist = datalist;
normval = getoption('normval');

for i=1:size(datalist.coordlist,1)
    ind = ( abs( datalist.coordlist(:,1) - datalist.coordlist(i,1) ) <= range(1) );
    for nd = 2:ndims
        ind = ind & (abs( datalist.coordlist(:,nd) - datalist.coordlist(i,nd) ) <= range(nd) );
    end
    if datalist.raw ==1 % For counting data: just add, normalize and obtain error as sqrt
        counts = datalist.valuelist(ind,1) .* datalist.monitorlist(ind,2) / normval;  % These are the real original counts!
        newlist.monitorlist(i,:) = sum( datalist.monitorlist(ind,:), 1 );
        newlist.valuelist(i,1) =       sum(counts) / newlist.monitorlist(i,2) * normval;
        newlist.valuelist(i,2) = sqrt(sum(counts)) / newlist.monitorlist(i,2) * normval;
    else % Otherwise weighted average
        [newlist.valuelist(i,1), newlist.valuelist(i,2)] = weightedmean ( datalist.valuelist(ind,1), datalist.valuelist(ind,2) ); 
    end
end       

% % % Kann man das ganze irgendwie optimieren ?
% % % z.B. wenn man stattdessen einfach die Nachbarn summiert ?
% % val = datalist.valuelist(:,1);
% % for i=1:3       % besser jede Kante nur einzeln zählen!
% %     intsum(:,i) = sum(val(datalist.delaunaytri(:,setdiff(1:3,i))),2); % intensity sum of the two connectected points
% % end


newlist.raw = false;