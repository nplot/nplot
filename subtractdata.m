function erg = subtractdata (list, minuslist, mode, opt)

% Subtract 'minuslist' from 'list'
%
% P. Steffens, 10/2008


erg = [];

% if isnumeric(minuslist) || numel(minuslist)==1
    

% First, check consistency of lists (coords, etc...) !!! 

%Check if  lists are of the same type
if ~strcmpi(list.type, minuslist.type), 
    fprintf('Error: tried to subtract data lists of different type! This is impossible. \n'); return; end
if ~strcmpi(list.coordtype, minuslist.coordtype), 
    fprintf('Error: tried to subtract data lists with different types of corrdinates! This is impossible. \n');  return;  end

%Check consistency of constants. However, do not alter constants of 'list', only give warnings
[eqset, eqvalues, deviate, notall] = checkconstants({list, minuslist});
%Warnings for missing values
for j=1:length(notall), fprintf('Warning: Value for "%s" does not exist in both lists! Continue... \n', notall{j}); end
%Warning for too large deviations
if ~isempty(deviate)
    errtxt=[]; for j=1:length(deviate), errtxt = [errtxt, deviate{j}, ' ']; end
    fprintf('Warning: too large deviation in %s! Continue... (Please check!)\n   (You may increase the maximum acceptance in the options file.) \n', errtxt);
end


if nargin < 3, mode = 'nearest'; end


switch upper(mode)
    
    case 'NEAREST'    
        if nargin>3, bindist = opt; else stdbindist = getoption('stdbindist'); bindist = stdbindist.(list.coordtype); end
        % Perform a binning to 'minuslist' 
        binnedpoints = bintogrid(list.coordlist, minuslist.coordlist, bindist);
        % ...and continue just with the points that are not too far away from any binpoint
        list = sublist(list, binnedpoints(:,1));
        minusval = minuslist.valuelist(binnedpoints(:,2), 1);
        minuserr = minuslist.valuelist(binnedpoints(:,2), 2);
        
    case 'INTERPOLATE'
        if nargin>3, dist = opt; else stdbindist = getoption('stdbindist'); dist = 5 * stdbindist.(list.coordtype); end % ** maybe [inf,inf]?
        % Interpolate minuslist to obtain values to be subtracted at the correct coordinates
        [minusval, minuserr] = linearinterpolation(minuslist, [list.coordlist(:,1),list.coordlist(:,2)], dist);
        % ** Attention, interpolation may return NaN's at the edges (coord within bin,
        % but outside convex hull of minuslist). Treat these separately: **
        %   ???
        
    case 'RANGE'
        if nargin>3, dist = opt; else dist = [inf,inf]; end 
        % Obtain subtracted value by averaging over a given range
        [minusval, minuserr] = lookaround(minuslist,[list.coordlist(:,1),list.coordlist(:,2)],dist);

end


% Now, perform the actual subtraction
list.valuelist(:,1) = list.valuelist(:,1) - minusval;
list.valuelist(:,2) = sqrt(list.valuelist(:,2).^2 + minuserr.^2);

% Delete NaN's
finites = ~isnan(list.valuelist(:,1));
list = sublist(list, find(finites));

list.raw = false;
if isfield(list, 'monitorlist'), list = rmfield(list,'monitorlist'); end
if isfield(list, 'taglist'), list = rmfield(list,'taglist'); end

if isfield(list, 'dataname') && isfield(minuslist, 'dataname'), list.dataname = ['(' list.dataname ') - (' minuslist.dataname ')']; end

erg = list;
