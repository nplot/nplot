function sums = sumsteps(scan,mode)

% Returns sum of all channels for each scan step (sum over channels)
% Use plot(sumsteps(..)) to plot
% scan can be filename or scan structure
% mode: 'original' (standard), 'monitor' normalized, 'relative' (mean=1)

% P. Steffens, 04/2008

%%
% If argument is filename, load scan first
if ~isstruct(scan)
    scan = tasread(scan);
end

sums = sum(scan.MULTI, 2);  % do simple summation...

if nargin < 2,  mode = 'original'; end

switch upper(mode)
    case 'ORIGINAL' ;
    case 'MONITOR'
        [normalizeto,normval] = getoption('normalizeto','normval');
        norm_measured = getvar(scan,normalizeto);
        multidat  = scan.MULTI ./ repmat(norm_measured,[1,31]) * normval ;
        sums = sum(multidat, 1);
    case 'RELATIVE'
        sums = sums / mean(sums);
    otherwise, fprintf('Argument not recognized: sue ''original'', ''monitor'' or ''relative''\n');
end
