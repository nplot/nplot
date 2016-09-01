function [data, errorstate] = calclinearcombination(data, calcstring, varargin)

% Calculate a linear combination of pal-states for a polarized data set
% calcstring like "2*pal1+pal3-3*pal4" etc.
% varargin can contain: 'output' (false,true) and 'dist' (for subtractdata)

% P. Steffens, 11/2009


errorstate = false;
nooutput = readinput('output',varargin);
if isempty(nooutput), nooutput = false; end
gridstep = readinput('dist',varargin);
if isempty(gridstep), try %#ok<ALIGN>
    stdbindist = getoption('stdbindist'); 
    gridstep = stdbindist.(data.coordtype); 
    catch fprintf('Error: Need a ''dist'' parameter in ''calclinearcombination''.\n'); errorstate = true; return; end
end

    
data.legend = calcstring;
dataname = [calcstring, ': ', data.dataname];

[st,en] = regexp(calcstring, 'PAL\d+'); %find all PALx entries in calcstring
for i=numel(st):-1:1
    % Convert calcstring in matlab readable string (PALx --> PAL(:,x))
    calcstring = [ calcstring(1:st(i)+2) '(:,' calcstring(st(i)+3:en(i)) ')' calcstring(en(i)+1:end) ];
end
% Find coefficients for each PAL-column
PAL = eye(max(data.pallist)); %#ok<NASGU>
pcoeff = [];
try
    eval(['pcoeff = ' calcstring ';']);
    if ~nooutput, fprintf(['** Calculating pal-combination. (Works only for linear combinations.) Coefficients are ' num2str(pcoeff(:)','%g  ') '\n']); end
catch
    em = lasterror;
    fprintf(['Error during calculation: ''' em.message '''  --> Check input! \n']); data = []; errrostate = true; return; 
end
% Use 'scaledata' and 'subtractdata' to do the linear combination
cind = find(pcoeff);
comb = scaledata( extractsubset(data,'pallist',cind(1)), pcoeff(cind(1)) );
if ~nooutput, fprintf('%s \n', comb.legend); end
for i = 2:numel(cind)
    nextdat = extractsubset(data,'pallist',cind(i));
    comb = subtractdata( comb, scaledata(nextdat , -pcoeff(cind(i)) ), 'nearest', abs(gridstep) );
    if ~nooutput, fprintf('%s \n', nextdat.legend); end
end
data = comb;
data.polarized = false;
if isfield(data,'taglist'), data = rmfield(data,'taglist'); end
if isfield(data,'dataname'), data.dataname = dataname; end

