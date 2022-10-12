function [val,rest] = readinput(name, arglist, varargin)

% function [val,rest] = readinput(name, arglist,varargin)
%
% Analyse the cell array arglist if the parameter (name) is existent and
% return its value (the following cell, or what comes after 'name' in a single string. 
% val is cell array if several occurences of 'name' are found.
% If vargargin = 'first' or 'last', return only first or last occurrence.
% rest is the part of arglist which does not contain 'name'.

% P. Steffens, 05/2018

val = [];
rest = [];

name = strtrim(name);
if ~iscell(arglist)
    m = arglist;
    for i = 1:length(m)
        arglist{i}=m(i); 
    end
    clear m; 
end % if not, make cell array


% Check for name/value pairs in single string, and convert to real pairs
ni = 0;
while ni < length(arglist)
    ni = ni+1;
    if ~ischar(arglist{ni}), continue; end
    argni = strtrim(arglist{ni});
    if length(argni)>length(name) && strcmpi(argni(1:length(name)+1), [name ' '])
        valstring = strtrim(argni(length(name)+1:end));
        arglist = {arglist{1:ni-1}, name, valstring, arglist{ni+1:end}};
    end
end


% rest = arglist;

% Go through list


vi = find(strcmpi(arglist,name));

if ~isempty(vi) && ~isempty(varargin) && any(strcmpi(varargin,'last')), vi = vi(end); end
if ~isempty(vi) && ~isempty(varargin) && any(strcmpi(varargin,'first')), vi = vi(1); end

for p=1:numel(vi)    % Allow for multiple occurences    
    if vi(p) >= length(arglist), continue; end
    valnew = arglist{vi(p)+1};

    if ischar(valnew) % try to convert to numeric
        valnew=strtrim(valnew);
        valn = str2num(arglist{vi(p)+1}); %#ok<ST2NM> 
        if ~isempty(valn) && (isnumeric(valn) || isobject(valn)) && ~isa(valn,'function_handle') && ~isa(valn,'iFunc'), valnew=valn; end
    end
    if ischar(valnew) && valnew(1)=='{' && valnew(end)=='}'
        % try to convert to cell array
        valnew = eval(valnew);
    end
    
    if isempty(val)
        val = valnew;     % first entry
    elseif ~iscell(val)
        vallist1{1} = val; % second entry, make cell array
        val = vallist1; 
        val{2} = valnew;
    else              % further entries
        val{length(val)+1} = valnew;
    end
end
    
restind = setdiff(1:length(arglist),[vi, vi+1]);
if restind
    rest={}; 
    for i=1:numel(restind), rest{i} = arglist{restind(i)}; end
end

