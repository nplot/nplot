function [message, rest, multiple] = checkoptions(optionstring, knownoptions, knownswitches, notwarnmultiple)

% Check if the cell array optionstring (usually the input to the calling
% function) contains only known pairs options (name + value) or switches.
% Known options and switches are defined by the cell arrays knownoptions
% and knownswitches. 
% [message, rest, multiple] are a message to display on screen, and two 
% cell arrays, containing the names of not interpretable parts, and of parts
% which occur more than once. 

% P. Steffens 07/2013

rest = optionstring; 
multiple = {};
reststring = [];
multstring = [];
message = [];

% Check if everything can be interpreted
for i=1:length(knownoptions)
    [o,rest] = readinput(knownoptions{i},rest);
%     if iscell(o) && length(o)>1, multiple = {multiple{:}, knownoptions{i}}; end
end
if ~isempty(rest)
    for i=1:length(rest), if ~strcmpi(rest{i},knownswitches), reststring = [reststring, rest{i}, ' ']; end; end
    if ~isempty(reststring), message = ['Unknown options: ', reststring, '\n']; end
end
% Check for multiple occurences
for i=setdiff(1:(length(optionstring)-1),notwarnmultiple)
    if ~ischar(optionstring{i}), continue; end
    optionstring{i}=strtrim(optionstring{i});
    name = optionstring{i}(1:min([end,(regexp(optionstring{i},' ','once')-1)]));
    if any(strcmpi(name,{knownoptions{:},knownswitches{:}}))
        for j=setdiff((i+1):length(optionstring),notwarnmultiple)
            if ischar(optionstring{j}) && strncmpi(name,strtrim(optionstring{j}),length(name)) && ~any(strcmpi(multiple,name)), multiple = {multiple{:},name}; end
        end
    end
end    
if ~isempty(multiple)
    for i=1:length(multiple), multstring = [multstring, multiple{i}, ' '];  end
    if ~isempty(multstring), message = [message, 'Options occuring several times: ', multstring, '\n']; end
end
