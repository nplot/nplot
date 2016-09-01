function [paldeflist, assignpal] = analyzepal(scan, paldeflist)

% function [paldeflist, assignpal] = analyzepal(scan, paldeflist)
% Analyze the POLAN section of the scan
% reteurn an extended paldeflist and the assignment table of PAL-values to
% this list

% P. Steffens, 07/2015


npollines = length(scan.POLAN);
paldef = [];
npal = 1;
currentvar = [];
currentval = [];
for p = 1:npollines
    if ~isempty(regexp(upper(scan.POLAN{p}),'^\s*DR', 'once' ))
        [st,en] = regexp(upper(scan.POLAN{p}),'(?<=DR.*\s+)\S+');  % all arguments of drive command
        for a = 1:numel(st)
            [s,e] = regexp(scan.POLAN{p}(st(a):en(a)),'[a-zA-Z]\w+'); % match variable name
            if ~isempty(s)
                if ~isempty(currentval) 
                    paldef.(currentvar) = currentval;
                end
                currentvar = upper(scan.POLAN{p}((st(a)-1+s):(st(a)-1+e))); 
                currentval = [];
            end
            [s,e] = regexp(scan.POLAN{p}(st(a):en(a)),'^-?(\d*\.?\d+|\d+\.?\d*)'); % match value
            if ~isempty(s) && ~strcmpi(currentvar, 'a3b')
                currentval = [currentval, str2double(scan.POLAN{p}((st(a)-1+s):(st(a)-1+e)))]; 
            end
        end
    end
    if ~isempty(regexp(upper(scan.POLAN{p}),'ON', 'once' ))
        [st,en] = regexp(upper(scan.POLAN{p}),'(?<=ON\s+(\S+\s+)*)\S+|(\S+)(?=\s+ON)'); % all arguments of ON-command
        for a = 1:numel(st)
            if ~isempty(currentval), paldef.(currentvar) = currentval; end
            currentval = [];
            paldef.(upper(scan.POLAN{p}(st(a):en(a)))) = 'on ';
        end
    end
    if ~isempty(regexp(upper(scan.POLAN{p}),'OF', 'once' ))
        [st,en] = regexp(upper(scan.POLAN{p}),'(?<=(OF|OFF)\s+(\S+\s+)*)\S+|(\S+)(?=\s+OFF)'); % all arguments of OFF-command
        for a = 1:numel(st)
            if ~isempty(currentval), paldef.(currentvar) = currentval; end
            currentval = [];
            paldef.(upper(scan.POLAN{p}(st(a):en(a)))) = 'off';
        end
    end
    if ~isempty(regexp(upper(scan.POLAN{p}),'FCU', 'once' ))
        [st,en] = regexp(upper(scan.POLAN{p}),'UP|DOWN'); % all arguments of OFF-command
        paldef.fcu = upper(scan.POLAN{p}(st(1):en(1)));
    end
    if ~isempty(regexp(upper(scan.POLAN{p}),'CO', 'once' ))
        if ~isempty(currentval), paldef.(currentvar) = currentval; end
        % Check if this paldef already exists in paldeflist...
        errorstate = false;

        assignpal(npal) = 0;
        for pd = 1:length(paldeflist)
            fn = fieldnames(paldeflist{pd});
            if length(fn) ~= length(fieldnames(paldef)), errorstate = true; end
            fitstate = true;
            for fni = 1:length(fn)
                if isfield(paldef, fn{fni})
                    fitstate = fitstate && length(paldef.(fn{fni})) == length(paldeflist{pd}.(fn{fni})) && all(paldef.(fn{fni}) == paldeflist{pd}.(fn{fni}));
                else
                    errorstate = true;
                end
            end
            if fitstate, assignpal(npal) = pd; end
        end
        if errorstate, fprintf('Error: pal-Definition %d in file %d is not compatible with those from previous files!\n', npal, str2num(scan.FILE)); avgdata = []; return; end 
        % if not found, append paldeflist        
        if assignpal(npal) == 0,
            paldeflist = {paldeflist{:}, paldef};
            assignpal(npal) = length(paldeflist);

        end
        npal = npal + 1;
    end            
end
assignpal = assignpal(:);