function steps = scansteps(command)

% extracts steps of scan from Command

% PS 06/23

steps = [];

vars = regexpmatch(upper(command),'(?<=\s+D)\w+(?=\s+\-?(\d*\.?\d+|\d+\.?\d*))');           %     \-?(\d*\.?\d+|\d+\.?\d*) is a number
% To recognize scanned variables, look at "D.." parts of COMND (the given steps)

vals = regexpmatch(upper(command),'(?<=\s+D\w+\s+)\-?(\d*\.?\d+|\d+\.?\d*)');
% corresponding values


% check for special case of variables starting with "d"
ind = true(1,length(vars));
for i=1:length(vars)
    if vars{i}(1)=='D' %special case for variable names starting with "d"
        ind(strcmpi(vars(1:i),vars{i}(2:end)))=false;
    end
end
vars = vars(ind);
vals = vals(ind);


for i=1:length(vars)
    steps.(vars{i}) = str2double(vals{i});
    if strcmp(vars{i},'QH')
        dqk = regexpmatch(upper(command), '(?<=\s+DQH\s+(\-?(\d*\.?\d+|\d+\.?\d*)\s+){1})\-?(\d*\.?\d+|\d+\.?\d*)');
        if isempty(dqk), steps.QK = 0; else, steps.QK = str2double(dqk{1}); end
        dql = regexpmatch(upper(command), '(?<=\s+DQH\s+(\-?(\d*\.?\d+|\d+\.?\d*)\s+){2})\-?(\d*\.?\d+|\d+\.?\d*)');
        if isempty(dql), steps.QL = 0; else, steps.QL = str2double(dql{1}); end
        den = regexpmatch(upper(command), '(?<=\s+DQH\s+(\-?(\d*\.?\d+|\d+\.?\d*)\s+){3})\-?(\d*\.?\d+|\d+\.?\d*)');
        if isempty(den), steps.EN = 0; else, steps.EN = str2double(den{1}); end      
    end
end