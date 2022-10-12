function output = fitfunctions

% call 'fitfunctions' to obtain a list of available functions

files = what('Functions');
if isempty(files), fprintf('Folder "Functions" not found.\n'); return; end
files = files.m;

output = {};

for nf = 1:length(files)
    fname = files{nf};
    [st,en] = regexp(fname,'\w+(?=\.m$)');
    fname = fname(st:en);
    if ~isempty(fname)
        try
            eval([' [val, paramnames, paramnum, description] = ' fname '([],[]);']);
            output = [output; {fname, description}];
        catch
        end
    end
end

% Print on screen if no output variable
if nargout == 0
    fprintf('----------------------------\nList of available functions:\n----------------------------\n');
    for li = 1:size(output,1)
            fprintf('%15s.m: %s \n', output{li,1}, output{li,2});
    end
    fprintf('Type ''showfunction [functionname]'' for help on particular function.\n');
    fprintf('Use file ''template.m'' to implement a new function.\n');
    clear output 
end
