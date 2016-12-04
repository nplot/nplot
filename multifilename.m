function filelist = multifilename(file)
% Returns cell array of full path filenames
% P. Steffens, 4/2015

function res = ssplit(str) % Replaces Matlab strsplit (introduced in MATLAB 8.1 (R2013a))
    res={};
    for cmind = fliplr([0,find(str==',')])
        res = [str(cmind+1:end), res];  %#ok<*AGROW>
        str = str(1:cmind-1);
    end
end

filelist = {};

if ~isempty(strfind(file,'[')) % if there are multiple scans
    
    checkstring = file( find(file=='[',1,'first') : find(file==']',1,'last') );

    % let Matlab evaluate string
    try
        specs = eval(checkstring);
    catch
        fprintf('String cannot be evaluated. Check input!\n'); return;
    end
    
    % check for preceding zeros in all elements
    z1 = regexp(checkstring,'\d+'); % starting ind. of all elements
    z2 = regexp(checkstring,'(?=0*)[1-9]\d*');  % -"- without preceding zeros
    zeroadd = min(z2-z1);
    % Treat part before []
    filefront = file(1:find(file=='[',1,'first')-1);
    useprefix = ~isempty(filefront) && filefront(end)~=',';
    filefront = ssplit(filefront);
    if useprefix,   prefix = filefront{end}; filefront = filefront(1:end-1); 
    else            prefix = ''; end
    % Treat part after []
    fileback = file(find(file==']',1,'last')+1:end);
    usesuffix = ~isempty(fileback) && fileback(1)~=',';
    fileback = ssplit(fileback);
    if usesuffix,   suffix = fileback{1}; fileback = fileback(2:end);
    else            suffix = ''; end
    % Treat part in []
    for i = 1:length(specs)
        multifile{i} = [prefix, num2str(specs(i),['%0' num2str(ceil(log10(max(specs)+1)))+zeroadd 'd']), suffix];
    end

    filelist = [filefront, multifile, fileback];
else
    filelist = ssplit(file);
end
end

