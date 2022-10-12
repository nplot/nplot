function varargout = getoption(varargin)

% Returns one or more variables from the options file
% varargin : Give variable names as strings
% varargout: Corresponding list of values
% If varargin terminates like 'check',input , where input is a string or
% cell array of strings that can be interpreted by readinput(), input is
% analyzed for a corresdponding value. If found, this value is returned
% for the corresponding variables, otherwise those from the options.m file.

% P. Steffens, 11/2009

if nargout == 0, return; end

% If there is some extra input, analyze this first
if nargin>1 && strcmpi('check',varargin{nargin-1}) && ~isempty(varargin{nargin})
    inputstring = varargin{nargin};
    if ~iscell(inputstring), inputstring = {inputstring}; end
    for i=1:min(nargin-2, nargout)
        varargout{i} = readinput(varargin{i},inputstring);
        found(i) =  ~isempty(varargout{i});
    end
    if all(found), return; end
    lookfor = find(~found(:)');
else
    lookfor = 1:nargout;
end


% If not all values have been assigned, go on:
% check if persistent variable exists, if not, load options file

persistent optionsvalues;

% 'ooptionsvalues' contains a field for every variable in options.m
% It contains also:
% .lastupdate:           serial date of last update of options.m on disk
% .lastversionchecktime: serial date(time) of last check if optionsvalues is up to date
% .optionsfilename:      name and path of options.m used for last update 

dnn = now;

if ~isempty(optionsvalues) && dnn - optionsvalues.lastversionchecktime > 5E-5  
% (do this check only once every 5sec. to avoid too many calls to the following lines)
    % Check if current file on disk is the same version as the values in memory
    % If no, clear memory.
    dinfo = dir('options.m');
    if isempty(dinfo), dinfo = dir(optionsvalues.optionsfilename); end
    optionsvalues.lastversionchecktime = dnn;
    if isempty(dinfo) || dinfo.datenum ~= optionsvalues.lastupdate || ~strcmpi(which('options'),optionsvalues.optionsfilename)
        optfile = which('options');
        if isempty(optfile) 
            error('Error: file options.m not found. Please create a local copy or include in Path.\n');
        end
        fprintf('File options.m has changed. Reload from %s.\n',which('options'));
        optionsvalues = [];
    end
end

if  isempty(optionsvalues)  % If not in memory, load options.m 
    name = which('options');
    dinfo = dir(name);
    optionsvalues.lastupdate = dinfo.datenum; % modification date of this options.m
    optionsvalues.optionsfilename = name;
    optionsvalues.lastversionchecktime = dnn;
    clear name dinfo
    varnamesold = who;
    
    options; % call to options.m to load all variables into workspace
    
    varnamesnew = who;
    for i=1:length(varnamesnew) % save all new variables in structure 'optionsvalues'
        if ~any(strcmpi(varnamesnew{i},varnamesold))
            eval(['optionsvalues.(varnamesnew{' num2str(i) '}) = eval(varnamesnew{' num2str(i) '});']);
        end
    end
end

% Assign the ouptput variables.

for i = lookfor
     
    try
        if all(varargin{i}~='.')
            varargout{i} = optionsvalues.(varargin{i});
        else
            istr = varargin{i};
            pp = find([istr,'.']=='.');
            varargout{i} = optionsvalues.(istr(1:pp(1)-1));
            for pi = 1:numel(pp)-1
                varargout{i} = varargout{i}.(istr(pp(pi)+1:pp(pi+1)-1));
            end
        end
        
%   * replaced the following         
%         [name,rem] = strtok(varargin{i},'.');
%         varargout{i} = optionsvalues.(name);
%         while ~isempty(rem)  % in case varargin{i} asks for the field of a struct, like aaa.bbb
%             [name,rem] = strtok(rem,'.'); %#ok<STTOK>
%             varargout{i} = varargout{i}.(name);
%         end
    catch
        varargout{i}=[];
        fprintf('Warning: Variable name  %s  not found in options file during call to getoption.\n',varargin{i});
        % Return [] if not found
    end
    
end