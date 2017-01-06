function [erg,optparam,dpa,paramoutput,outputmessage,chisqN] = funcfit(func, x, y, err, startparam, varindex, varargin)

% Syntax: [erg,optparam,dpa] = funcfit(func,x,y,err,startparam,varindex [,opt])
%
% Can handle N-dim. x's if function is adapted (f: N-dim --> 1D)
%
% Can fit multiple data sets, if x,y,err are cell arrays (i.e. x{i},y{i},err{i} = i-th dataset)
% In this case startparam can be length=m or simple array as for one dataset (->startparams equal for all)
% "common" parameter in varargin contains indices of common params,
% "constraint" can contain p2:3 for "p2 of data 3" (Recommended use for multiple datasets)
%
% options:
%   'nooutput' suppresses all screen output
%   'fullparamoutput' returns optparam and dpa in nd*np matrix containing all parameters (including common; for case of nd datasets)

% P. Steffens, 10/2016


erg=[]; optparam=[]; dpa=[]; paramoutput=[]; outputmessage=[]; chisqN=[];

%% Check input parameters 
allx=[]; ally=[]; allerr=[];
if iscell(y) || iscell(err) || iscell(x)
    if iscell(err) && iscell(y) && length(y)==length(err) && ((iscell(x) && length(x)==length(y)) || ~iscell(x))
        if ~iscell(x)
            for i=1:length(y), xx{i}=x; end; x=xx; clear xx;
        end
        if ~any(strcmpi(varargin,'NoOutput')),  fprintf('Simultaneous fitting of %d datasets.\n',length(y)); end
    else
        fprintf(2,'Error in ''funcfit.m'': inconsistent cell array sizes of x,y,err data for fitting of multiple data sets.\n'); return;
    end
    numdata = length(x);
    [~,datadim] = max(size(y{1})); % row or col vectors?
    try
        for i=1:numdata 
            allx=cat(datadim, allx, x{i}); ally=cat(datadim, ally, y{i}); allerr=cat(datadim, allerr, err{i}); 
            datalength(i) = numel(y{i}); 
        end
    catch
        fprintf(2,'Error: Dimension mismatch in x,y,err data passed to ''funcfit.m''.\n'); return;
    end
else
    numdata=1; [~,datadim] = max(size(y));
    allx=x; ally=y; allerr=err; datalength=numel(y);
end
% Check lengths of each
if size(allx,datadim)~= size(ally,datadim) || size(ally,datadim)~=size(allerr,datadim), fprintf(2,'Error: inconsistent sizes of x,y,err data on call to ''funcfit.m''\n'); return; end;


% Test call to function
try
    [erg, names, paramnum, des] = func([],[]); % Get function information
catch errmsg
    fprintf(2,'Error in ''funcfit.m'': Cannot call fit function. (%s)\n',errmsg.message); return;
end

% Sorting of parameters
paramassign = repmat(1:paramnum,numdata,1);
commonparams = readinput('common',varargin,'last');
if ~isnumeric(commonparams), fprintf(2,'Error in ''funcfit.m'': Value for ''common'' parameter is not numeric.\n'); return; end
paramassign(2:end,setdiff(1:paramnum,commonparams))=paramnum + reshape(1:(numdata-1)*(paramnum-numel(commonparams)), paramnum-numel(commonparams), numdata-1)';
% this means paramassign(n,m) = index of m-th parameter for n-th data set
allparamnum = max(paramassign(:));


% Assign startvalues
sp = nan(1,allparamnum);
if iscell(startparam)
    for i=1:length(startparam)
        if numel(startparam{i}) > paramnum, warning(['Too many start values for fit of dataset ',num2str(i)]); readlen=paramnum; else readlen = numel(startparam{i}); end
        sp( paramassign(i,1:readlen)) = startparam{i}(:)';
    end
elseif numel(startparam)<=paramnum
    for i=1:numdata
        sp(paramassign(i,1:numel(startparam))) = startparam(:)';
    end
elseif numel(startparam) == paramnum*numdata
    startparam = reshape(startparam',paramnum,numdata)';
    equal = true; for i=commonparams(:)', equal = equal && (numel(unique(startparam(:,i)))==1); end
    if ~equal, fprintf('Warning: non-equal start values for common parameters. (On call to ''funcfit.m'')\n'); end
    sp(paramassign) = startparam;
else
    if numel(startparam)>allparamnum, warning('Too many start values for fitting. Array must not have more entries than free parameters.'); end
    sp = startparam(:)'; sp = sp(1:allparamnum);
end
if any(isnan(sp)), warning('Start values have not been set for all fit parameters. (Zeros will be used.)'); end
sp(isnan(sp)) = 0;
startparam = sp;

optparam = zeros(1,allparamnum); dpa = optparam;
if ~any(strcmpi(varargin,'NoOutput')),  fprintf('%s . \n', des); end
% if (nargin>5) && (max(size(startparam))<paramnum), fprintf('Not enough parameters!\n'); return; end
varindex(end+1:numel(startparam)) = true; % Ensure correct size, fill 1's
varindex = logical(varindex);

%% Look for constraints
constring = readinput('constraint',varargin);

% translate input of form p2:3 using "real" parameter numbering
if ~isempty(constring)
    [st,en] = regexp(upper(strtrim(constring)), 'P\d+:\d+'); %find all Px entries in calcstring
    try for i=numel(st):-1:1
        pnum = strsplit(constring(st(i)+1:en(i)),':');
        pnum = [str2double(pnum{1}), str2double(pnum{2})]; 
        constring = [ constring(1:st(i)), num2str(paramassign(pnum(2),pnum(1))), constring(en(i)+1:end) ];
        end
    catch, fprintf(2,'Error in ''funcfit.m'' on interpretation of ''constraint'' string for multiple data sets.\n'); return; 
    end
end

if any(~varindex)
    % treat fixed variables as additional constraint
    if ~isempty(constring) && ~strcmp(constring(end),';'), constring(end+1) = ';'; end   % Ensure ; at the end
    ind = find(~varindex);
    for i = intersect(ind, 1:allparamnum)
        constring = [constring, 'p' num2str(i,'%d') '=' num2str(startparam(i)) ';' ]; %#ok<*AGROW>
    end
end    
    
A = []; b = zeros(0,1);
if ~isempty(constring)
    [A, b] = constraintmatrix(constring);  % Obtain matrix form of contraint equations
    if isempty(A), fprintf('** Stop fitting due to problem with constraint string. \n'); return; end
    if rank(A) < rank([A,b]), fprintf('** Constraints cannot be fulfilled simultaneously! Stop fitting. \n'); return; end
    if size(A,2) > allparamnum, fprintf('** Constraint string contains inexistent parameters! Stop fitting. \n'); return; end;
end

A(:,end+1:allparamnum) = 0;
% First, remove unneccesary lines
for i=size(A,1):-1:1
    if rank(A([1:i-1,i+1:end],:)) == rank(A)
        A = A([1:i-1,i+1:end],:);
        b = b([1:i-1,i+1:end]);
    end
end
% R = rank of matrix is number of dependent variables
% find R colums of A that make a square matrix of full rank
vi = []; vd = [];
Ain = []; Adep = [];
for i=size(A,2):-1:1
    if rank(A(:,1:i-1)) == rank(A)
        Ain = [Ain, A(:,i)];
        vi = [vi, i];       % vi is the set of independent variables
    else
        Adep = [Adep, A(:,i)];
        vd = [vd, i];       % vd is the set of dependent variables
    end
    A = A(:,1:i-1);
end
% Matrix A has been split into to matrices
% One can now rewrite A * x = b as Ain * x_i + Adep * x_d = b
% (where x_i is the vector of independent variables, etc.)
% -->  x_d = Adep^-1 * (b - Ain*x_i)


freeparam = startparam(vi);


%% Helper function
% Calculates the dependent parameters in order to fulfill the constraints
    function param = fulfillconstraints(param)
        xi = param(vi);
        if isempty(xi)
            xd = Adep \ b;
        else
            xd = Adep \ (b - Ain * xi');
        end
        param(vd) = xd;
    end


%% Helper function
% combines fixed, free and constrained variables and calls the fit function
% returns squared deviation (--> is the one to be minimized)
    function aw=abw(varpar)
        allpar(vi)   = varpar;
        allpar  = fulfillconstraints(allpar);
        aw = sum(((funceval(func,allpar,allx)-ally)./allerr).^2);
    end

%% Helper function
% Calls the function, piecewise for every dataset
    function f=funceval(func,param,x)
        for ii=1:numdata
            ind = sum(datalength(1:(ii-1))) + (1:datalength(ii));
            if datadim==1, xx=x(ind,:); else xx=x(:,ind); end
            f(ind) = func(param(paramassign(ii,:)),xx);
        end
        if datadim==1, f=f'; end
    end

%% Do Fit

[optfreeparam,fval,exitflag] = fminsearch(@abw,freeparam, optimset('display','off')); %#ok<ASGLU>
startparam(vi) = optfreeparam;
optparam = fulfillconstraints(startparam);

if exitflag~=1, outputmessage='Convergence NOT reached! '; else outputmessage = []; end
    

%% Calculate Errors

% Errors of independent variables
% Jacobian
if ~isempty(vi)
    
    smallval = 1e-6;
    for n=1:length(vi) %loop through varying parameters
        h = zeros(size(optparam)); %initialize h
        h(vi(n)) = smallval * optparam(vi(n)); 
        dfp(:,n) = (funceval(func,optparam + h, allx) - funceval(func,optparam - h, allx)) / 2/h(vi(n));
    end
    alpha = dfp ./ repmat(allerr(:),1,numel(vi)); 
    alpha = alpha' * alpha; 
    if rcond(alpha) > 1e-15
        C = inv(alpha); % raw correlation matrix
        dpa(vi) = sqrt(diag(C)); % error in a(j) is sqrt(Cjj)
        % Errors of dependent variables
        dpa(vd) = abs( Adep \ Ain * dpa(vi)' );
    else 
        dpa(1:end)=inf;
        outputmessage = [outputmessage, 'Matrix close to singular, cannot calculate errors! ']; 
    end
else
    dpa(vd) = 0;
end

erg = funceval(func,optparam,allx);
      
chisq  = sum(((ally - erg)./allerr).^2); % chi_squared  
chisqN = chisq / (length(ally) + 1e-8 - length(vi));  %normalized chi-squared

%% Output

paramoutput=[];
o = ~any(strcmpi(varargin,'NoOutput')); % flag for output on screen
if o && ~isempty(outputmessage), warning(outputmessage); end
if o, fprintf('Chi2 = %10.5g \n',chisqN); end
if o, fprintf('Result:\n'); end
outline=0;
for d=1:numdata
    if numdata>1 
        if o, fprintf('Dataset %d\n',d); end
        outline=outline+1; paramoutput{outline} = sprintf('Dataset %d\n',d);
    end
    for i=1:1:paramnum
        paramindex = paramassign(d,i); 
        if d>1 && paramindex<=paramnum, continue; end %(common param, do not print again)
        if o, fprintf('%20s :  %10.6g  ',names{i},optparam(paramindex)); end
        outline=outline+1; paramoutput{outline} = sprintf('%20s : ', names{i});
        if ~varindex(paramindex)
            if o, fprintf(' fixed '); end
            paramoutput{outline} = [paramoutput{outline}, num2str(optparam(paramindex),'%10.6g'), ' (fixed)'];
        else
            if o, fprintf(' +/-   %10.5g',dpa(paramindex)); end
    %             paramoutput{i} = [paramoutput{i}, ' \pm ', num2str(dpa(i),'%10.5g')];
            paramoutput{outline} =  [paramoutput{outline}, converterror(optparam(paramindex),dpa(paramindex),'pm',2)];
            if ismember(paramindex,vd)
                if o, fprintf('  (constrained)'); end
                paramoutput{outline} = [paramoutput{outline}, ' (constrained) ']; 
            end
        end  
        if ismember(i,commonparams)
            if o, fprintf(' (common)'); end
            paramoutput{outline} = [paramoutput{outline}, ' (common)'];
        end
        if o, fprintf('\n'); end
    end
end
if any(strcmpi(varargin,'printchi2')),  fprintf('Chi2 = %10.5g \n', chisqN); end %in case of only chi2 output

% parameter in matrix form if desired
if any(strcmpi(varargin,'fullparamoutput'))
    optparam = optparam(paramassign);
    dpa = dpa(paramassign);
end

end  % funcfit
