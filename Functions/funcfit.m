function [erg,optparam,dpa,paramoutput,outputmessage] = funcfit(func, x, y, err, startparam, varindex, varargin)

% Syntax: [erg,optparam,dpa] = funcfit(func,x,y,err,startparam,varindex [,opt])

% can handle N-dim. x's if function is adapted (f: N-dim --> 1D)
% opt = 'nooutput' suppresses all screen output

[erg names paramnum des] = func([],[]); % Get function information
optparam = zeros(1,paramnum); dpa = optparam;
if ~any(strcmpi(varargin,'NoOutput')),  fprintf('%s . \n', des); end
if (nargin>5) && (max(size(startparam))<paramnum), fprintf('Not enough parameters!\n'); return; end
varindex(end+1:numel(startparam)) = true; % Ensure corect size, fill 1's
varindex = logical(varindex);
startparam = startparam(:)';

%% Look for constraints
constring = readinput('constraint',varargin);
if any(~varindex)
    % treat fixed variables as additional constraint
    if ~isempty(constring) && ~strcmp(constring(end),';'), constring(end+1) = ';'; end   % Ensure ; at the end
    ind = find(~varindex);
    for i = intersect(ind, 1:paramnum)
        constring = [constring, 'p' num2str(i,'%d') '=' num2str(startparam(i)) ';' ];
    end
end    
    
A = []; b = zeros(0,1);
if ~isempty(constring)
    [A, b] = constraintmatrix(constring);  % Obtain matrix form of contraint equations
    if isempty(A), fprintf('** Stop fitting due to problem with constraint string. \n'); return; end
    if rank(A) < rank([A,b]), fprintf('** Constraints cannot be fulfilled simultaneously! Stop fitting. \n'); return; end
    if size(A,2) > paramnum, fprintf('** Constraint string contains inexistent parameters! Stop fitting. \n'); return; end;
end

A(:,end+1:paramnum) = 0;
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
        aw = sum(((func(allpar,x)-y)./err).^2);
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
        dfp(:,n) = (func(optparam + h, x) - func(optparam - h, x)) / 2/h(vi(n));
    end
    alpha = dfp ./ repmat(err(:),1,numel(vi)); 
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

erg = func(optparam,x);

chisq  = sum(((y - erg)./err).^2); % chi_squared  
chisqN = chisq / (length(y) + 1e-8 - length(vi));  %normalized chi-squared

%% Output

paramoutput=[];
if ~any(strcmpi(varargin,'NoOutput'))
    if ~isempty(outputmessage), fprintf('Warning: %s\n',outputmessage); end
    fprintf('Chi2 = %10.5g \n',chisqN);
    fprintf('Result:\n');
    for i=1:1:paramnum
        fprintf('%20s :  %10.6g  ',names{i},optparam(i));
        paramoutput{i} = sprintf('%20s : ', names{i});
        if ~varindex(i)
            fprintf(' fixed ');
            paramoutput{i} = [paramoutput{i}, num2str(optparam(i),'%10.6g'), ' (fixed)'];
        else
            fprintf(' +/-   %10.5g',dpa(i));
%             paramoutput{i} = [paramoutput{i}, ' \pm ', num2str(dpa(i),'%10.5g')];
            paramoutput{i} =  [paramoutput{i}, converterror(optparam(i),dpa(i),'pm',2)];
            if ismember(i,vd), fprintf('  (constrained)'); paramoutput{i} = [paramoutput{i}, ' (constrained) ']; end
        end          
        fprintf('\n');
    end
end
if any(strcmpi(varargin,'printchi2')),  fprintf('Chi2 = %10.5g \n', chisqN); end


end  % funcfit
