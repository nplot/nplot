function [func, msg, paramnames, paramnum, description] = makefunction(func)

% Transform input into function handle and test validity
% Returns valid function handle or []

% P. Steffens 2/2015


% ** Treat "+" etc.

msg = [];
paramnames=[]; paramnum=[]; description=[];

if ischar(func) % Ensure "@" for function handles ** überlegen wg fsum, fmult
%    summands = strsplit(func,'+'); 
   % Replaces Matlab strsplit (not in older versions)
    plind = [0, find(func=='+'), numel(func)+1];
    for f=1:numel(plind)-1, summands{f} = func(plind(f)+1:plind(f+1)-1); end %#ok<AGROW>
    
    for f=1:length(summands)
        if summands{f}(1) ~= '@' && isempty(strfind(func,'fsum')) && isempty(strfind(func,'fmult')) 
            func = ['@',summands{f}]; 
        else func = summands{f}; end
        try eval(['fc{' num2str(f) '}=' func ';']);
        catch, msg = 'Error: problem on calling function'; func = []; return; end
    end
    if f==1, func = fc{1}; else func=fsum(fc{:}); end
end
try % function test call 
    [~, paramnames, paramnum, description] = func([],[]); %#ok<RHSFN>
catch
    func = []; msg = 'Error: problem on calling function';
end