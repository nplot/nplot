function [func, msg, paramnames, paramnum, description] = makefunction(func)

% Transform input into function handle and test validity
% Returns valid function handle or []

% P. Steffens 2/2015


% ** Treat "+" etc.

msg = [];
paramnames=[]; paramnum=[]; description=[];
if ischar(func) % Ensure "@" for function handles ** überlegen wg fsum, fmult
    if func(1) ~= '@' && isempty(strfind(func,'fsum')) && isempty(strfind(func,'fmult')), func = ['@',func]; end  
    eval(['func=' func ';']);
end
try % function test call 
    [~, paramnames, paramnum, description] = func([],[]); %#ok<RHSFN>
catch
    msg = 'Error: problem on calling function'; 
    func = [];
end