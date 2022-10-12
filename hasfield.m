function tf = hasfield(S, C)

% Overrides Matlab function 'isfield' for reasons of compatibility with
% former Matlab versions (which did not permit C as cell array)

tf = [];
if iscellstr(C)
    for i=1:length(C)
        tf(i) = isfield(S, C{i});
    end
else
    tf = isfield(S,C);
end
