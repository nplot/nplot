function matches = regexpmatch(text,expr)
% Simulates the behavior of regexp(..,..,'match'), which is not implemented
% in older Matlab-Versions
[st,en] = regexp(text,expr);
matches=[];
for ii=1:numel(st)
    matches{ii}=text(st(ii):en(ii));
end