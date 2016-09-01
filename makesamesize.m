function [ok,varargout] = makesamesize(varargin)

% Check argument list and make arrays the same size (columns), if possible

varargout = varargin;
z=length(varargin);
for i=1:z
    gr(i) = numel(varargin{i});
end
if any(gr==0) || numel(unique(gr(gr>1)))>1
    fprintf('Inconsistent array size!\n');
    ok=false;
    return;
end

gg = max(gr);
for i=1:z
    arr = varargin{i};
    arr =arr(:);
    if gr(i)==1
        arr = arr * ones(gg,1);
    end
    varargout(i) = {arr};
end
    
ok=true;