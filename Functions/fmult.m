function handle=fmult(varargin)

% fsum(varargin)
% Creates the product function of the function handles in varargin

handle=@(a,b)prodfunc(a,b,varargin{:});