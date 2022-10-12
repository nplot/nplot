function handle=fsum(varargin)

% fsum(varargin)
% Creates the sum function of the function handles in varargin

handle=@(a,b)sumfunc(a,b,varargin{:});