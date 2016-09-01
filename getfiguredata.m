function [list,ps] = getfiguredata(fig)
% Return data list of fcplot figure window

try
    if nargin==0, fig = gcf; end
    ps = guidata(fig);
    list = ps.datalist;
catch
    if nargin<1, fprintf('Error: The active window does not seem to contain valid data for getfiguredata.\n'); 
    else fprintf('Error: The input argument is not a window which contains valid data.\n'); end
    list = [];
    ps =[];
end