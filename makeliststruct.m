function data = makeliststruct(coords, values, varargin)

% Create a standard compatible list structure manually
% coords :  linear list of coordinates, dim(coords)=[npoints, nDimensions]
% values :  linear list of values and errors, dim(values)=[npoints, 2]
% varargin: pairs of field names and values


data.coordlist = coords;
data.valuelist = values;
data.type = 'GENERAL';
data.coordtype = 'GENERAL';
data.dataname = '';
data.expname = '';
data.raw = 0;
data.constants = {};
try
    for f=1:2:nargin-2
        data.(varargin{f}) = varargin{f+1};
    end
catch
    fprintf('Error during assignment of input to makeliststruct.\n');
end

