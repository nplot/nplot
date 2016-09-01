function scaled = scaledata(data,factor)

% Returns a multiple of 'data'
%
% P. Steffens, 08/2008

scaled = data;

if isfield(scaled, 'monitorlist'), scaled = rmfield(scaled,'monitorlist'); end
scaled.raw = false;

if isfield(scaled,'dataname'), scaled.dataname = [num2str(factor) '*(' scaled.dataname ')']; end

scaled.valuelist(:,1) =     factor  * scaled.valuelist(:,1);
scaled.valuelist(:,2) = abs(factor) * scaled.valuelist(:,2);

