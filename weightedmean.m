function [mval,msig] = weightedmean(val,sig,weight)

% Weighted average
% 
% P. Steffens 01/2008

if nargin==2 ; weight= max(1E-20, 1./sig.^2); end  % ** Max to be revised

sumwgt = sum(weight(:));

mval = val(:)' * weight(:) / sumwgt;

msig = sqrt ( weight(:)'.^2 * sig(:).^2 ) /  sumwgt;