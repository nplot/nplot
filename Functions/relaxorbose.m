function [val paramnames paramnum description]=relaxorbose(param,x,opt)

% relaxor(param,x) .* bose
% "Single relaxor"*bose f(x) = c * G*x / (x-w0)^2+G^2 * bose(w,T)
% parameters: chi0, Gamma, w0, T

description='"Single Relaxor"*Bose';
paramnames = {'chi0', 'Gamma', 'w0', 'T'}; paramnum=4;

if isempty(param), val=[]; return; end

%-------------------------------------------------------
% Ab hier Definition der Funktion
%-------------------------------------------------------

chi0=param(1);
Gamma=param(2);
w0=param(3);
T = param(4);

val = chi0 * x * Gamma ./ ((x-w0).^2+Gamma^2) .* bose(x,T) ;

val(isnan(val))=0;
