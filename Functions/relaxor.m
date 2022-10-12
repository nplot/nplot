function [val paramnames paramnum description]=relaxor(param,x,opt)

% relaxor(param,x)
% "Single relaxor" f(x) = c * G*x / (x-w0)^2+G^2
% parameters: chi0, Gamma, w0

description='"Single Relaxor"';
paramnames = {'chi0', 'Gamma', 'w0'}; paramnum=3;

if isempty(param), val=[]; return; end

%-------------------------------------------------------
% Ab hier Definition der Funktion
%-------------------------------------------------------

chi0=param(1);
Gamma=param(2);
w0=param(3);

val = chi0 * x * Gamma ./ ((x-w0).^2+Gamma^2) ;
