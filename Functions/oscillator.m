function [val paramnames paramnum description]=oscillator(param,x,opt)

% oscillator(param,x)
% Harmonic oscillator, f(x) = a* G*x*w0^2 / (x^2-w0^2)^2+(G*x)^2
% parameters: Amplitude, Gamma, w0


description='Harmonic Oscillator';
paramnames = {'Ampl.', 'Gamma', 'w0'}; paramnum=3;

if isempty(param),  val=[]; return; end

%-------------------------------------------------------
% Ab hier Definition der Funktion
%-------------------------------------------------------

chi0=param(1);
Gamma=param(2);
w0=param(3);

val = chi0 * x * Gamma *w0^2 ./ ((x.^2-w0^2).^2+Gamma^2*x.^2) ;
