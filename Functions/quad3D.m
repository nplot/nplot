function [val paramnames paramnum description]=quad3D(param,x)
% quad3D(param,x)
% Quadratic function in 3D f(x) = y0 + ...
% Parameters: y0, a1-3, b1-3, c1-3


description = 'quad3D';
paramnames = {'y0', 'a1', 'a2', 'a3', 'b1', 'b2', 'b3', 'c1', 'c2', 'c3'}; paramnum = 10;

if isempty(param),  val=[]; return; end

%-------------------------------------------------------
% Ab hier Definition der Funktion
%-------------------------------------------------------

val = param(1) + x(:,1)*param(2) + x(:,2)*param(3) + x(:,3)*param(4) + x(:,1).^2*param(5) + x(:,2).^2*param(6) + x(:,3).^2*param(7) + x(:,1).*x(:,2)*param(8) + x(:,1).*x(:,3)*param(9) + x(:,2).*x(:,3)*param(10);
