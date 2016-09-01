function [val paramnames paramnum description]=powfit(param,x,opt)

% powfit (param, x)
% calc a4 position of powder peaks

description = 'powfit';
paramnames  = {'da1','da4'}; paramnum=2;

if isempty(param),  val=[]; return; end


%-------------------------------------------------------
% Ab hier Definition der Funktion
%-------------------------------------------------------

da1 = param(1);
da4 = param(2);


val = 2 * asind(3.355./x(:,2) .* sind(x(:,1) + da1)) + da4;

