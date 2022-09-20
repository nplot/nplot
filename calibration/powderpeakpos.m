function [val, paramnames, paramnum, description] = powderpeakpos(param,x,opt)

% function used during fitting
% calculate expected a4 value of powder Bragg peak including angular offsets

% powderpeakpos(param,x,opt)

% Parameters: da2, da4
% x = [a4theor,kfix,dm,a2]


description = 'powder Bragg peaks';
paramnames = {'da2', 'da4'}; paramnum = length(paramnames);

if isempty(param),  val=[]; return; end

%-------------------------------------------------------
% Function definition
%-------------------------------------------------------


a2 = x(:,4);
dm = x(:,3);
kfix = x(:,2);
a4theor = x(:,1);

da2 = param(1);
da4 = param(2);

lambda = abs(2*dm.*sind((a2+da2)./2)); % "real" wavelength

% d-value of powder reflection, deduced from nominal a4 and lambda
d = pi./kfix./sind(a4theor/2);  % (is negative if a4theor negative)

val = 2 * asind(lambda./(2*d)) - da4; % expected measured a4-pos of peak
