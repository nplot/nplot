function [xx,yy,ee] = rawchannelplot(scannr,chan,var,ax)

% function [xx,yy,ee] = rawchannelplot(scannr,chan,var,ax)
% Plot one or several channels of data scannr
% If var is given, use this column as x-values

% P.Steffens, 10/2009

if nargin<4
    figure
else
    axes(ax);
end

scan = tasread(scannr,'download');

if isstr(chan), chan=str2num(chan); end
xoffset = zeros(size(chan));

if nargin>2 && isfield(scan.DATA,var)
    xvals = scan.DATA.(var);
    xlab = var;
elseif nargin>2 && strcmpi(var,'a4real')
    xvals = getvar(scan,'a4');
    xoffset = (chan-16) * 2.5;
    xlab = 'Scattering angle';
else
    xvals = scan.DATA.PNT;
    xlab = 'PNT';
end

style = {'ok','or','ob','og','oc','.k','.r','.b','.g','.c'};

xx = []; yy = []; ee = [];

for i=1:numel(chan)
    errorbar(xvals + xoffset(i), scan.MULTI(:,chan(i)), sqrt(scan.MULTI(:,chan(i))), style{mod(i,length(style))+1});
    xx = [xx; xvals(:)+xoffset(i)];
    yy = [yy; scan.MULTI(:,chan(i))];
    ee = [ee; sqrt(scan.MULTI(:,chan(i)))];
    leg{i} = ['Channel ' num2str(chan(i),'%d')];
    hold on
end
hold off

legend(leg);
xlabel(xlab);

if isfield(scan.PARAM,'MN'), ylabel(['counts (MN ' num2str(scan.PARAM.MN) ')']);
elseif isfield(scan.PARAM,'TI'), ylabel(['counts (t=' num2str(scan.PARAM.TI) 's)']);
else ylabel('counts');
end
title(scannr);

if nargout<1, clear xx; end
