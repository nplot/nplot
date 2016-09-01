function [col, x,y,e] = rawplot(varargin)

% function rawplot(varargin)
% Plot raw data of a single scan file as color plot (.MULTI), lin or log
% with possibility to plot selected channels
% Examples:     rawplot 067353
%               rawplot log 067353      for log-colorplot
%               rawplot 067353 A4       to use a4-coordinate
%               rawplot 067353 A4 18:20 to show channels 18 to 20 separately

% P. Steffens, 10/2009


scannr=[];
logplot=false;
var = [];
channel = [];
offs=0;

if any(strcmpi(varargin{1},{'lin','log'}))
    logplot = strcmpi(varargin{1},'log');
    offs=1;
end

try
    scannr = varargin{1+offs};
    var = varargin{2+offs};
    channel = varargin{3+offs};
catch
end
    

scan = tasread(scannr,'download');
figure

if ~isempty(channel), subplot(1,2,1); end

if ~isempty(var) && isfield(scan.DATA,var)
    xvals = scan.DATA.(var);
    xlab = var;
elseif ~isempty(var) && strcmpi(var,'a4real')
    xvals = getvar(scan,'a4');
    xlab = 'A4';
else
    xvals = scan.DATA.PNT;
    xlab = 'PNT';
end

if logplot
    col = log(max(1,scan.MULTI([1:end,end],[1:31,31])));
else
    col = scan.MULTI([1:end,end],[1:31,31]);
end

pcolor(.5:31.5, (xvals([1,1:end])+xvals([1:end,end]))/2, col);
xlabel('channel');
ylabel(xlab);
title(scannr);

if ~isempty(channel)
    subplot(1,2,2);
    [x,y,e] = rawchannelplot(scannr, channel, var, gca);
end

if nargout==0, clear col; end
