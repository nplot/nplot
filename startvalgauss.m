function param = startvalgauss(x,y,n)

% function param = startval(x,y,n)
% Guess parameters for a fit to n Gaussians
% param = [ BG, x1, A1, Fwhm1, x2, A2, Fwhm2, ... ]

% P. Steffens, 10/2008


if numel(x)~=numel(y) || numel(y)<2, param=nan; return; end
np = numel(x);
[x,i] = sort(x,'ascend');
y = y(i); x=x(:);
if nargin<3, n=1; end

% First, smooth the data a bit
for ns = 1:3 % how often?
    yn(1) = (y(1) + y(2)/2) / 1.5;
    yn(2:(np-1)) = ( 0.5*y(1:(np-2)) + y(2:(np-1)) + 0.5*y(3:np) ) / 2;
    yn(np) = (y(np) + y(np-1)/2) / 1.5;
    y = yn;
end
y=y(:);

% Take lowest point as Background
param(1) = min(y);

subt = param(1)*ones(size(y));  % the already fitted portion
allowedx = true(size(x));       % where to look for a Max.

for ng = 1:n
    % Highest point defines x0 and Amplitude
    [A, i] = max(y(allowedx)-subt(allowedx));
    ind = find(allowedx);
    i = ind(i);
    x0 = x(i);
    allowedx(i) = false;
    % look 2 points left and right to get curvature (-->fwhm)
    ind = [max(i-2,1):(i-1), (i+1):min(i+2,np)];
    ind = ind(y(ind)-subt(ind)>0);
    fwhm = sqrt( 4*log(2) / sum( (log(A)-log(y(ind)-subt(ind)))./((x(ind)-x0).^2) ) * numel(ind) );
    
    param = [param, x0, A, fwhm];
    
    subt = subt + gaussA([x0,A,fwhm],x);
    exrange = min(fwhm/2, (max(x)-min(x))/10); % exclude a certain x-range from further peak search to avoid excessive overlap
    allowedx = allowedx & ((x > x0+exrange) | (x < x0-exrange)); 
    if ~any(allowedx), break; end
end
