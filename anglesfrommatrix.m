function [a3,gu,gl,a3p,varargout] = anglesfrommatrix(varargin)

% Calculate a3,gu,gl,a3p from given rotation matrix
%
% Input arguments must be:
% M :                 rotation matrix
% opt (optional) :    define type of solution (e.g. 'A3p=10' etc.)
% zeroval (optional): angles zeros
% config:             spectrometer details (sign_gu etc.)
%
% P. Steffens 02/2008

%options;

M = varargin{1};

if nargin < 2 , opt = getoption('stdscanmode');  else opt = varargin{2}; end  % if scan mode not specified, use standard
opt = upper(opt);

if nargin < 3 , zerovals = getoption('stdzeros'); else zerovals = varargin{3}; end % if not specified, take standard zeros

if nargin < 4
    [config.za3_0, config.zgu_0, config.zgl_0, config.GUsign, config.GLsign] = getoption('za3_0', 'zgu_0', 'zgl_0', 'GUsign', 'GLsign');
else
    config = varargin{4};
end
za3_0 = config.za3_0; zgu_0 = config.zgu_0; zgl_0 = config.zgl_0;
GUsign = config.GUsign; GLsign = config.GLsign;

U0 = zeromatrix(zerovals, config);

if nargout>4, varargout = {opt,zerovals,config}; end

% Now, it depends which solution is desired (use of A3/A3')

if strcmp(opt(1:4),'A3P=')
    % A3' is at fixed value, solve for a3,gu,gl
    a3p = str2double(opt(5:end));
    a3p_phys = pi/180 * (a3p - zerovals.a3p);
    RA3P= [cos(a3p_phys),  -sin(a3p_phys),  0; sin(a3p_phys),  cos(a3p_phys),  0; 0, 0, 1];

    R = M * U0 * RA3P';

    a3_phys = 180/pi * atan2( R(2,1), R(1,1));
    gu_phys = 180/pi * atan2( R(3,2), R(3,3));
    gl_phys = 180/pi * atan2(-R(3,1), sqrt(R(3,2)^2+R(3,3)^2));
    a3 =  a3_phys + (zerovals.a3 - za3_0);
    gu = (gu_phys + (zerovals.gu - zgu_0)*GUsign) * GUsign;
    gl = (gl_phys + (zerovals.gl - zgl_0)*GLsign) * GLsign;  
end

if strcmp(opt(1:3),'A3=')
    % A3 is at fixed value, solve for a3',gu,gl
    a3 = str2double(opt(4:end));
    a3_phys = pi/180 * (a3 - zerovals.a3 - za3_0);
    RA3= [cos(a3_phys),  -sin(a3_phys),  0; sin(a3_phys),  cos(a3_phys),  0; 0, 0, 1];

    R = RA3' * M * U0 ;

    a3p_phys= 180/pi * atan2( R(2,1), R(2,2));
    gu_phys = 180/pi * atan2(-R(2,3), sqrt(R(2,1)^2+R(2,2)^2));
    gl_phys = 180/pi * atan2( R(1,3), R(3,3));
    gu = (gu_phys  + (zerovals.gu - zgu_0)*GUsign) * GUsign;
    gl = (gl_phys  + (zerovals.gl - zgl_0)*GLsign) * GLsign;  
    a3p=  a3p_phys +  zerovals.a3p;
end
