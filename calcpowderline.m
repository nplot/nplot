function [vertices, connection] = calcpowderline(data, dval, type)

% Calculate the cut of a powder line through a 2D data slice in n-dim coordinate space 
% dval:  d-spacing of the powder reflection
% type:  1 incoherent on Ana (kf'=ki)
%        2 incoherent on Mono (ki'=kf)

% P. Steffens, 04/2009 - 08/2014


% helper function
function c = getconstant(field)
    if ~isfield(data,field)
        fprintf('Error in powder line calculation: field %s is needed, but not existent in given data.\n', field);
        error('Field not existent');
    else
        c = data.(field);
    end
end



% main function

vertices = [];
connection = [];

if type == 2
    try kf = getconstant('KF'); catch return; end 
    % (Make use of fact that no kf variation is ever possible)
    lambda = 2*pi/kf;
end

if strcmpi(data.coordtype, 'angles')
    % need Ki, Qvert
    try ki = getconstant('KI');  Qvert = getconstant('QVERT'); catch return; end
    if type==1, lambda = 2*pi/ki; end
    
    pos = abs(data.vertexlist(:,1)) - acosd( (1 - lambda^2/dval^2/2) / sqrt(1 - Qvert^2/ki^2) );
    
elseif strcmpi(data.coordtype, 'anglesqz')
    % need Ki
    try ki = getconstant('KI');  catch return; end
    if type==1, lambda = 2*pi/ki; end
    
    pos = abs(data.vertexlist(:,1)) - acosd( (1 - lambda^2/dval^2/2) ./ sqrt(1 - data.vertexlist(:,3).^2/ki^2) );
    
elseif strcmpi(data.coordtype, 'qxy')
    % need Ki, Kf, Qvert
    try ki = getconstant('KI'); kf = getconstant('KF');  Qvert = getconstant('QVERT'); catch return; end
    if type==1, lambda = 2*pi/ki; end
    
    % The scattering angle 2th is Q^2 = ki^2+kf^2-2kikf*cos(2th). Note: 2th~=a4 in 3D
    % Use cos(2th) = 1-2sin^2(th).
    % The |q| for ki=kf=k is q=2k*sin(th). Equals this to Q(pow)=pi/d ?
    
    pos = data.vertexlist(:,1).^2 + data.vertexlist(:,2).^2 + Qvert^2 - ki^2 - kf^2 + 2*ki*kf*(1 - lambda^2/dval^2/2);
    
elseif strcmpi(data.coordtype, 'qxyz')
    % need Ki, Kf
    try ki = getconstant('KI'); kf = getconstant('KF'); catch return; end
    if type==1, lambda = 2*pi/ki; end
    
    pos = data.vertexlist(:,1).^2 + data.vertexlist(:,2).^2 + data.vertexlist(:,3).^2 - ki^2 - kf^2 + 2*ki*kf(1 - lambda^2/dval^2/2);
    
elseif any(strcmpi(data.coordtype, {'anglesenergy','a4energy'}))
    % use here again that kf is always constant (otherwise would have to distinguish ki/kf const. scans)
    % need Kf, Qvert
    try kf = getconstant('KF'); Qvert = getconstant('QVERT');  catch return; end
    [hbar,mass_n,meVJ] = getoption('hbar', 'mass_n', 'meVJ');
    ki = sqrt(2*mass_n/meVJ/hbar^2 * data.vertexlist(:,2) + kf^2*1e20) *1e-10;
    if type == 1, lambda = 2*pi./ki; end
    
    pos = abs(data.vertexlist(:,1)) - acosd( (1 - lambda.^2/dval^2/2) ./ sqrt(1 - Qvert^2./ki.^2) );

else
    
    fprintf('Error in powder line calculation: coordinate type is not recognized.\n');
    
end

[vertices, connection] = createintersection(data.vertexlist, data.faces, pos, 2);

end
    
    
    