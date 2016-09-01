function [xang,yang]=XYangles(a4,gfc,ch,psi)

% Calculate angles for FC data treatment
% xang is pseudo detector scattering angle (replaces a4)
% yang is pseudo sample rotation
%
% P. Steffens 10/2008

% a4,gfc,psi must have same number of elements


[alpha,alphaP]=Scatt_Angle(a4,gfc,ch);

[alpha31,alphaP31] = Scatt_Angle(a4,gfc,31);

xang = alphaP;

yang = (psi(:) + alphaP31) * ones(size(ch(:)')) - 90;

% The following is transferred to bringtorange.m, and called from outside

% % Problem: In a scan, the psi-value might "jump" from -180 to 180 degrees
% % and thus give a discontinuous set of yang, which is a problem for the
% % Voronoi cells, interpolation, etc.
% % If this is the case, add 360 degrees to a part of the yang's.
% if max(yang(:))-min(yang(:)) > 360-4 %(4 is arbitrary, might take stdcell size)
%     % Try to find the 'hole' in the yang data
%     yangsort = sort(unique(yang));
%     [m,i] = max(yangsort(2:end)-yangsort(1:end-1));
%     if m>6 %(6 again arbitrary)
%         % A hole exists
%         hole = (yangsort(i)+yangsort(i+1))/2;
%         % Add 360 to values left of the hole
%         yang(yang<hole) = yang(yang<hole)+360;
%     end
% end






