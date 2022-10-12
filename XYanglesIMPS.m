function [xang,yang] = XYanglesIMPS(a4, a3, a5, roi, D, L)

% Calculate the twotheta angles

% P. St. 10/2011



a4=a4(:); a3=a3(:); a5=a5(:); roi=roi(:);
xang = [];
yang = [];


delta = atand( -(roi-5) * D .* cosd(a5)  ./ ( (roi-5) * D .*sind(a5) + L) ) ;
xang = a4 + delta;
yang = a3;


