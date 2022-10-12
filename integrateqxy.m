function [x,y,dy,H,K,L] = integrateqxy(data, startpoint, endpoint, corner, npoints)

% [x,y,dy,H,K,L] = integrateqxy(data, startpoint, endpoint, corner, npoints)
%
% integrate the data set "data" over a Q-region
% startpoint:   starting point of the scan
% endpoint:     end point of the scan
% corner:   	third point, defining one of the corners of the integration region
%               if a single number is given, interprete it as the total thickness of a rectangular integration region in recipr. Angstrom
% npoints:      number of points
%
% Convention for planes in Q-Sapce:
% For startpoint, endpoint and corner, a three-element vector will be interpreted as HKL, 
% while a two-element vector will be interpreted as (qx,qy)
%
% Return values:
% x:    linear coordinate of scan in A^-1 (or other if not in Q-plane) in direction of largest variation
% y,dy: values and errors
% H,K,L:coordinates of scan points

% P. Steffens, 05/2009 - 01/2017

qztolerance = 0.05;  %Tolerance for vertical momentum component
qzwarning = false;

startpoint = startpoint(:); endpoint = endpoint(:); corner = corner(:);

if ~any(strcmpi(data.coordtype,{'QXY','QPLANE','ANGLES'}))
    qplane = false; 
else
    data = coordtransform(data,'qxy');
    qplane=true;
   

    if any([numel(startpoint),numel(endpoint),numel(corner)]==3) && any([numel(startpoint),numel(endpoint),numel(corner)]==2)
        fprintf('Two-element vectors are interpreted as (qx,qy) and three-element vectors are interpreted as (H,K,L).\n');
    end

    try
        UB = UBmatrix( data.sampleinfo.lattice, data.sampleinfo.ax, data.sampleinfo.bx);
        if numel(startpoint)==3  % interprete as HKL
            [startpoint(1),startpoint(2),qz] = calcQS(startpoint(1),startpoint(2),startpoint(3),UB);
            if isfield(data,'QVERT'),  if abs(qz-data.QVERT) > qztolerance, qzwarning=true; end, end %check qvert
        end
        if numel(endpoint)==3  % interprete as HKL
            [endpoint(1),endpoint(2),qz] = calcQS(endpoint(1),endpoint(2),endpoint(3),UB);
            if isfield(data,'QVERT'),  if abs(qz-data.QVERT) > qztolerance, qzwarning=true; end, end %check qvert
        end
        if numel(corner)==3  % interprete as HKL
            [corner(1),corner(2),qz] = calcQS(corner(1),corner(2),corner(3),UB);
            if isfield(data,'QVERT'),  if abs(qz-data.QVERT) > qztolerance, qzwarning=true; end, end %check qvert
        end
        hklok= true;
    catch
        fprintf('The transformation of H,K,L-values is not possible.\n');
        hklok=false;
    end

    if qzwarning, fprintf(['Warning: At least one of the given H,K,L is more than %g A^-1 away from the Qvert of the data set.\n' ...
                           'This is ignored in the integration (use projection on data plane). Please check.\n'], qztolerance); end

end
                   
pathvector = endpoint(1:2)-startpoint(1:2);
pv = [-pathvector(2), pathvector(1)];
pv = pv / sqrt(pv*pv'); 
% unit vector perpendicular to pathvector

if numel(corner)==1
    sidevector = corner/2 * pv';
else
    sidevector = corner(1:2)-startpoint(1:2);
end

if qplane && hklok
    [starthkl(1),starthkl(2),starthkl(3)] = calcHKL(startpoint(1),startpoint(2),data.QVERT,UB);
    [endhkl(1),endhkl(2),endhkl(3)] = calcHKL(endpoint(1),endpoint(2),data.QVERT,UB);
end
    
fprintf('Integration in %d steps from (x,y)=(%g,%g) to (%g,%g) A^-1', npoints, round(1E4*startpoint(1:2))/1E4, round(1E4*(startpoint(1:2)+pathvector))/1E4);
if qplane && hklok, fprintf(', resp. (H,K,L)=(%g,%g,%g) to (%g,%g,%g)', round(1E3*starthkl)/1E3, round(1E3*endhkl)/1E3); end
fprintf('.\nThickness of integration range: %g A^-1',abs(2*pv*sidevector));   
if abs(sidevector(:)'*pathvector) >1E-5, fprintf(' (non-rectangular)');end
fprintf('.\n');

integ = integratepatch(data.faces, data.vertexlist, data.valuelist(:,1), data.valuelist(:,2),  ... 
                       startpoint(1:2)'-sidevector(:)', [pathvector(:)'; 2*sidevector(:)'], [npoints,1], data.raw);
if isempty(integ), fprintf('Integration was not successful.\n'); x=[]; y=[]; dy=[]; return; end
result = integ{1};

xdim = 1;
if pathvector(2)>pathvector(1), xdim = 2; end
% if y-component of pathvector larger than x-component, return qy-values, else qx

x  = startpoint(xdim) + ((1:npoints)'-0.5) * pathvector(xdim)/npoints;
y  = result(:,1);
dy = result(:,2);

if nargout>3
    if qplane && hklok
        [H,K,L] = calcHKL(startpoint(1) + ((1:npoints)'-0.5) * pathvector(1)/npoints, startpoint(2) + ((1:npoints)'-0.5) * pathvector(2)/npoints, data.QVERT, UB);
    else
        fprintf('Cannot calculate H,K,L for output!\n'); H=[]; K=[]; L=[];
    end
end

