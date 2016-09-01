function [zi,errzi] = linearinterpolation2D(datalist,xi,yi,maxsize)

% Triangle-based linear interpolation in 2D
%
% Uses Delaunay triangulation
% Slightly modified MATLAB code from "griddata":
% Uses additionally the errorbars of the given points x,y as wheigt upon averaging 
% maxsize (optional) determines max. size of delaunay triangles in the respective dimension (standard from options file if not given)
%
% P. Steffens, 06/2008


siz = size(xi);
zi=[]; errzi=[];
xi = xi(:); yi = yi(:); % Treat these as columns
[ok,xi,yi] = makesamesize(xi,yi); if ~ok, fprintf('Wrong input format!\n'); return; end
x = datalist.coordlist(:,1); y  = datalist.coordlist(:,2); 
z = datalist.valuelist(:,1); dz = datalist.valuelist(:,2); mo = datalist.monitorlist(:,2);


if isfield(datalist,'delaunaytri')
    tri = datalist.delaunaytri;
else
    % Triangularize the data
    tri = delaunayn([x y]);
end
    
% ******* added code (not in griddata) *****************************
% Delete those triangles that are too large (unphysically large)
stdcell = getoption('stdcell');
if nargin<4
    maxsize = stdcell.(upper(datalist.coordtype));
end
% i=1;
% while i<=size(tri,1)  % 
%     if    (max(datalist.coordlist(tri(i,:),1)) - min(datalist.coordlist(tri(i,:),1))) > maxsize(1) || ...
%           (max(datalist.coordlist(tri(i,:),2)) - min(datalist.coordlist(tri(i,:),2))) > maxsize(2)
%         tri = tri([1:(i-1),(i+1):end],:);
%     else
%         i=i+1;
%     end
% end

toolarge =   ( max(reshape(datalist.coordlist(tri,1),size(tri)),[],2) - min(reshape(datalist.coordlist(tri,1),size(tri)),[],2) > maxsize(1) ) ...
           | ( max(reshape(datalist.coordlist(tri,2),size(tri)),[],2) - min(reshape(datalist.coordlist(tri,2),size(tri)),[],2) > maxsize(2) );
% ******************************************************************
    
if isempty(tri),
  warning('linearinterp:CannotTriangulate','Data cannot be triangulated.');
  zi = repmat(NaN,size(xi));
  return
end

% Find the nearest triangle (t) for each xi,yi
t = tsearch(x,y,tri,xi,yi);

badts  = ismember(t, find(toolarge)) | isnan(t) ; % logical index into t for those triangles that are too large
goodts = ~badts;

% Only keep the relevant triangles.
out = find(badts);
if ~isempty(out), t(out) = ones(size(out)); end
tri = tri(t,:);

% Compute Barycentric coordinates (w).  P. 78 in Watson.
del = (x(tri(:,2))-x(tri(:,1))) .* (y(tri(:,3))-y(tri(:,1))) - ...
      (x(tri(:,3))-x(tri(:,1))) .* (y(tri(:,2))-y(tri(:,1)));
w(:,3) = ((x(tri(:,1))-xi).*(y(tri(:,2))-yi) - ...
          (x(tri(:,2))-xi).*(y(tri(:,1))-yi)) ./ del;
w(:,2) = ((x(tri(:,3))-xi).*(y(tri(:,1))-yi) - ...
          (x(tri(:,1))-xi).*(y(tri(:,3))-yi)) ./ del;
w(:,1) = ((x(tri(:,2))-xi).*(y(tri(:,3))-yi) - ...
          (x(tri(:,3))-xi).*(y(tri(:,2))-yi)) ./ del;
w(out,:) = zeros(length(out),3);


% Until here, code is essentially identical to original. Now, use 
% weighted means with error bars instead.

for i=1:size(tri,1)
    if datalist.raw ==1
     %   [zi(i), errzi(i)] = weightedmean ( z(tri(i,:)), dz(tri(i,:)),  w(i,:) .* z(tri(i,:))' ./ dz(tri(i,:))'.^2 ); %#ok<AGROW>
        [zi(i), errzi(i)] = weightedmean ( z(tri(i,:)), dz(tri(i,:)),  w(i,:) .* mo(tri(i,:))' ); %#ok<AGROW>
    else
        [zi(i), errzi(i)] = weightedmean ( z(tri(i,:)), dz(tri(i,:)),  w(i,:) ./ dz(tri(i,:))'.^2  ); %#ok<AGROW>
    end
end

zi = reshape(zi,siz);
errzi = reshape(errzi,siz);

if ~isempty(out), zi(out) = NaN; errzi(out) = NaN; end
%------------------------------------------------------------