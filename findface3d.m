function [nface,nslice,nfaceonslice] = findface3d (ax, plotcoordlist, plotvertexlist, datalist)

% Identify face and slice on which a mous click in 3d is lying.
% nan if empty
% nface: index in global list, 
% nfaceonslice: face index on slice nslice

% replaces findface.m

% P. Steffens, 08/2014



nface = nan; nslice = nan; nfaceonslice = nan;
if size(plotcoordlist,2)==2,  plotcoordlist(:,3)=0;  end
if size(plotvertexlist,2)==2, plotvertexlist(:,3)=0; end
ndata = length(datalist);
if ndata==1 && ~iscell(datalist), m=datalist; clear datalist; datalist{1}=m; clear m; end % if not, make cell array

% count nb of points and vertices in each dataset
for sl = 1:ndata, nbpoints(sl) = size(datalist{sl}.faces,1); end %#ok<AGROW>
for sl = 1:ndata, nbverts(sl) = size(datalist{sl}.vertexlist,1); end %#ok<AGROW>

% Transform all coordinates to the viewing plane system

point     = get(ax, 'CurrentPoint');    % mouse click position
camPos    = get(ax, 'CameraPosition');  % camera position
rot = rotationtoviewingplane(ax);       % Get the rotation matrix

% apply the rotation matrix to all coordpoints, 
rotcoord = rot * plotcoordlist';
rotpoint = rot * point' ;

% loop over slices
faceindex = nan(ndata,1);
for sl = 1:ndata
    % find nearest point(s) (using only points in slice sl)
    coordind = dsearchn(rotcoord(1:2, (1:nbpoints(sl))+sum(nbpoints(1:sl-1)) )', rotpoint(1:2));
    for ci=1:numel(coordind) % (there may be several, if distances equal)
        % test if the point is in the polygon
        rotpoly = rot * plotvertexlist(sum(nbverts(1:sl-1))+datalist{sl}.faces(coordind(ci), isfinite(datalist{sl}.faces(coordind(ci),:))), :)';
        if inpolygon(rotpoint(1),rotpoint(2),rotpoly(1,:)', rotpoly(2,:)')
            faceindex(sl) = coordind(ci);  % if yes, done
        else % otherwise, check neighbors (because it's not a strict Voronoi-diagram in transformed coord's, that's possible)
            if isfield(datalist{sl},'delaunaytri')
                neighbors =  findneighbors(coordind, datalist{sl}.coordlist, datalist{sl}.delaunaytri);
            else
                neighbors =  findneighbors(coordind, datalist{sl}.coordlist);
            end
            for ni = 1:numel(neighbors) % do the same as before for ci
                rotpoly = rot * plotvertexlist(sum(nbverts(1:sl-1))+datalist{sl}.faces(neighbors(ni), isfinite(datalist{sl}.faces(neighbors(ni),:))), :)';
                if inpolygon(rotpoint(1),rotpoint(2),rotpoly(1,:)', rotpoly(2,:)')
                    faceindex(sl) = neighbors(ni);  break;   % if yes, done
                end
            end
        end
    end
end

if all(isnan(faceindex)), return; end   % no intersection with any face. Stop.

% With faceindex(), we know for each slice if and through which face the viewing line crosses.
% Now it has to be determined which among these is nearest to the viewpoint.
        
if sum(isfinite(faceindex))==1  % only one slice -> done.
    nslice = find(isfinite(faceindex));
    nface = sum(nbpoints(1:nslice-1)) + faceindex(nslice);
    nfaceonslice = faceindex(nslice);
    return;
end

% Otherwise: determine the points where each face is intersected by the
% viewline and choose nearest one to camera

dist=nan(1,ndata);
for sl = find(isfinite(faceindex'))
    faceverts = plotvertexlist(sum(nbverts(1:sl-1))+datalist{sl}.faces(faceindex(sl), isfinite(datalist{sl}.faces(faceindex(sl),:))), :);
    [normal,c] = fithyperplane(faceverts);  % obtain a plane to describe this face
    lambda = (c-point(1,:)*normal)/((point(2,:)-point(1,:))*normal); % determine intersection place--viewline
    spoint = point(1,:) + lambda*(point(2,:)-point(1,:));
    dist(sl) = norm(camPos-spoint); %#ok<AGROW>
end

[~,nslice] = min(dist);
nface = sum(nbpoints(1:nslice-1)) + faceindex(nslice);
nfaceonslice = faceindex(nslice);

        
        
        
        
        