function valsum = integratepatch(faces, vertices, values, errors, corner, sides, steps, avgopt)

% Performs integration along different directions (projection) of a defined region as a weighted average.
% The intersection of the patch cells with the slices to be integrated is calculated in each step. 
% - faces, vertices, values, errors define the data
% - corner defines one corner of the parallelepiped to be considered
% - sides give N edges as vectors
% - steps (N components) defines the number of slices to be taken along each direction.
%
% Prepared to work in arbitrary dimensions.

% P. Steffens 08/2008


%sides = [side1(:)'; side2(:)'];
valsum = [];
if det(sides)==0, fprintf('Collinear vectors! Exit integratepatch.\n'); return; end
ndims = size(sides,1);
if ndims ~= size(vertices,2) || ndims~=numel(steps), fprintf('Dimensional mismatch in integratepatch.\n'); return; end

c=[];
count1 = 1;

corner = corner(:)';
values = values(:);
errors = errors(:);
if nargin < 8 , avgopt=1; end

%******************************************
%helper function
function cutaway
        % Obtain plane normal vector and C
        if any(sum(abs(points),1) <= 1E-10)
            normal=zeros(ndims,1); 
            normal(sum(abs(points),1) <= 1E-10)=1; 
        else
            normal = null(round(points * 1E6));  
            if ~isempty (normal) && any(normal(1:ndims)~=0)
                normal = normal(:,1);
            else
                normal = sum(inv(points),2);  
            end
        end
        normal = normal ./ sqrt(normal'*normal);
        C = points(1,:)*normal;

        % use this hyperplane to cut
        [newfaces, newvertices, cellassignment, partition] = cutpatch(faces, vertices, normal, C);
        % Which is the good side?
        goodside = (sign(sides(ns,:) * normal) > 0) ~= lowerupper;
        
        % First, integrate bad cells 
        outvalues = values(cellassignment(partition~=goodside),1);
        outerrors = errors(cellassignment(partition~=goodside),1);
        outfaces = newfaces(partition~=goodside,:);
        
        if calcarea
%            ar = abs(polygonarea(newvertices(:,1),newvertices(:,2), outfaces));
            ar = polyedervolume(newvertices, outfaces);     % ** Here, one could save time by not recalculating volume of unchanged cells!!
            if avgopt ==1                                                                       % ** !
                [wmean,werr] = weightedmean(outvalues, outerrors, ar .* abs(outvalues) ./ outerrors.^2);
            else
                [wmean,werr] = weightedmean(outvalues, outerrors, ar ./ outerrors.^2);
            end
            c = [c, size(outfaces,1)];
        end
        
        % Keep only good cells
        values = values(cellassignment(partition==goodside),1);
        errors = errors(cellassignment(partition==goodside),1);
        faces = newfaces(partition==goodside,:);
        % Try to delete unused vertices...
        % ... but do this not every time, but only every fifth (better performance)
        if count1 == 5
            usedvert = faces(isfinite(faces));
            [usedvert,m,index]=unique(usedvert); 
            vertices = newvertices(usedvert,:);
            faces(isfinite(faces)) = index;
            count1 = 1;
        else
            vertices = newvertices;
            count1 = count1 + 1;
        end
end
%******************************************


% First, cut away all outlying parts of the patch
calcarea = false;
for ns = 1:ndims;
    for lowerupper = 0:1 % 0 lower, 1 upper boundary
        % Setup hyperplane...
        % Setup matrix of ndims points defining the plane
        points = repmat(corner,ndims,1) + lowerupper * repmat(sides(ns,:),ndims,1) + repmat((1:ndims ~= ns)',1,ndims) .* sides;
        % Cut away everything beyond this hyperplane
        cutaway;
    end
end

% Keep this patch;
% (everything outside the parallelepiped has now been cut away)
faces0      = faces;
vertices0   = vertices;
values0     = values;
errors0     = errors;

lowerupper = 0;
valsum = cell(ndims,1);
calcarea = true;

for ns = 1:ndims
    % Restore patch
    faces      = faces0;
    vertices   = vertices0;
    values     = values0;
    errors     = errors0;
    
    % Now, cut stepwise
    for st=1:steps(ns)
        % Setup hyperplane...
        % Setup matrix of ndims points defining the plane
        % (use tiny offset to avoid doing exactly the same cut twice)
        points = repmat(corner,ndims,1) + (st+1E-8)/steps(ns) * repmat(sides(ns,:),ndims,1) + repmat((1:ndims ~= ns)',1,ndims) .* sides;
        % Cut portion below first step
        cutaway;
        valsum{ns} = [valsum{ns}; wmean, werr];
    end  
end

end
