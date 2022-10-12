function np = estimateNP (points, vertices, faces)

% Estimate number of scan points that is reasonable for the scan from points(1,:) to points(2,:)
% 
% P. Steffens, 01/2012


% Cut away portions of patch beyond start and end point
% follow a procedure similar to the one in 'cutaway' in 'integratepatch'
% Therefore, cut through planes perpendicular to the line that pass through
% start and end point.
normal = points(2,:) - points(1,:);
if max(abs(normal))<1E-4, np=1; return; end
normal = normal' ./ sqrt(normal*normal');

C = points(1,:)*normal;

[newfaces, newvertices, cellassignment, partition] = cutpatch(faces, vertices, normal, C);
newfaces = newfaces(partition==1,:);

C = points(2,:)*normal;

[newfaces, newvertices, cellassignment, partition] = cutpatch(newfaces, newvertices, normal, C);
newfaces = newfaces(partition==0,:);

% Now, cut along scan line and count number of faces before and after
nf_before = size(newfaces,1);

normal = [-normal(2); normal(1); normal(3:end)];
normal = normal ./ sqrt(normal'*normal);
C = points(1,:)*normal;

newfaces = cutpatch(newfaces, newvertices, normal, C);  % cut
nf_after = size(newfaces,1);  %count again. These are more, because each face through which the line runs is cut into two.

np = max(1, nf_after - nf_before);

