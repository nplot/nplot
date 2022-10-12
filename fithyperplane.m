function [normal,c] = fithyperplane ( points, orthvecs )

% Fit hyperplane to point set points by minimizing the sum of squared distances
% Works in N dimensions
% [normal,c] = fithyperplane ( points, orthvec )
% points:          (NP x dim)
% orthvec(opt.):   vectors (nvec * dim) to which normal must be orthogonal

% P. Steffens, 03/2009

[npoints, ndims] = size(points);

M = zeros(ndims);
b = zeros(ndims,1);

for i=1:ndims
    for j=1:ndims
        M(i,j) = sum(points(:,i).*points(:,j));
    end
    b(i)= sum(points(:,i));
end

% if nargin>1
%     % if orthvec(s) given, reduce dimension of problem, modify M
%     useddims = true(1,ndims);   
%     for nov=1:size(orthvec,1)
%         [m,nnind(nov)] = max(abs(orthvec(nov,useddims))); %nnind(nov) is non-null entry of orthvec(nov)
%         useddims(nnind(nov)) = false;
%         for sp = find(useddims)
%             M(:,sp) = M(:,sp) - orthvec(nov, sp) / orthvec(nov, nnind) * M(:,nnind);
%         end
%     end
%     M = M(useddims,useddims);
%     b = b(useddims);
% end

% If orthvec(s) given, introduce Lagrange multipliers, one for each
% condition orthvec(i,:)*normal = 0

if nargin<2, orthvecs = []; end
M = [M, orthvecs'; orthvecs, zeros(size(orthvecs,1))];
b = [b; zeros(size(orthvecs,1),1)];



normal = null(round(M * 1E6));  % See if solution for c=0 exists (round sets accuracy)
if ~isempty (normal) && any(normal(1:ndims)~=0)
    normal = normal(1:ndims,1);
    c = 0;
else
    normal = inv(M) * b;
    normal = normal(1:ndims);
    c = 1/sqrt(sum(normal.^2));
    normal = normal*c;    
end


% % Now, treat the possibility that the input forms a point set that has not sufficient dimension
% % This shows up as a nearly singular matrix M
% deletedim = false(ndims,1);
% while cond(M) > 1E8
%     % determine dimension mindim which to delete
%     % set mindim-th component of normal to 1
%     mindim = find(all(M==0,1),1);
%     if isempty(mindim)
%         mincond = inf;
%         mindim = 1;
%         for i=1:size(M,1)
%             icond = cond(M([1:(i-1),(i+1):end],[1:(i-1),(i+1):end]));
%             if  icond < mincond
%                 mincond = icond;
%                 mindim = i;
%             end
%         end
%     end
%     b = b - M(:,mindim);
%     b = b([1:(mindim-1),(mindim+1):end]);
%     M = M ( [1:(mindim-1),(mindim+1):end], [1:(mindim-1),(mindim+1):end] );
%     mindim = mindim + sum(deletedim(1:mindim)); % (original dimension index) 
%     deletedim(mindim) = true;
% end            
%         
% 
% normal(~deletedim) = inv(M) * b;
% normal(deletedim) = 1;
% c = 1/sqrt(sum(normal.^2));
% normal = normal*c;
