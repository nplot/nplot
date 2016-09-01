function P = cellclip(Qu, Vo)

% Spezialized version of Polygonclip
% (assumes Qu as a rectangle)
%
% P. Steffens, 04/2008

faces = 1:numel(Vo.x);
vertices = [Vo.x(:), Vo.y(:)];

center = [mean(Qu.x); mean(Qu.y)];

for i=1:4
    dy = Qu.y(mod(i,4)+1) - Qu.y(i);
    dx = Qu.x(mod(i,4)+1) - Qu.x(i);
    n = [ -dy, dx ];
    n = n / sqrt(n*n');
    c = n * [Qu.x(i) ; Qu.y(i)];
    
    part0 = ( n * center > c);
    
    [newfaces, newvertices, cellassign, part] = cutpatch(faces, vertices, n, c);
    
    vlist = newfaces(logical(part==part0), :);  % These are the numbers of the vertices of the good polygon
    vlist = vlist(isfinite(vlist));
    
    vertices = newvertices( vlist , : );
    
    faces = 1:size(vertices,1);
    
end

P.x = vertices(:,1);
P.y = vertices(:,2);