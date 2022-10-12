function [vertexlist,faces,colors,area]=calcpatch(xdat,ydat,cdat,xsize,ysize,ratio)

%Calculates the optimum cells for representation of the data at xdat,ydat
%as a coloured patch.
%
%To achieve an optimum representation in case of overlapping cells and
%non-equally spaced cell centers, the patch is based on a Voronoi diagram.
%If necessary, the Voronoi cells are cut to not exceed the maximum cell
%size. Note that in most standard cases (regular mesh xdat,ydat), this
%algorithm results in the usual rectangular patch.
%
%xdat, ydat:    cell center coordinates (measured points) in arbitrary format and order
%xsize,ysize:   maximum (rectangular) extension of the cell
%ratio:         how to scale the coordinate axis - necessary if different
%               dimensions (energy, angle) to obtain a nice Voronoi diagram
%
%vertexlist:    rows are the x,y coordinates of all vertices
%faces:         rows are lists of the vertices (in vertexlist) defining each face
%colors:        colors of the faces
%area:          surface area of the cells
%
%to plot, use e.g. patch('Faces',faces,'Vertices',vertexlist,'FaceColor','flat','FaceVertexCData',colors)
%
%Paul Steffens, 07/2008
%
%

xdat=xdat(:);
ydat=ydat(:);
cdat=cdat(:);


xs=max(xdat)-min(xdat); ys=max(ydat)-min(ydat);
xext=[mean(xdat);           max(xdat)+3*xsize+xs; mean(xdat);            min(xdat)-3*xsize-xs];
yext=[max(ydat)+3*ysize+ys; mean(ydat);           min(ydat)-3*ysize-ys;  mean(ydat)          ];
%add artificial points to make all cells finite


if nargin<6, ratio = [1,1]; end

%compute Voronoi vertices V and cells C

[V,C] = voronoin([ ratio(1)*[xdat;xext], ratio(2)*[ydat;yext] ],{'Qbb','Qz'});
V(:,1) = V(:,1)/ratio(1);
V(:,2) = V(:,2)/ratio(2);


% count vertices to preallocate memory for vertexlist
nvert = 0;
for i=1:length(xdat)
    nvert = nvert + numel(C{i});
end    
vertexlist=zeros(nvert,2);

faces = zeros(length(xdat),4);
area  = zeros(length(xdat),1);

nvert = 1;
for i=1:length(xdat)
    Qu.x=[xdat(i)-xsize/2, xdat(i)+xsize/2, xdat(i)+xsize/2, xdat(i)-xsize/2];
    Qu.y=[ydat(i)-ysize/2, ydat(i)-ysize/2, ydat(i)+ysize/2, ydat(i)+ysize/2];
    Qu.hole=0;
    Vo.x=V(C{i},1);
    Vo.y=V(C{i},2);
    Vo.hole=0;
    %Qu is a rectangle around the central point ("maximum cell size"), Vo is the Voronoi cell
    %Intersect these two polygons, if necessary
    if (max(Vo.x) > xdat(i) + xsize/2) || (min(Vo.x)  < xdat(i)-xsize/2 ) || (max(Vo.y) > ydat(i) + ysize/2) || (min(Vo.y)  < ydat(i)-ysize/2 )
       P=cellclip(Qu,Vo);
%      P=PolygonClip(Qu,Vo,1);  % This is significantly more efficient than the previous line, but works only in Windows       
    else
        P = Vo;
    end
    
    vertexlist( nvert : (nvert+numel(P.x)-1), 1:2)  = [P.x(:), P.y(:)];
    nvert = nvert + numel(P.x);
    
    
    faces(i,1:length(P.x))= (nvert-length(P.x)) : (nvert-1);
    
    % Calculate the area of this cell
    if nargout > 3
        area(i) = 1/2 * sum( (P.x(:)' - P.x([2:end,1])') .* (P.y(:)' + P.y([2:end,1])') );
    end
end

faces(faces==0)=NaN;

[vertexlist, faces] = optimizepatch(vertexlist, faces);

colors=cdat;
