function decorateQplot(ax,bx,lattice,qv,type,axhandle,opt)

% Plots :
% a) the orienting vectors, as given in the scan (more precisely: their projection onto the scattering plane)
% b) a grid consisting of these vectors
% c) constant-psi lines
% into the coordinate system (for qx,qy-axes only!)
%
% type: "arrows" for (a), "grid" for (b), "psi=90" (or other val.) for (c)
% ax, bx:  orientation vectors (HKL)
% lattice: structure containing lattice constants etc. (use getlattice)
% qv: vertical component of Q
% opt=true fixes the range of the axes

% P. Steffens 07/2008

%%
function setequallimits(xlimits,ylimits)
    range = max( abs( [ xlimits(2)-xlimits(1), ylimits(2)-ylimits(1) ]));
    xlimits = (xlimits(2)+xlimits(1))/2 + [-range/2, range/2];
    ylimits = (ylimits(2)+ylimits(1))/2 + [-range/2, range/2];
    xlim(xlimits);  
    ylim(ylimits); 
end

%%

if nargin<6
    axhandle=gca;
end

if nargin<5
    type = 'arrows';
end

if nargin>6 
    keepview = opt;
else
    keepview = false;
end

font='Verdana';

holdstate=ishold(axhandle);

UB = UBmatrix( lattice,ax,bx );

[q1x, q1y] = calcQS(ax(1), ax(2), ax(3), UB);
[q2x, q2y] = calcQS(bx(1), bx(2), bx(3), UB);


%Consider eventual offset of center in non-orthogonal lattices (**to be checked!**)
cx = [ax(2)*bx(3)-ax(3)*bx(2), ax(3)*bx(1)-ax(1)*bx(3), ax(1)*bx(2)-ax(2)*bx(1)];
%third (perpendicular) vector  in HKL coord.
[q3x,q3y,q3z] = calcQS(cx(1), cx(2), cx(3), UB);
centerx =  qv/q3z * q3x;    % ** '-'
centery =  qv/q3z * q3y;    % ** '-'
c=[centerx,centery];


hold(axhandle,'on');

xlimits = get(axhandle,'xlim');
ylimits = get(axhandle,'ylim');
    
if strcmpi(type,'ARROWS')
    quiver(axhandle, centerx+[0,0], centery+[0,0], [q1x,q2x], [q1y,q2y], 0, 'k');
    text(c(1)+q2x+.1*q1x,c(2)+q2y,['[' num2str(bx(1)) ',' num2str(bx(2)) ',' num2str(bx(3)) ']' ], 'Fontname',font);    
    if sign(q1x)==-1; al='right'; else al='left'; end
    text(c(1)+1.1*q1x,c(2)+q1y,['[' num2str(ax(1)) ',' num2str(ax(2)) ',' num2str(ax(3)) ']' ], 'Fontname',font,'horizontalalignment',al);
    xlimits = [min(xlimits(1),-.2), max(xlimits(2),q1x+.2)];
    ylimits = [min(ylimits(1),-.2), max(ylimits(2),q2y+.2)];   
    if ~keepview
        setequallimits(xlimits, ylimits);
    end
elseif strcmpi(type,'GRID')
%     xlimits2 = [ floor(min(plotstruct.vertexlist(:,1))*2)/2, ceil(max(plotstruct.vertexlist(:,1))*2)/2] ; % "Normal" axis limits (full range)
%     ylimits2 = [ floor(min(plotstruct.vertexlist(:,2))*2)/2, ceil(max(plotstruct.vertexlist(:,2))*2)/2] ; % use these for plotting to ensure that grid covers full range (even when zooming out)
    xlimits2 = [-8,8];
    ylimits2 = [-8,8];
    for yi=ceil((ylimits2(1)-c(2))/abs(q2y)):floor((ylimits2(2)-c(2))/abs(q2y))
        line([xlimits2(1),xlimits2(2)],yi*abs(q2y)*[1,1]+c(2),'LineStyle',':','Color','r');
    end
    if q2x<0; q2x=-q2x; q2y=-q2y; end
    for xi = ceil((xlimits2(1)-c(1) - q2x*max(ylimits2/q2y)) / abs(q1x)) : floor((xlimits2(2)-c(1) + q2x*abs(min(ylimits2/q2y))) / abs(q1x))
        q1x=abs(q1x);
        x1 = max(xlimits2(1), q2x/q2y*ylimits2(1) + xi*q1x);
        x2 = min(xlimits2(2), q2x/q2y*ylimits2(2) + xi*q1x);
        if q2x < 1E-10
            y1 = ylimits2(1); y2 = ylimits2(2);
        else
            y1 = max(ylimits2(1), q2y/q2x*(x1 - xi*q1x));
            y2 = min(ylimits2(2), q2y/q2x*(x2 - xi*q1x));
        end
        line([x1,x2]+c(1),[y1,y2]+c(2),'LineStyle',':','Color','r');
    end
    if ~keepview
        setequallimits(xlimits, ylimits);
    end
elseif strcmpi(type(1:3),'PSI')
    psival=str2num(type(5:end));
    [lx,ly]=inplaneQ(0:100,psival*ones(1,101),ki,kf,qv);
    plot(lx,ly,'-r');
    text(lx(end),ly(end),['\psi=' num2str(psival)]);
end


    
if holdstate==0
    hold(axhandle,'off');
end


end %function