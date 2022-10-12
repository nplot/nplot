function exportpov


%% Subfunction coordinate transformation
    function co = coord(co)
        co = [co(1)/da(1), co(3)/da(3), co(2)/da(2)];
        % scale with aspect ratio and exchange y,z
    end


%% ask for filename
directories = getoption('directories');
if ~ isempty(directories), stddir = directories.povoutput; end
[file,path] = uiputfile('*.pov','Save .pov Skript', stddir); 
if file==0, return; end

%% get plot window

[d,plotstruct] = getfiguredata;
da = daspect;


if isfield(plotstruct,'sectionlist')
    section = plotstruct.sectionlist;
    nsecs = numel(unique(plotstruct.sectionlist));
else
    [nsecs,section] = finddistinctsections(plotstruct.datalist.faces);
    plotstruct.sectionlist = section;
    guidata(gcf,plotstruct);
end

clims = caxis;

%% some defs
smallval = .01;
textscale  = 0.18;  % Text size
textalignout = -1;  % +1 or -1,  Text inside or outside plot box ?
lineradius = 0.01;  % Radius of cyclindric elements that make the lines

%% write the file

file = fopen([path,file],'w');

try  %% in case of error, close output stream correctly


fprintf(file,'#include "colors.inc"\n#include "textures.inc"\n#include "transforms.inc"\n\n');

% Material definitions
fprintf(file,'#declare Side_Finish=\n  finish {\n    specular 1\n     roughness 0.001\n     ambient 0\n     diffuse 0\n     reflection 0.05\n}\n\n');
fprintf(file,'#declare SideMaterial =\n   texture {\n      pigment { rgbf<1.0, 1.0, 1.0, 0.8> }\n       finish  { Side_Finish }\n}\n\n');


%Background color
fprintf(file,'background { White }\n\n');

% camera
fprintf(file,'camera {\n   location <%f,%f,%f>\n   look_at <%f,%f,%f>\n\n   angle %f\n}\n',coord(campos),coord(camtarget),2*camva);

% light
fprintf(file,'light_source { <%f,%f,%f>\n   color White\n   //area_light <15, 0, 0>, <0, 0, 15>, 5, 5\n   //adaptive 1\n}\n',coord(campos+[0,0,20]));
fprintf(file,'light_source { <%f,%f,%f>\n   color Gray50\n   //area_light <15, 0, 0>, <0, 0, 15>, 5, 5\n   //adaptive 1\n}\n',coord(campos+[0,0,0]));


% Koord-system, Achsen und Ebenen
xl=xlim; yl=ylim; zl=zlim; 
viewdir = (camtarget - campos)>0; % 1 for each axis if looking in positive direction
xback = xl(viewdir(1)+1);  yback = yl(viewdir(2)+1);  zback = zl(viewdir(3)+1); % far end of axes from viewpoint. Draw the planes here
xfront = xl(2-viewdir(1)); yfront = yl(2-viewdir(2)); zfront = zl(2-viewdir(3));% near end... Draw the axes labels

fprintf(file,'polygon {\n    5, <%f,%f,%f> <%f,%f,%f> <%f,%f,%f> <%f,%f,%f> <%f,%f,%f> texture{SideMaterial}}\n\n', ...
    coord([xback,yl(1),zl(1)]), coord([xback,yl(1),zl(2)]), coord([xback,yl(2),zl(2)]), coord([xback,yl(2),zl(1)]), coord([xback,yl(1),zl(1)]));
fprintf(file,'polygon {\n    5, <%f,%f,%f> <%f,%f,%f> <%f,%f,%f> <%f,%f,%f> <%f,%f,%f> texture{SideMaterial}}\n\n', ...
    coord([xl(1),yback,zl(1)]), coord([xl(1),yback,zl(2)]), coord([xl(2),yback,zl(2)]), coord([xl(2),yback,zl(1)]), coord([xl(1),yback,zl(1)]));
fprintf(file,'polygon {\n    5, <%f,%f,%f> <%f,%f,%f> <%f,%f,%f> <%f,%f,%f> <%f,%f,%f> texture{SideMaterial}}\n\n', ...
    coord([xl(1),yl(1),zback]), coord([xl(1),yl(2),zback]), coord([xl(2),yl(2),zback]), coord([xl(2),yl(1),zback]), coord([xl(1),yl(1),zback]));
crad = min(abs([xl(2)-xl(1),yl(2)-yl(1),zl(2)-zl(1)]))/300;
% Koord-linien
rad = min([1/50, min(abs([xl(2)-xl(1),yl(2)-yl(1),zl(2)-zl(1)]))/200 ])  / 3;
tickparam = 'texture{pigment { rgbt<.3, .3, .3, 0.1> }}';
stepopt = [.01,.02,.05,.1,.2,.5,1,2,5];
for x = get(gca,'xtick')
    [step,i] = min(abs(stepopt / ((yl(2)-yl(1))/30) - 1)); step = stepopt(i); % which of the poss. step sizes gives closest to 30 points?
    for y=ceil(yl(1)/step)*step : step : floor(yl(2)/step)*step
        fprintf(file,'cylinder { <%f,%f,%f>,  <%f,%f,%f>, %f %s}\n', coord([x,y,zback-smallval]),coord([x,y,zback+smallval]),rad,tickparam);
    end
    [step,i] = min(abs(stepopt / ((zl(2)-zl(1))/30) - 1)); step = stepopt(i);
    for z=ceil(zl(1)/step)*step : step : floor(zl(2)/step)*step
        fprintf(file,'cylinder { <%f,%f,%f>, <%f,%f,%f>, %f %s}\n', coord([x,yback-smallval,z]),coord([x,yback+smallval,z]),rad,tickparam);
    end
    textstr = num2str([x,smallval,textscale], 'ttf "verdana.ttf", "%g", %f, 0  scale %f rotate<90,0,0>');
    fprintf(file,'text {%s texture{pigment{color Gray20}} Align_Trans( text{%s}, <-1,1,%d>, <%f,%f,%f>)}\n', textstr, textstr, (2*viewdir(2)-1)*textalignout, coord([x,yfront,zback]));

end
for y = get(gca,'ytick')
    [step,i] = min(abs(stepopt / ((xl(2)-xl(1))/30) - 1)); step = stepopt(i);
    for x=ceil(xl(1)/step)*step : step : floor(xl(2)/step)*step
        fprintf(file,'cylinder { <%f,%f,%f>, <%f,%f,%f>, %f %s}\n', coord([x,y,zback-smallval]),coord([x,y,zback+smallval]),rad,tickparam);
    end
    [step,i] = min(abs(stepopt / ((zl(2)-zl(1))/30) - 1)); step = stepopt(i);
    for z=ceil(zl(1)/step)*step : step : floor(zl(2)/step)*step
        fprintf(file,'cylinder { <%f,%f,%f>, <%f,%f,%f>, %f %s}\n', coord([xback-smallval,y,z]),coord([xback+smallval,y,z]),rad,tickparam);
    end
    textstr = num2str([y,smallval,textscale], 'ttf "verdana.ttf", "%g", %f, 0  scale %f rotate<90,0,0>');
    fprintf(file,'text {%s texture{pigment{color Gray20}} Align_Trans( text{%s}, <%d,1,-1>, <%f,%f,%f>)}\n', textstr, textstr, (2*viewdir(1)-1)*textalignout, coord([xfront,y,zback]));
end
for z = get(gca,'ztick')
    [step,i] = min(abs(stepopt / ((xl(2)-xl(1))/30) - 1)); step = stepopt(i);
    for x=ceil(xl(1)/step)*step : step : floor(xl(2)/step)*step
        fprintf(file,'cylinder { <%f,%f,%f>, <%f,%f,%f>, %f %s}\n', coord([x,yback-smallval,z]),coord([x,yback+smallval,z]),rad,tickparam);
    end
    [step,i] = min(abs(stepopt / ((yl(2)-yl(1))/30) - 1)); step = stepopt(i);
    for y=ceil(yl(1)/step)*step : step : floor(yl(2)/step)*step
        fprintf(file,'cylinder { <%f,%f,%f>, <%f,%f,%f>, %f %s}\n', coord([xback-smallval,y,z]),coord([xback+smallval,y,z]),rad,tickparam);
    end
    textstr = num2str([z,smallval,textscale], 'ttf "verdana.ttf", "%g", %f, 0  scale %f ');
    fprintf(file,'text {%s texture{pigment{color Gray20}} Align_Trans( text{%s}, <%d,-1,1>, <%f,%f,%f>)}\n', textstr, textstr, (2*viewdir(1)-1)*textalignout, coord([xfront,yback,z]));
% fprintf(file,'text {ttf "verdana.ttf", "%g", %f, 0  scale .2 texture{pigment{color Gray20}} translate <%f,%f,%f>}\n', z, smallval, coord([xfront,yback,z]));
end   

% Achsen (Kanten)
fprintf(file,'cylinder {<%f,%f,%f>, <%f,%f,%f>, %f texture {Chrome_Metal}}\n\n', coord([xback,yback,zl(1)]), coord([xback,yback,zl(2)]), crad);
fprintf(file,'cylinder {<%f,%f,%f>, <%f,%f,%f>, %f texture {Chrome_Metal}}\n\n', coord([xback,yl(1),zback]), coord([xback,yl(2),zback]), crad);
fprintf(file,'cylinder {<%f,%f,%f>, <%f,%f,%f>, %f texture {Chrome_Metal}}\n\n', coord([xl(1),yback,zback]), coord([xl(2),yback,zback]), crad);


% Create nearly transparent planes just below/above every slice with
% constant z-value
 for nsl=1:nsecs
     zval = plotstruct.coordlist(section==nsl,3);
     if max(zval)-min(zval) > .1, continue; end
     zval = mean(zval) + (viewdir(3)-.5)*smallval;
     fprintf(file,'polygon {\n    5, <%f,%f,%f> <%f,%f,%f> <%f,%f,%f> <%f,%f,%f> <%f,%f,%f> texture{pigment { rgbft<1.0, 1.0, 1.0, .5, 0.8> }}}\n\n', ...
          coord([xl(1),yl(1),zval]), coord([xl(1),yl(2),zval]), coord([xl(2),yl(2),zval]), coord([xl(2),yl(1),zval]), coord([xl(1),yl(1),zval]));
 end

% determine colors (index into colormap)
if strcmpi(plotstruct.linlog,'LIN')
    cvalue = plotstruct.datalist.valuelist(:,1);
elseif strcmpi(plotstruct.linlog,'LOG')
    m = min(plotstruct.datalist.valuelist(:,1));
    cvalue = log10(plotstruct.datalist.valuelist(:,1) + max(0,-m) + 1);
end
cmap = colormap; numcol = size(cmap,1); 
cl = caxis;
cvalue = round ((cvalue-cl(1)) / (cl(2)-cl(1)) * (numcol-1) );
cvalue(cvalue>=numcol) = numcol-1;
cvalue(cvalue<0) = 0;


% Colored surfaces (as one object)
fprintf(file,'mesh2 {\n   vertex_vectors {\n      %d,\n',size(plotstruct.vertexlist,1));
% list of nodes
for i=1:size(plotstruct.vertexlist,1)
    fprintf(file,'      <%f,%f,%f>,\n',coord(plotstruct.vertexlist(i,1:3)));  
end
% list of colors
fprintf(file,'   }\n   texture_list {\n      %d,\n',numcol);
for i=1:numcol
    fprintf(file,'      texture {pigment { rgbt <%f,%f,%f,0.15> }  finish{diffuse 1 specular .8 roughness .1 } },\n',cmap(i,1:3)); 
end
% list of faces
fprintf(file,'   }\n   face_indices {\n      %d,\n', sum(sum(isfinite(plotstruct.datalist.faces(:,3:end)))));
for i=1:size(plotstruct.datalist.faces,1)
    for j=3:sum(isfinite(plotstruct.datalist.faces(i,:)))
        fprintf(file,'      <%d,%d,%d>, %d,\n',plotstruct.datalist.faces(i,[1,j-1,j])-1, cvalue(i) );       % "-1" because povray-indices 0-based
    end
end
fprintf(file,'   }\n}\n ');

%% Lines (HKL, powder etc.)
if isfield(plotstruct,'linedef')
    for nline = 1:length(plotstruct.linedef)
        switch readinput('edgecolor',plotstruct.linedef{nline}.properties)
            case 'k', colorstr = 'Black';
            case 'r', colorstr = 'Red';
            case 'g', colorstr = 'Green';
            case 'b', colorstr = 'Blue';
            case 'w', colorstr = 'White';
            otherwise colorstr = 'Gray60'; %#ok<SEPEX>
        end   
        for n = 1:size(plotstruct.linedef{nline}.connection,1)
            fprintf(file,'cylinder{ <%f,%f,%f>, <%f,%f,%f,>, %f texture{pigment{color %s}}}\n', ...
                coord(plotstruct.linedef{nline}.points(plotstruct.linedef{nline}.connection(n,1),1:3)), ...
                coord(plotstruct.linedef{nline}.points(plotstruct.linedef{nline}.connection(n,2),1:3)), lineradius, colorstr ) ;
        end
        fprintf(file,'\n');
    end
end

%% Plot evtl. planes (volume cuts etc.)
if isfield(plotstruct, 'planedef')
    if ~iscell(plotstruct.planedef), planedef{1}=plotstruct.planedef; else planedef = plotstruct.planedef; end  % ensure cell array
    for npl = 1:length(planedef)
        corners = [xl(1),yl(1),zl(1);xl(1),yl(1),zl(2);xl(1),yl(2),zl(1);xl(1),yl(2),zl(2);xl(2),yl(1),zl(1);xl(2),yl(1),zl(2);xl(2),yl(2),zl(1);xl(2),yl(2),zl(2)];
        [planevectors, origin] = getplaneparameter(planedef{npl}.normal(:)', planedef{npl}.c);
        % intersection
        [cutvertices, cutorder] = createmesh (corners, 1:8, planevectors', origin);
        % transform projected ccordinates back to full system
        cutvertices = cutvertices * planevectors' + repmat(origin(:)',size(cutvertices,1),1);
        
%         patch('Vertices',cutvertices, 'Faces',cutorder, planedef{npl}.properties{:});
        
        fprintf(file,'polygon {\n   %d   ',size(cutvertices,1));
        for i=1:size(cutvertices,1), fprintf(file,', <%f,%f,%f>', coord(cutvertices(cutorder(i),1:3))); end
        fprintf(file,'\n   texture {pigment{rgbft <.8,0,0,.2,.8>}}\n}\n');
    end
end

catch 
    fclose(file);
    rethrow(lasterror);
    return;
end

fclose(file);

end