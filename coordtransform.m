function [newlist, optout] = coordtransform(datalist, coordtype, varargin)

% newlist = coordtransform(datalist, coordtype, varargin)
% Transform to a different coordinate system
%
% optout set according to case (pairs name/value)
%
% coordtype:    New coordinate type
% if coordtype='PROJECTION': do a projection on a subspace:
%   first argument in varargin are vectors (n rows are n vectors)
%   second argument in varargin is origin, third (optional) becomes newlist.coordtype
%
% if coordtype==datalist.coordtype --> do nothing
% if transformation not explicitly implemented: 
%   can use option 'knowncoord' to interpolate 

% P. Steffens, 01/2012


    function trygeneraltransform
        % if target coordinates of some points are known, try to find the
        % others by simple linear method. Good approximation for points
        % lying close to known ones, and far from singularities
        % Typical case: back-transformation in coordinates defined by one
        % scan, but unknown for others
        knowncoord = [];
        if strcmp(coordfield,'coordlist'),      knowncoord = readinput('knowncoord',varargin); 
        elseif strcmp(coordfield,'vertexlist'), knowncoord = readinput('knownvertices',varargin); end
        if ~isempty(knowncoord)
            try
                calcind = any(isnan(knowncoord),2); % Index of points to be calculated
                for cd=1:size(knowncoord,2)
                    newlist.(coordfield)(calcind,cd) = griddatan(datalist.(coordfield)(~calcind,:), knowncoord(~calcind,cd), datalist.(coordfield)(calcind,:));
                end
                newlist.(coordfield)(~calcind,:) = knowncoord(~calcind,:);
            catch
                fprintf('Problem during calculation of coordinate transform.\n');
                newlist = [];
                return
            end
        else
            if strcmp(coordfield,'coordlist')
                fprintf('Transformation not implemented or impossible or coordinate type not recognized.\n');
                newlist = [];
                return
            else
                if hasfield(newlist,coordfield), newlist = rmfield(newlist,coordfield); end
                % for supplementary fields: simply delete them if they cannot be transformed
            end
        end
    end


optout = {};
newlist = datalist;
if isempty(datalist) || (nargin<3 && strcmpi(datalist.coordtype,coordtype))
    return;
end
newlist.coordtype = upper(coordtype);

[mass_n, meVJ, hbar] = getoption('mass_n','meVJ','hbar');


if isfield(datalist,'vertexlist'), coordfieldlist = {'coordlist','vertexlist'}; else coordfieldlist = {'coordlist'}; end  
% These are the names of the fields that contain coordinates to be transformed
% ** could do this more general, allowing more fields
% * (do always coordlist first!)


for coordfield = coordfieldlist
    
    coordfield = coordfield{1}; %#ok<FXSET>

    if strcmpi(coordtype,'QXY') || strcmpi(coordtype,'QPLANE')
        if strcmpi(datalist.coordtype,'ANGLES')
            % Transform angles into qx,qy
            [newlist.(coordfield)(:,1), newlist.(coordfield)(:,2)] = inplaneQ(datalist.(coordfield)(:,1), datalist.(coordfield)(:,2), datalist.KI, datalist.KF, datalist.QVERT);
        else
            trygeneraltransform;
        end

    elseif strcmpi(coordtype,'QXYZ')
        if strcmpi(datalist.coordtype,'ANGLESQZ')
            % Transform 3d (angles,qz) into qx,qy,qz
            [newlist.(coordfield)(:,1), newlist.(coordfield)(:,2)] = inplaneQ(datalist.(coordfield)(:,1), datalist.(coordfield)(:,2), datalist.KI, datalist.KF, datalist.(coordfield)(:,3));
            newlist.(coordfield)(:,3) = datalist.(coordfield)(:,3);

        elseif any(strcmpi(datalist.coordtype, {'QXY','QPLANE','ANGLES'})) && isfield(datalist,'QVERT')
            % Transform 2d (qxy or angles) into 3d qx,qy,qz
            newlist = coordtransform(datalist,'QXY');
            newlist.(coordfield)(:,3) = datalist.QVERT;

        else
            trygeneraltransform;
        end

    elseif strcmpi(coordtype,'hklvectors')
        if strcmpi(datalist.coordtype,'ANGLESENERGY'), datalist = coordtransform(datalist,'QXQYEN'); end
        if strcmpi(datalist.coordtype,'ANGLESQZ'), datalist = coordtransform(datalist,'QXYZ'); end
        if strcmpi(datalist.coordtype,'ANGLES'), datalist = coordtransform(datalist,'QXY'); end
        if any(strcmpi(datalist.coordtype,{'QXY','QPLANE','QXYZ','QXQYEN'}))
           try
                if strcmpi(datalist.coordtype,'QXYZ') 
                    Qz = datalist.(coordfield)(:,3); 
                else
                    Qz = datalist.QVERT * ones(size(datalist.(coordfield),1),1);
                end
                sample = datalist.sampleinfo;  UB = UBmatrix(sample.lattice,sample.ax,sample.bx);
                % Determine basis vectors for plane
                perphkl = cross(sample.ax,sample.bx); 
                if all(perphkl~=0)
                    basis1 = sample.ax; % none of 100, 010 or 001 is in the plane, stick to ax as basis vector
                else
                    basis1 = [0,0,0]; basis1(find(perphkl==0,1,'first')) = 1; %choose first of 100,010,001 as first basis
                end
                basis1q = UB*basis1(:);
                % second basis vec must be orthogonal to first in Q(!)-space
                perpvec = cross(UB*sample.ax(:),UB*sample.bx(:));
                basis2 = cross(perpvec,basis1q); % in Q
                basis2 = UB\basis2;  % in HKL
                basis2 = round(basis2*1E10)*1E-10;  % precision
                basis2 = basis2' / basis2(find(abs(basis2)>.001,1,'first')); basis2(basis2==0)=0; % ensure first significant nonzero entry is 1
                basis2q = UB*basis2(:);
                
                perphkl = cross(basis1,basis2);
                perphkl = round(perphkl*1E10)*1E-10; perphkl = perphkl' / perphkl(find(abs(perphkl)>.001,1,'first'));  perphkl(perphkl==0)=0; % like above
                perphklq = UB*perphkl(:);
                p = Qz/perphklq(3);
                originq = p * perphklq';
                if any(max(originq(:,1:2),[],1) - min(originq(:,1:2),[],1) > 1e-4) % are the q-planes shifted (due to non-orth. axes)?
                    fprintf('HKL coordinates of different vertical momenta are shifted with respect to each other\n');
                    fprintf('(This is typically due to non-orthogonal axes.) Cannot determine well-defined transformation. Exit.\n');
                    newlist = []; return;
                end
                newlist.(coordfield)(:,1:2) = datalist.(coordfield)(:,1:2) - originq(:,1:2);
                if strcmpi(datalist.coordtype,'QXYZ'), newlist.(coordfield)(:,3) = p; 
                elseif strcmpi(datalist.coordtype,'QXQYEN'), newlist.(coordfield)(:,3) = datalist.(coordfield)(:,3); 
                end
                % projection on basis1 and basis2
                newlist.(coordfield)(:,1:2) = [newlist.(coordfield)(:,1:2)*basis1q(1:2)/(basis1q'*basis1q), newlist.(coordfield)(:,1:2)*basis2q(1:2)/(basis2q'*basis2q)];
                
                originhkl = round(UB\mean(originq,1)'*1E5)'/1E5;
                % axes labels etc.
                axname1 = {num2str(basis1,'units of [%g,%g,%g]'), num2str(originhkl,'(origin at [%g,%g,%g])' ) };
                axname2 = num2str(basis2,'units of [%g,%g,%g]');
                basisvectors.vector1 = basis1; basisvectors.type1 = 'hkl';
                basisvectors.vector2 = basis2; basisvectors.type2 = 'hkl';
                basisvectors.origin = originhkl;
                if strcmpi(datalist.coordtype,'QXYZ')
                    basisvectors.vector3 = basis3; basisvectors.type3 = 'hkl';
                    axname1 = {num2str(basis1,'units of [%g,%g,%g]'), '(origin at [0,0,0])' };
                    axname3 = num2str(basis3,'units of [%g,%g,%g]');
                elseif strcmpi(datalist.coordtype,'QXQYEN')
                    basisvectors.vector3 = [0,0,1]; basisvectors.type3 = 'energy';
                    axname3 = 'Energy';     
                end
                % assign output if not yet done    
                if isempty(readinput('basisvectors',optout)), optout = {optout{:}, 'basisvectors', basisvectors}; end
                if isempty(readinput('axesnames',optout))
                    if any(strcmpi(datalist.coordtype,{'QXYZ','QXQYEN'})), optout = {optout{:}, 'axesnames', {axname1,axname2,axname3}};
                    else optout = {optout{:}, 'axesnames', {axname1,axname2}}; end
                end
                
            catch
                fprintf('Error on transformation into HKL coordinates.\n');
                newlist = []; return;
            end
        else
            trygeneraltransform;
        end
                

    elseif strcmpi(coordtype,'QEPLANE')
        if strcmpi(datalist.coordtype,'A4ENERGY') || strcmpi(datalist.coordtype,'ANGLESENERGY')
            Es = datalist.(coordfield)(:,2);
            kis = sqrt(2*mass_n * Es / meVJ / hbar^2 + datalist.KF^2 * 1E20) * 1E-10;
            %Do the transformation to inplaneQ with beta=0; take then absolute values
            [newlist.(coordfield)(:,1), qp] = inplaneQ(datalist.(coordfield)(:,1), zeros(size(datalist.(coordfield)(:,1))), kis, datalist.KF, datalist.QVERT);
            newlist.(coordfield)(:,1) = sqrt (newlist.(coordfield)(:,1).^2 + qp.^2);  % ** ??
            newlist.(coordfield) = newlist.(coordfield)(:,1:2);
        elseif strcmpi(datalist.coordtype,'LINEARQ')
            % Q = Q0 + lambda*[-n(2),n(1)]
            Q0 = datalist.normal * datalist.c;
            Qs = [Q0(1) - datalist.normal(2) * datalist.(coordfield)(:,1), Q0(2) + datalist.normal(1) * datalist.(coordfield)(:,1)];
            Es = datalist.(coordfield)(:,2);
            [h,k,l] = calcHKL(Qs(:,1), Qs(:,2), 0, UBmatrix(datalist.sampleinfo.lattice,datalist.sampleinfo.ax,datalist.sampleinfo.bx));
            hkl = [h,k,l];  [~,whichx] = max(max(hkl,[],1)-min(hkl,[],1)); % which of h,k,l varies the most?
            newlist.(coordfield) = [hkl(:,whichx), Es];
            vn = {'H (r.l.u.)','K (r.l.u.)','L (r.l.u.)'};
            newlist.variables = {vn{whichx},'E (meV)'};
        else
            trygeneraltransform;
        end


    elseif strcmpi(coordtype,'QXQYEN')
        if strcmpi(datalist.coordtype,'ANGLESENERGY')
            Es = datalist.(coordfield)(:,2);
            kis = sqrt(2*mass_n * Es / meVJ / hbar^2 + datalist.KF^2 * 1E20) * 1E-10;
            %Do the transformation to inplaneQ 
            newlist.(coordfield)(:,3) = Es;
            [newlist.(coordfield)(:,1), newlist.(coordfield)(:,2)] = inplaneQ(datalist.(coordfield)(:,1), datalist.(coordfield)(:,3), kis, datalist.KF, datalist.QVERT);
        else
            trygeneraltransform;
        end

    elseif strcmpi(coordtype,'PROJECTION')
        % first argument in varargin are vectors (n rows are n vectors)
        % second argument in varargin is origin
        if nargin < 4, fprintf('Error: Missing parameters for Transformation PROJECTION.\n'); newlist = []; return; end
        vectors = varargin{1};
        origin = varargin{2}; origin = origin(:)';
        if ~(size(vectors,2)== numel(origin)) || ~(size(vectors,2) == size(datalist.(coordfield),2)), fprintf('Error: Dimension mismatch in transformation PROJECTION.\n'); newlist=[]; return; end
        newdims = size(vectors,1);
        stdratio = getoption('stdratio');
        ndims = size(datalist.(coordfield),2);
        if isfield(datalist,'coordtype') && isfield(stdratio,upper(datalist.coordtype))  % apply scale corrections, if existent
            ratio = stdratio.(upper(datalist.coordtype));
            for d = 1:ndims
                datalist.(coordfield)(:,d) = datalist.(coordfield)(:,d) * ratio(d);
                vectors(:,d) = vectors(:,d) * ratio(d);
                origin(d) = origin(d) * ratio(d);
            end
        end
        newlist.(coordfield) = zeros(size(datalist.(coordfield),1),newdims);
        for d = 1:newdims
            vi = vectors(d,:)';
            vi2 = vi' * vi;
            for np = 1:size(datalist.(coordfield),1)
                newlist.(coordfield)(np,d) = (datalist.(coordfield)(np,:) - origin) * vi / vi2;
            end
        end
        if nargin>4, newlist.coordtype = varargin{3}; end

    % 

    else
        trygeneraltransform;
    end
    
end %for (coordfield)


end %coordtransform

