function liststruct = makelist(scanlist, type, varargin)

% Construct a linear list of all data points with respective coordinates and counts rates
% 'type' determines the coordinates
%
% P. Steffens, 01/2013


% Use standard options (for normalization etc.):
[maxdeviate, normalizeto, normval, channels, vanacorr] = getoption('maxdeviate','normalizeto','normval','channels','vanacorr','check',varargin);
if nargin<2, type='Angles'; end

%%
% Determine if list contains the scans or only the filenames
if isstruct(scanlist) || iscell(scanlist)
    nscans = length(scanlist);
    nodata = [];
elseif ischar(scanlist)
    % Load scans
    [scanlist, nscans, nodata] = tasread(scanlist,'download');
else
    fprintf('Could not interprete input argument "scanlist".\nGive either the filenames as string, or variable containing scan structure.\n');
    nscans = 0;
end

if nscans==0
    fprintf('No data loaded! \n');
    liststruct=[];
    return;
end


%% Check if scan files contain data

nd=numel(nodata);
if nd>0
    fprintf('Warning: No data in file(s): ');
    for i=1:nd, fprintf([scanlist(nodata(i)).FILE ' ']); end
    fprintf('\n');
end

scanlist = scanlist(setdiff(1:nscans, nodata));
nscans = nscans-nd;

%% Loop over all scans

% Which case of Detector do we use ??
if isfield(scanlist(1),'MULTI'), flatconemode = true; else flatconemode = false; end

if isfield(scanlist(1).DATA,'ROI'), impsmode = true; za = getoption('impsanadist','check',varargin); else impsmode = false; end
if impsmode && any(channels>9), channels = 1:9; end  % if value in options.m is wrong (likely for FC)

if ~isempty(readinput('multiflexxanalyzer',varargin))
    flatconemode = true;
    anachannel = readinput('multiflexxanalyzer',varargin);
    multiflexxkfvalues = [1.0984, 1.2032, 1.2996, 1.3894, 1.4737]; % = sqrt(2*m*E/hbar^2) with E=[2.5, 3, 3.5, 4, 4.5] meV;
    kf = multiflexxkfvalues(anachannel);
else
    anachannel = 0;
end
    

singlemode = ~flatconemode && ~impsmode; % if neither FC nor IMPS, assume single detector


polarized = false;
paldeflist = {};

% Convert scans into lists with coordinates of desired type
for i=1:nscans
    
    if isstruct(scanlist), scan = scanlist(i); else scan = scanlist{i}; end
    
%     %artificially set a3p=0  **to be changed!!**
    scan.VARIA.A3P=0;
    
    
    if (~isfield(scan.VARIA,'GU') || ~isfield(scan.VARIA,'GL')) && ~any(strcmpi(varargin,'nogoniomode'))
        scan.VARIA.GL=0; scan.ZEROS.GL=0;
        scan.VARIA.GU=0; scan.ZEROS.GU=0;
        fprintf('No gonios found in file %s: gu, gl, zgu, zgl are set to zero (nogoniomode).\n',scan.FILE);
    % Set gonios to zero if in 2D mode
    elseif any(strcmpi(varargin,'2Dmode'))  
        % this does not set the zgl and zgu to 0; this allows for use of a3p (inclined) axis
        scan.ZEROS.GL = scan.ZEROS.GL - scan.VARIA.GL;   scan.VARIA.GL=0;
        scan.ZEROS.GU = scan.ZEROS.GU - scan.VARIA.GU;   scan.VARIA.GU=0;
        if isfield(scan.DATA,'GU'), scan.DATA.GU = zeros(scan.DATA.GU); end
        if isfield(scan.DATA,'GL'), scan.DATA.GL = zeros(scan.DATA.GL); end
    elseif any(strcmpi(varargin,'nogoniomode'))  
        % set really everything to zero 
        % (this is needed when a3p is used and no gonio installed; zgu and zgl then must be zero)
        scan.ZEROS.GL = 0;   scan.VARIA.GL=0;
        scan.ZEROS.GU = 0;   scan.VARIA.GU=0;
        if isfield(scan.DATA,'GU'), scan.DATA.GU = zeros(scan.DATA.GU); end
        if isfield(scan.DATA,'GL'), scan.DATA.GL = zeros(scan.DATA.GL); end
    elseif isfield (scan.PARAM,'GONIO') && scan.PARAM.GONIO == 0
        % if Gonios deactivated in MAD (in newer Mad versions only), like 2Dmode
        scan.ZEROS.GL = scan.ZEROS.GL - scan.VARIA.GL;   scan.VARIA.GL=0;
        scan.ZEROS.GU = scan.ZEROS.GU - scan.VARIA.GU;   scan.VARIA.GU=0;
        if isfield(scan.DATA,'GU') || isfield(scan.DATA,'GL')
            fprintf('GU/GL columns in %s are set to zero, because Gonios are deactivated.\n',scan.FILE);
            if isfield(scan.DATA,'GU'), scan.DATA.GU = scan.DATA.GU*0; end
            if isfield(scan.DATA,'GL'), scan.DATA.GL = scan.DATA.GL*0; end
        end
    else
    % Do nevertheless a plausibility test for Gonios 
        dv = datevec(scan.DATE);
        if dv(1)<2008 && all(isfield(scan.VARIA,{'GU','GL'})) && any(abs([scan.VARIA.GU,scan.VARIA.GL])>0.2)
            % old data, there was no 3D mode yet
            fprintf('Scan %s is old data with non-zero Gonios. You might want to use the switch "2DMode" or "nogoniomode" to ignore these values.\n',scan.FILE);
        elseif (isfield(scan.VARIA,'GL') && abs(scan.VARIA.GL)>30) || (isfield(scan.VARIA,'GU') && abs(scan.VARIA.GU)>30) || ...
               (isfield(scan.DATA,'GU') && max(abs(scan.DATA.GU))>30) || (isfield(scan.DATA,'GU') && max(abs(scan.DATA.GU))>30)         
            fprintf('In file %s the values for the Gonios are unrealistically large. Please check.\nYou can use the switch "2DMode" or "nogoniomode" to ignore them. Continue...\n',scan.FILE);
        end
    end
    
   
    liststruct.sampleinfo.lattice = getlattice(scan);
    [liststruct.sampleinfo.ax, liststruct.sampleinfo.bx] = getorientation(scan);
    
    % Determine the detector efficiency correction factors to be applied
    if flatconemode
        if vanacorr == 2, correctionfactors = getoption('det_eff', 'check', varargin); correctionfactors = correctionfactors / mean(correctionfactors); 
        else correctionfactors = detectorcorrection; % **
        end
    elseif impsmode 
        if vanacorr==2, correctionfactors = getoption('impsdet_eff', 'check', varargin); correctionfactors = correctionfactors /sum(correctionfactors)*9;
        elseif vanacorr == 1, correctionfactors = []; fprintf('Correction on Vana not yet iplemented for Imps\n');
        else correctionfactors = [];    
        end
    else correctionfactors = []; 
    end
    
    if impsmode, lthree = scan.PARAM.LTHREE * 100; end % in cm
   
    
    % Get the coordinates of the data points in the desired form, eventually
    % after suitable coordinate transform, and construct list accordingly
    
    if anachannel==0, kf = getvar(scan,'kf'); end
    if flatconemode && anachannel ==0, gfc = getvar(scan,'gfc'); else gfc=0; end
    
       
    if any(strcmpi(type, {'ANGLES', 'QPLANE', 'QXY'}))
        % Coordinate transformation into the angle-plane
        if flatconemode
            [xang,yang] = XYangles(getvar(scan,'a4'), gfc , channels, getvar(scan,'psi'));
        elseif impsmode
%             xang = getvar(scan,'a4');
%             yang = getvar(scan,'a3');
%             if ~hasfield(scan.DATA,'QH') % no Q-information from MAD
                [xang,yang] = XYanglesIMPS(getvar(scan,'a4'), getvar(scan,'a3'), getvar(scan,'a5'), getvar(scan,'roi'), za, lthree);
%             end
        end
        % Determine and set parameters
        liststruct.coordlist = [xang(:), yang(:)]; 
        liststruct.type = 'Const-Energy, Const-Qvert';
        liststruct.raw = 1;
        liststruct.coordtype = 'Angles';
        % [ki,kf,Qvert] = kikfqv(scan);
        ki = getvar(scan,'ki'); %kf = getvar(scan,'kf'); 
        if impsmode, Qv= 0; elseif flatconemode, Qv = Qvert(scan, 1); end
        if (max(ki)-min(ki) > maxdeviate.KI) || (max(kf)-min(kf) > maxdeviate.KF) || (max(Qv)-min(Qv) > maxdeviate.QVERT)
            fprintf('Warning: ki, kf, or Qvert not constant for all points in scan %s !! Going on with averages... (check if correct!!)\n The tolerances can be set in the options file.\n',scan.FILE);
        end
        liststruct.KI = mean(ki); liststruct.KF = mean (kf); liststruct.QVERT = mean(Qv);
        liststruct.constants = {'KI', 'KF', 'QVERT'}; 
        
        liststruct = coordtransform(liststruct,type); % if necessary, here the coord-transformation into QXY

    elseif any(strcmpi(type, {'ANGLESQZ','QXYZ'}))
        [xang,yang] = XYangles(getvar(scan,'a4'), gfc, channels, getvar(scan,'psi'));
        Qv = Qvert(scan, true);
        liststruct.coordlist = [xang(:), yang(:), repmat(Qv,numel(channels),1)];
        liststruct.type = 'Const-Energy';
        liststruct.raw = 1;
        liststruct.coordtype = 'AnglesQZ';
        ki = getvar(scan,'ki'); % kf = getvar(scan,'kf');
        if (max(ki)-min(ki) > maxdeviate.KI) || (max(kf)-min(kf) > maxdeviate.KF)
            fprintf('Warning: ki or kf not constant for all points in scan %s !! Going on with averages... (check if correct!!)\n The tolerances can be set in the options file.\n',scan.FILE);
        end
        liststruct.KI = mean(ki); liststruct.KF = mean (kf);
        liststruct.constants = {'KI', 'KF'}; 
        liststruct = coordtransform(liststruct,type); % if necessary, here the coord-transformation into QXYZ
        

    elseif impsmode && any(strcmpi(type,{'QxQyEn','linearQ'}))
        % kf = getvar(scan,'kf');
        ki = getvar(scan,'ki');
        [xang,yang] = XYanglesIMPS(getvar(scan,'a4'), getvar(scan,'a3'), getvar(scan,'a5'), getvar(scan,'roi'), za, lthree); % "a4,a3"
        % Q in Lab system
        qlx = - kf .* sind(xang);
        qly =   kf .* cosd(xang)  - ki ;
        qlz =   zeros(size(kf));
        % Transform to crystal system
        zerovals.gu = scan.ZEROS.GU; zerovals.gl = scan.ZEROS.GL; zerovals.a3 = scan.ZEROS.A3;zerovals.a3p = scan.ZEROS.A3P; 
        [qx,qy] = QSampleA3(getvar(scan,'a3'),getvar(scan,'gu'),getvar(scan,'gl'),getvar(scan,'a3p'),qlx,qly,qlz,zerovals);
        
        liststruct.type = 'IMPS general data: qx qy En';
        liststruct.raw = 1;
        liststruct.constants = {'QVERT'};
        liststruct.QVERT = 0;
        liststruct.coordlist = [qx,qy,getvar(scan,'EN')];
        liststruct.coordtype = 'QXQYEN';
        
        if strcmpi(type,'linearQ')  % try to go to 2 dimensions : Q, E
            % Check if all Q lie on a line
            % For this, fit a line:
            [normal,c] = fithyperplane(liststruct.coordlist(:,1:2));
            % Check distances
            if any(normal(1)*liststruct.coordlist(:,1) + normal(2)*liststruct.coordlist(:,2) > 0.02)  % There are some deviations
                % Check if a projection is explicitly given
                prvec = readinput('projectionnormal',varargin); prdis = readinput('projectiondistance',varargin);
                if isempty(prvec) || isempty(prdis)
                    fprintf(['The Q-vectors do not fall on one line as required for linear-Q mode !\n', ...
                             'Go on with a projection to a fitted line (which may be meaningless - please check!!)\n', ...
                             'You can explicitly provide a projection if desired.\n']);
                else
                    normal = prvec / sqrt(sum(prvec.^2));
                    c = prdis;
                end                         
            end 
            liststruct.constants = [liststruct.constants, 'normal', 'c'];
            liststruct.normal = normal;
            liststruct.c = c;
            lambdas = -normal(2) * liststruct.coordlist(:,1) + normal(1)*liststruct.coordlist(:,2); 
            % Q = Q0 + lambda * [-n(2),n(1)]
            liststruct.coordlist = [lambdas, liststruct.coordlist(:,3)];
            liststruct.coordtype = 'LINEARQ';
        end
                
        
   
    elseif any(strcmpi(type,{'ENERGY3D','AnglesEnergy','EnergyProj','A4Energy','QxQyEn'}))
        [xang,yang] = XYangles(getvar(scan,'a4'), gfc, channels, getvar(scan,'psi'));
        energy = repmat(getvar(scan,'en'),[1,numel(channels)]);
        if anachannel>0 % getvar(..en) does not work in this special case... **
            energy = repmat(getvar(scan,'ei') - 2.072128*kf.^2, [1,numel(channels)]); %ei-ef
        end
        liststruct.type = 'Energyscan, Const-Qvert';
        liststruct.raw = 1;
        if any(strcmpi(type,{'ENERGY3D','AnglesEnergy','QxQyEn'})), liststruct.coordlist = [xang(:), energy(:), yang(:)]; liststruct.coordtype = 'AnglesEnergy';
        else liststruct.coordlist = [xang(:), energy(:)]; liststruct.coordtype = 'A4Energy'; end

        Qv = Qvert(scan, 1); % kf = getvar(scan,'kf'); 
        if (max(kf)-min(kf) > maxdeviate.KF) || (max(Qv)-min(Qv) > maxdeviate.QVERT)
            fprintf('Warning: Qvert not constant for all points in scan %s !! Going on with averages... (check if correct!!)\n The tolerances can be set in the options file.\n',scan.FILE);
        end
        liststruct.KF = mean(kf); liststruct.QVERT = mean(Qv);
        liststruct.constants = {'KF', 'QVERT'};
        
        if strcmpi(type,'QXQYEN')
             liststruct = coordtransform(liststruct,'QXQYEN');
        end

    elseif strcmpi(type,'direct')
        variables = readinput('variables',varargin,'last');
        if isempty(variables), fprintf('Error: names of variables to be used not provided.\n'); liststruct=[]; return; end
        if ~iscell(variables), fprintf('Error: names of variables not correctly given.\n'); liststruct=[]; return; end
        liststruct.coordlist=[];
        liststruct.coordtype = 'General';
        liststruct.type = 'General';
        liststruct.raw = 1;
        try
            for v=1:length(variables)
                liststruct.coordlist = [liststruct.coordlist, getvar(scan,variables{v})];
            end
        catch
            fprintf('Error on reading variable %s.\n',variables{v}); liststruct=[]; return; 
        end
        liststruct.constants = {};    %** maybe add sth ??
        
    else
        fprintf('Error: Type not recognized. Could not create list.\n');
        liststruct = [];
        return;   
    end
    
    % Filename, Title, etc.
    liststruct.dataname = scan.FILE; 
    liststruct.expname  = scan.TITLE;
    
    % Temperature
    if any(fieldcheck(scan.DATA,{'TT','TRT'}))  % Check if temperature written in file
        liststruct.constants = [liststruct.constants, 'TEMP'];
        if fieldcheck(scan.DATA,'TT')  % if TT present, take TT. Else TRT
            liststruct.TEMP = mean(scan.DATA.TT);
        else
            liststruct.TEMP = mean(scan.DATA.TRT);
        end
    end


    % Check for Polarization
    if ~isfield(scan,'POLAN') && ~isempty(paldeflist) 
        assignpal = readinput('setpal',varargin);
        scan.DATA.PAL = ones(size(scan.DATA.PNT));
        if isempty(assignpal)
            fprintf('Error: File %d does not contain Polarization info. Use "setpal" option to combine with the others.\n',scannr);
            liststruct = []; return;
        end        
    elseif isfield(scan,'POLAN')
        if (i>1) && ~polarized, fprintf('Error: Trying to combine non-polarized with polarized data.\n'); liststruct = []; return; end
        polarized = true;
        % Analyze the information in POLAN and create (append) the list of
        % PAL-Definitions (paldeflist)
        [paldeflist, assignpal] = analyzepal(scan, paldeflist);
    end
    
    

    % Counting data of all points (linear list), with the desired normalization
    norm_measured = getvar(scan,normalizeto); 
    
    if any(norm_measured == 0), fprintf(['Warning: During normalization on ' normalizeto ' some zero values in ' normalizeto ' occured. Please check the settings for normalization!\n']); end
    if isempty(norm_measured), fprintf(['Error: Could not perform normalization: ' normalizeto ' column not found in scan file.\n']);  liststruct = []; return; end
    
    if flatconemode
        norm_measured = repmat(norm_measured,[1,numel(channels)]); 
        if ~isempty(correctionfactors), correctionfactors = repmat(correctionfactors(channels), size(norm_measured,1), 1); end
    elseif impsmode
        if ~isempty(correctionfactors), correctionfactors = correctionfactors(getvar(scan,'ROI'))'; end
    end

    monitorvalues = norm_measured;                            % The values of the desired 'monitor' (i.e. of M1, TIME, ..)
    if any(monitorvalues==0), fprintf([ 'Warning: Zeros detected in Monitor column (' normalizeto ')! This will cause trouble in the following...\n']); end
    if ~isempty(correctionfactors)
        eff_monvalues = monitorvalues .* correctionfactors;
        % Because of the different detector efficiencies, some channels are effectively measured "longer" than others, i.e. with different
        % effective monitor count rates (times, etc.). In this way the efficiency correction is considered.
    else
        eff_monvalues = monitorvalues;
        vanacorr = 0;
    end
    
    if polarized   %** only valid for FC...!
        pallist = getvar(scan,'PAL');
        pallist = repmat(pallist,[1,numel(channels)]);
        liststruct.pallist = assignpal(pallist(:));
    end
    liststruct.polarized = polarized;
        
    if flatconemode && anachannel==0 % Flatcone
        multidat  =                scan.MULTI(:,channels)     ./ eff_monvalues * normval ;       %Count rates
        multierr  =       sqrt(max(scan.MULTI(:,channels),1)) ./ eff_monvalues * normval ;       %Error bars, the max(..,1) avoids zero error bars... **
        
    elseif flatconemode % hypothetical "Multi-kf Flatcone", read each 5th line of MULTI block
        multidat  =                scan.MULTI(5*(0:numel(scan.DATA.PNT)-1)+anachannel, channels)     ./ eff_monvalues * normval ;       %Count rates
        multierr  =       sqrt(max(scan.MULTI(5*(0:numel(scan.DATA.PNT)-1)+anachannel, channels),1)) ./ eff_monvalues * normval ;       %Error bars, the max(..,1) avoids zero error bars... **
       
    elseif impsmode
        newimpsrois = readinput('reintegrateimps',varargin);
        multidetector = [];
        if ~isempty(newimpsrois) % Re-integrate the imps-multidetector
            fprintf('Re-Integrating the IMPS-Multidetector for file %s\n',scan.FILE);
            [scan,multidetector] = reintegrateimps(scan,newimpsrois);
        end    
        
        impsbgdef = readinput('impsbackground',varargin);
        if ~isempty(impsbgdef) && ~isempty(multidetector) % Subtract a background deduced from certain regions of the IMPS-multidetector
            fprintf('Subtracting background from IMPS-Multidetector for file %s\n',scan.FILE);
            [bgcnts,bgerr] = impsbackground(multidetector,newimpsrois,impsbgdef);
            multidat = zeros(size(scan.DATA.CNTS)); multierr = multidat;
            for nr=1:length(bgcnts)
                multidat(scan.DATA.ROI==nr) = (scan.DATA.CNTS (scan.DATA.ROI==nr) - bgcnts{nr}(:) ) ./ eff_monvalues(scan.DATA.ROI==nr) * normval ;
                multierr(scan.DATA.ROI==nr) = sqrt(max(scan.DATA.CNTS(scan.DATA.ROI==nr),1) + bgerr{nr}(:).^2) ./ eff_monvalues(scan.DATA.ROI==nr) * normval ;
            end
            liststruct.raw = 0;
        else
            multidat  =          scan.DATA.CNTS     ./ eff_monvalues * normval ;
            multierr  = sqrt(max(scan.DATA.CNTS,1)) ./ eff_monvalues * normval ;        
        end
        
    elseif singlemode
        multidat  =          scan.DATA.CNTS     ./ eff_monvalues * normval ;
        multierr  = sqrt(max(scan.DATA.CNTS,1)) ./ eff_monvalues * normval ; 
    
    
    end
        
    liststruct.valuelist = [multidat(:), multierr(:)];
    liststruct.monitorlist = [eff_monvalues(:), monitorvalues(:)];  % effective and true values of monitor

    liste{i} = liststruct; %#ok<AGROW>
end

if polarized, liste{1}.paldeflist = paldeflist; end




%% Make a single list (without averaging)

if ~isempty(i)
    % Combine
    liststruct = cmbavg(liste, 'noAvg');
else liststruct = []; return;
end

%%
switch vanacorr
    case 0, liststruct.properties.correction = 'No detector efficiency correction';
    case 1, liststruct.properties.correction = ['Detector efficiency correction on vanadium scan file', getoption('vanafile')];
    case 2, liststruct.properties.correction = 'Detector efficiency correction on given values';
end  

liststruct.properties.normalization = [normalizeto,'=',num2str(normval)];



