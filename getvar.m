function erg = getvar(scan,varname,varargin)

% Extract a Scan variable
% If it is listed in the file, just returns the listing (trivial).
% Otherwise, calculate it from the available information.
% Works currently for ki,kf,ei,ef,a3,a3p,a4,gfc,psi,gu,gl,qh,qk,ql,en.
%
% (Attention: the solution is sometimes not unique if the scan mode (e.g.
% constant a3p etc.) is not known. As this can at present not be inferred
% from file, the behavior is controlled by the stdscanmode parameter in the
% options-file!! This may concern a3,a3p,gu,gl.)
%
% Suitable not only for FC data (if no gfc/chan is found, set gfc=0 and
% chan=16 to make appropriate for standard TAS). 
% 
% If varargin contains  'includeparam' or 'readheader', then look also in PARAM, ZEROS and VARIA section
%
% ** what happens in a pure dEN - Scan?? **
% ** Zero-handling to be verified (should be OK) **
% ** make IMPS-compatible!!!! **
%
% P. Steffens, 11/2011





varname = upper(varname);

%If the field already exists:  nothing to do, be happy and exit 
if fieldcheck(scan.DATA,varname)
    erg = scan.DATA.(varname);
    return;
end

impsmode = isfield(scan,'impsmode') && scan.impsmode; % is this an IMPS data set?
impsrefchannel = 5;  %% This may be generalized later **
impsanadist = 4; %cm

UB = UBmatrix(getlattice(scan),[scan.PARAM.AX,scan.PARAM.AY,scan.PARAM.AZ], [scan.PARAM.BX,scan.PARAM.BY,scan.PARAM.BZ]); % ** Do not always need
erg = [];

[hbar, mass_n, meVJ,config.za3_0, config.zgu_0, config.zgl_0, config.GUsign, config.GLsign] = ...
    getoption('hbar', 'mass_n', 'meVJ','za3_0', 'zgu_0', 'zgl_0', 'GUsign', 'GLsign');

npoints = 0;
if ~isempty(scan.DATA) && isfield(scan.DATA,'PNT'), npoints = length(scan.DATA.PNT); end


if fieldcheck(scan.ZEROS,'A3')
    zerovals.a3  = scan.ZEROS.A3;
else
    zerovals.a3=0;
    scanmode = 'A3=0';
end
    
zerovals.gu  = scan.ZEROS.GU;
zerovals.gl  = scan.ZEROS.GL;

if fieldcheck(scan.ZEROS,'A3P') 
    % a3p is listed among variables (though not sure if or how it is used)
    zerovals.a3p = scan.ZEROS.A3P;
    scanmode = getoption('stdscanmode');  %This can at present not be inferred from file...!
else
    % in this case we can be sure that no a3p is installed!
    zerovals.a3p = 0;
    scanmode = 'A3P=0';
    if strcmp(varname,'A3P'), erg = zeros(npoints,1); end
end
   


%***************************** 
%helper function
    function CalcAnglesFromQs
        for n=1:npoints
            %Q in sample cartesian system
            [vx(n), vy(n), vz(n)] = calcQS(scan.DATA.QH(n), scan.DATA.QK(n), scan.DATA.QL(n), UB); %#ok<AGROW>
        end
        [a3,gu,gl,a3p,a4,gfc] = spectroangles(vx, vy, vz, getvar(scan,'KI'), getvar(scan,'KF'), getvar(scan,'CHAN'), scan.PARAM.SS, [], scanmode, zerovals);
    end  
%*****************************


%GFC-----------------------------------------------------------------------
if strcmp(varname,'GFC')
    if ~fieldcheck(scan.STEPS, 'QH') %gfc constant (gfc changes only in q-scan or gfc-scan)
        if fieldcheck(scan.VARIA,'GFC')
            erg = scan.VARIA.GFC * ones(npoints,1);
        else
            erg = zeros(npoints,1);
        end
    else
        CalcAnglesFromQs;
        erg = gfc(:);
    end
    return;
end

%A3P-----------------------------------------------------------------------
if strcmp(varname,'A3P')
    % changes in a3p, q- and psi
    if ~any(fieldcheck(scan.STEPS, {'QH','PSI'})) || strcmpi(scanmode(1:4),'A3P=') %a3p constant 
        if fieldcheck(scan.VARIA,'A3P') && isfinite(scan.VARIA.A3P)
            erg = scan.VARIA.A3P * ones(npoints,1);
        else
            erg = zeros(npoints,1);
        end
    else
        CalcAnglesFromQs;
        erg = a3p(:);
    end
    return;
end

%PSI-----------------------------------------------------------------------
if strcmp(varname,'PSI')
    %psi changes in psi-Scans or as fct. of a3,gu,gl,a4,gfc
    erg = fcpsi(getvar(scan,'A3'), getvar(scan,'GU'), getvar(scan,'GL'), getvar(scan,'A3P'), getvar(scan,'A4'), getvar(scan,'GFC'), zerovals);
    return;
end

%A4------------------------------------------------------------------------
if strcmp(varname,'A4')
    %A4 changes only in a4- or Q-scans
    if ~fieldcheck(scan.STEPS, 'QH')
        erg = scan.VARIA.A4 * ones(npoints,1);
        return;
    end
    CalcAnglesFromQs;
    erg = a4(:);
    return;
end

%A3------------------------------------------------------------------------
if strcmp(varname,'A3')
    %a3 changes only in a3-, psi- or Q-Scans
    if ~any(fieldcheck(scan.STEPS, {'QH','PSI'}))
        erg = scan.VARIA.A3 * ones(npoints,1);
        return;
    end
    if fieldcheck(scan.STEPS, 'PSI')
        [a3,gu,gl,a3p]= a3gugl(getvar(scan,'A4'), getvar(scan,'GFC'), scan.DATA.PSI, scanmode, zerovals, config);
        erg = a3(:);
    else % (Q-Scan)
        CalcAnglesFromQs;
        erg = a3(:);
    end
    return;
end    

   
%GU------------------------------------------------------------------------
if strcmp(varname,'GU')
    %GU changes only in gu-, psi-, or Q-scans
    if ~any(fieldcheck(scan.STEPS, {'QH','PSI'}))
        erg = scan.VARIA.GU * ones(npoints,1);
        return;
    end
    if fieldcheck(scan.STEPS, 'PSI')
        [a3,gu,gl,a3p]= a3gugl(getvar(scan,'A4'), getvar(scan,'GFC'), scan.DATA.PSI, scanmode, zerovals, config);
        erg = gu(:);
    else % (Q-Scan)
        CalcAnglesFromQs;
        erg = gu(:);
    end
    return;
end

%GL------------------------------------------------------------------------
if strcmp(varname,'GL')
    %GL changes only in gl-, psi-, or Q-scans
    if ~any(fieldcheck(scan.STEPS, {'QH','PSI'}))
        erg = scan.VARIA.GL * ones(npoints,1);
        return;
    end
    if fieldcheck(scan.STEPS, 'PSI') % (psi-Scan)
        [a3,gu,gl,a3p]= a3gugl(getvar(scan,'A4'), getvar(scan,'GFC'), scan.DATA.PSI, scanmode, zerovals, config);
        erg = gl(:);
    else        % (Q-Scan)
        CalcAnglesFromQs;
        erg = gl(:);
    end
    return;
end



%Q-------------------------------------------------------------------------
if any(strcmp(varname, {'QH','QK','QL'}))
    if impsmode
        [qlx,qly,qlz] = QLabIMPS(getvar(scan,'A4'),getvar(scan,'A5'),getvar(scan,'KI'), getvar(scan,'KF'),scan.DATA.ROI,scan.PARAM.LTHREE);
    else
        %need a4,gfc,ch,ki,kf to calculate QLab
        [qlx, qly, qlz] = QLab(getvar(scan,'A4'), getvar(scan,'GFC'), getvar(scan,'CHAN'), getvar(scan,'KI'), getvar(scan,'KF'));
    end
    %now, need sample orientation: either (psi) or (a3,gu,gl,a3p)
    [qx,qy,qz] = QSampleA3(getvar(scan,'A3'), getvar(scan,'GU'),  getvar(scan,'GL'), getvar(scan,'A3P'), qlx, qly, qlz);
    %Convert from sample holder system to HKL by inv. UB-Matrix
    for j=1:npoints
        [QH(j), QK(j), QL(j)] = calcHKL(qx(j),qy(j),qz(j),UB); %#ok<AGROW,NASGU>
    end
    erg = eval(varname);
    erg = erg(:);
    return;    
end

%KF------------------------------------------------------------------------
if strcmp(varname,'KF')
    if ~impsmode && scan.PARAM.FX == 2
        erg = scan.PARAM.KFIX * ones(npoints,1); return;
    end
    
    if ~impsmode && ~any(fieldcheck(scan.STEPS, {'EN','EF'}))
        scan.DATA.EN = scan.POSQE.EN * ones(npoints,1);
    end
    
    if fieldcheck(scan.DATA,'EF')
        erg = sqrt(scan.DATA.EF / meVJ * 2 * mass_n) / hbar * 1E-10; return
    end
    
    if fieldcheck(scan.DATA,'EN')
        scan.DATA.EF = hbar^2 * getvar(scan,'KI').^2 * 1E20 / 2 / mass_n * meVJ  -  scan.DATA.EN;
        erg = sqrt(scan.DATA.EF / meVJ * 2 * mass_n) / hbar * 1E-10; return
    end
        
    if impsmode % (a priori, this must be the case at this point)
        % There is neither KF, EF, EN 
        % kf's change only if a5 and/or cry's are scanned
        beta = getvar(scan,'a5'); 
        nAna = getvar(scan,'ROI'); % ** To be seen... !!
        delta    = atand( - (nAna-5) * impsanadist .* cosd(beta) ./ ( (nAna-5)* impsanadist .* sind(beta) + scan.PARAM.LTHREE*100 ) );
        erg = nan(size(scan.DATA.PNT,1),1);
        for i=1:9
            psi = getvar(scan,['CRY',num2str(i)]);
            theta = psi + beta - delta;
            erg(nAna==i) = pi/scan.PARAM.DA ./ sind(abs(theta(nAna==i)));
        end
        return
    end
       
%     if scan.PARAM.FX == 2
%         erg = scan.PARAM.KFIX * ones(npoints,1); return;
%     end
%     if ~any(fieldcheck(scan.STEPS, {'EN','EF'}))
%         scan.DATA.EN = scan.POSQE.EN * ones(npoints,1);
%     end
%     if fieldcheck(scan.DATA,'EN')
%         scan.DATA.EF = hbar^2 * scan.PARAM.KFIX^2 * 1E20 / 2 / mass_n * meVJ  -  scan.DATA.EN;
%     end
%     if fieldcheck(scan.DATA,'EF')
%         erg = sqrt(scan.DATA.EF / meVJ * 2 * mass_n) / hbar * 1E-10;
%     end
%     return;
end
    
%KI------------------------------------------------------------------------
if strcmp(varname,'KI')
    if scan.PARAM.FX == 1
        erg = scan.PARAM.KFIX * ones(npoints,1); return;
    end
    if ~any(fieldcheck(scan.STEPS, {'EN','EI'}))
        % constant, POSQE.EN + kf^2....
        scan.DATA.EI = hbar^2 * scan.PARAM.KFIX^2 * 1E20 / 2 / mass_n * meVJ + scan.POSQE.EN * ones(npoints,1);
    end
    if fieldcheck(scan.DATA,'EN')
        scan.DATA.EI = hbar^2 * scan.PARAM.KFIX^2 * 1E20 / 2 / mass_n * meVJ  +  scan.DATA.EN;
    end
    if fieldcheck(scan.DATA,'EI')
        erg = sqrt(scan.DATA.EI / meVJ * 2 * mass_n) / hbar * 1E-10;
    end     
    return;
end

%EF------------------------------------------------------------------------
if any(strcmp(varname, 'EF'))
    erg = getvar(scan,'KF').^2 * hbar^2 * 1E20 / 2 / mass_n *meVJ;
    return;
end

%EI------------------------------------------------------------------------
if any(strcmp(varname, 'EI'))
    erg = getvar(scan,'KI').^2 * hbar^2 * 1E20 / 2 / mass_n *meVJ;
    return;
end

%EN------------------------------------------------------------------------
if any(strcmp(varname, 'EN'))
    erg = (getvar(scan,'KI').^2 - getvar(scan,'KF').^2) * hbar^2 * 1E20 / 2 / mass_n *meVJ;
    return;
end

%QVERT-----(treat as var...)-----------------------------------------------
if any(strcmp(varname, 'QVERT'))
    erg = Qvert(scan, true);
    return;
end

%A5------------------------------------------------------------------------
if any(strcmp(varname, 'A5'))
    if fieldcheck(scan.VARIA,'A5')
        erg = ones(npoints,1) * scan.VARIA.A5;
    return;
    end
end

%A6------------------------------------------------------------------------
if any(strcmp(varname, 'A6'))
    if fieldcheck(scan.VARIA,'A6')
        erg = ones(npoints,1) * scan.VARIA.A6;
    return;
    end
end

%CRY-----------------------------------------------------------------------
if any(strcmp(varname(1:min(end,3)), 'CRY'))
    if any(fieldcheck(scan.STEPS, {'EN','EF','KF','QH'}))
        % Here, CryX may move 
        crynum = str2double(varname(4));
        theta = scan.PARAM.SA * asind(pi/scan.PARAM.DA./getvar(scan,'KF'));
        for i=1:max(scan.DATA.PNT)
            theta(scan.DATA.PNT==i) = theta(scan.DATA.PNT==i & scan.DATA.ROI==crynum);
        end
        beta  = getvar(scan,'a5'); 
        delta = atand( (crynum-5) * impsanadist * sind(beta) ./ ( (crynum-5)* impsanadist * cosd(beta) - scan.PARAM.LTHREE*100 ) );
        erg = delta + theta -beta;
        return
    end
    if fieldcheck(scan.VARIA,varname)
        erg = ones(npoints,1) * scan.VARIA.(varname);
        return;
    end
    if strcmpi(varname,'CRY') && fieldcheck(scan.VARIA,'CRY1') % give all nine CRY's
        for i=1:9
            erg =  [erg, ones(npoints,1) * scan.VARIA.(['CRY',num2str(i)]) ]; 
        end
        return
    end
end



%CHAN---------------------------(though no 'variable' in strict sense)-----
if strcmp(varname,'CHAN')
    if fieldcheck(scan.PARAM,'CHAN')
        erg = scan.PARAM.CHAN;
    else
        erg = 16;
    end
    return;
end



% --- Look through Params and Zeros -----
    
if nargin > 2 && (any(strcmpi(varargin,'includeparam')) || any(strcmpi(varargin,'readheader')))
    if fieldcheck(scan.VARIA,varname)
        erg = scan.VARIA.(varname);
        return;
    end    
    if fieldcheck(scan.PARAM,varname)
        erg = scan.PARAM.(varname);
        return;
    end
    if varname(1)=='Z' && fieldcheck(scan.ZEROS,varname(2:end))
        erg = scan.ZEROS.(varname(2:end));
        return;
    end
end


disp(['Attention! Could not extract variable ' varname]);

end

      