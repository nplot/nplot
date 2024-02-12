function data = nxsTAS(nxs)
% Translate nxs file structure in traditional TAS structure

% PS 06/23 - 02/24

    function setdata(target, source, warn)
        nxslevel = strsplit(source,'.');
        src = nxs;
        for li=1:length(nxslevel)
            if hasfield(src,nxslevel{li})
                src = src.(nxslevel{li});
            else
                if nargin>2 && warn
                    warning([source ' not found in data file.']);
                end
                return;
            end
        end
        targetlevel = strsplit(target,'.');
        if length(targetlevel)==1
            data.(target)=src;
        elseif length(targetlevel)==2
            data.(targetlevel{1}).(targetlevel{2}) = src;
        else
            error('Only two nested levels allowed.'); 
        end        
    end

data = [];
if isempty(nxs), return; end

if ~hasfield(nxs,'instrument_name')
    error('Error: instrument_name not found in data file.');
end

setdata('INSTR','instrument_name');

if ~hasfield(nxs,nxs.instrument_name) && hasfield(nxs,'instrument')
    nxs.instrument_name = 'instrument'; % in some nexus files this field is called "instrument" rather than "IN20" etc.
end
    

setdata('USER','user.name');
setdata('LOCAL','user.namelocalcontact');
setdata('FILE','file_name');
if hasfield(data,'FILE'), data.FILE = regexpmatch(data.FILE,'[^/\s]+$'); if iscell(data.FILE), data.FILE = data.FILE{1}; end, end
setdata('DATE','start_time');
setdata('TITLE','title');
setdata('TYPE','instrument_mode');
setdata('COMND',[nxs.instrument_name '.command_line.actual_command']);

setdata('POSQE.QH','sample.qh'); 
setdata('POSQE.QK','sample.qk'); 
setdata('POSQE.QL','sample.ql'); 
setdata('POSQE.EN','sample.en');
setdata('CURVE.MONO',[nxs.instrument_name '.Monochromator.automatic_curvature']); if hasfield(data,'CURVE') && hasfield(data.CURVE,'MONO'), if data.CURVE.MONO, data.CURVE.MONO='AUTO'; else, data.CURVE.MONO='MANU'; end, end
setdata('CURVE.ANA', [nxs.instrument_name '.Analyzer.automatic_curvature']);      if hasfield(data,'CURVE') && hasfield(data.CURVE,'ANA'),  if data.CURVE.ANA,  data.CURVE.ANA ='AUTO'; else, data.CURVE.ANA ='MANU'; end, end

if hasfield(data,'COMND')
    data.STEPS = scansteps(data.COMND);
end

setdata('PARAM.FX','sample.fx');
if hasfield(data,'PARAM') && hasfield(data.PARAM,'FX')
    if data.PARAM.FX == 2, data.PARAM.KFIX = nxs.(nxs.instrument_name).Analyser.kf; end
    if data.PARAM.FX == 1, data.PARAM.KFIX = nxs.(nxs.instrument_name).Monochromator.ki; end
end

setdata('PARAM.GONIO','sample.automatic_gonio]');
setdata('PARAM.DA',[nxs.instrument_name '.Analyser.d_spacing']);
setdata('PARAM.DM',[nxs.instrument_name '.Monochromator.d_spacing']);
setdata('PARAM.SM',[nxs.instrument_name '.Monochromator.sens']);
setdata('PARAM.SS','sample.sens');
setdata('PARAM.SA',[nxs.instrument_name '.Analyser.sens']);
setdata('PARAM.ALF1',[nxs.instrument_name '.Distance.alf1']); setdata('PARAM.ALF2',[nxs.instrument_name '.Distance.alf2']); 
setdata('PARAM.ALF3',[nxs.instrument_name '.Distance.alf3']); setdata('PARAM.ALF4',[nxs.instrument_name '.Distance.alf4']); 
setdata('PARAM.BET1',[nxs.instrument_name '.Distance.bet1']); setdata('PARAM.BET2',[nxs.instrument_name '.Distance.bet2']); 
setdata('PARAM.BET3',[nxs.instrument_name '.Distance.bet3']); setdata('PARAM.BET4',[nxs.instrument_name '.Distance.bet4']); 
setdata('PARAM.ETAM',[nxs.instrument_name '.Monochromator.mosaic']);
setdata('PARAM.ETAA',[nxs.instrument_name '.Analyser.mosaic']);
setdata('PARAM.AS','sample.unit_cell_a');       setdata('PARAM.BS','sample.unit_cell_b');       setdata('PARAM.CS','sample.unit_cell_c'); 
setdata('PARAM.AA','sample.unit_cell_alpha');   setdata('PARAM.BB','sample.unit_cell_beta');    setdata('PARAM.CC','sample.unit_cell_gamma'); setdata('PARAM.ETAS','sample.mosaic');
setdata('PARAM.AX','sample.ax'); setdata('PARAM.AY','sample.ay'); setdata('PARAM.AZ','sample.az');
setdata('PARAM.BX','sample.bx'); setdata('PARAM.BY','sample.by'); setdata('PARAM.BZ','sample.bz');
setdata('PARAM.REACTOR',[nxs.instrument_name '.source.power']);

setdata('VARIA.USESELECTOR',[nxs.instrument_name '.Monochromator.use_selector']);
setdata('VARIA.CHGUIDE',[nxs.instrument_name '.source.guide']);

for f = fieldnames(nxs.(nxs.instrument_name))'
    f = f{1}; %#ok<FXSET>
    if isfield(nxs.(nxs.instrument_name).(f),'offset_value')
        data.ZEROS.(f) = nxs.(nxs.instrument_name).(f).offset_value;
    end
    if isfield(nxs.(nxs.instrument_name).(f),'value')
        data.VARIA.(f) = nxs.(nxs.instrument_name).(f).value;
    end
end

data.FORMT= '';

setdata('POLAN',[nxs.instrument_name '.pal.pal_contents']);
if hasfield(data,'POLAN')
    if regexp(data.POLAN,'\w+'), data.POLAN = strsplit(data.POLAN,'|'); else data = rmfield(data,'POLAN'); end
end

if ~hasfield(nxs,'data_scan') || ~hasfield(nxs.data_scan,'scanned_variables') || ~hasfield(nxs.data_scan.scanned_variables,'variables_names') || ~hasfield(nxs.data_scan.scanned_variables,'data') || ~hasfield(nxs.data_scan.scanned_variables.variables_names,'label')
    warning('Warning: could not extract data from nxs file.');
else
    datanames = nxs.data_scan.scanned_variables.variables_names.label;
    ndat = 0;
    for i = 1:length(datanames)
        if strcmpi(datanames{i}, 'Monitor1'), dname='M1';
        elseif strcmpi(datanames{i}, 'Monitor2'), dname='M2';
        elseif strcmpi(datanames{i}, 'Time'), dname='TIME';
        elseif strcmpi(datanames{i}, 'Detector'), dname='CNTS';
        else, dname = datanames{i};
        end
        data.DATA.(dname) = nxs.data_scan.scanned_variables.data(:,i);
        ndat = numel(data.DATA.(dname));
    end
    % "Guess" PNT and PAL   !** to be made more stable when real PAL will exist in nexus !!
    datanames{end+1} = 'PNT';
    if hasfield(data,'POLAN') % polarized
        npal = numel(cell2mat(strfind(join(data.POLAN(:)),'CO '))); % number of counts in polarization Loop (deduced by counting "co"'s in POLAN)
        pnt =  repmat(1:ceil(ndat/npal),npal,1);
        data.DATA.PNT = pnt(1:ndat)';
        pal =  repmat(1:npal,ceil(ndat/npal),1)';
        data.DATA.PAL =pal(1:ndat)';
        datanames{end+1} = 'PAL';
    else
        data.DATA.PNT =  (1:ndat)' ;
    end
    data.DATA.columnames = datanames;
end


% setdata('MULTI','data_scan.detector_data.data'); ** how to detect if Flatcone ?


end
