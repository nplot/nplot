function [det_eff, scan] = detectorcorrection(scan)

% Perform corrections on the scan data (Flatcone detector efficiency)
% Returns the correction factors to be applied. If input argument scan is
% given, corrected scan structure is returned.
%
% P. Steffens 10/2009


[det_eff, vanacorr, vanafile] = getoption('det_eff','vanacorr','vanafile');
det_eff = det_eff(:)';    
scan = [];

% Vanadium correction for FC
if vanacorr == 1
    vanascan = tasread(vanafile);
    if isempty(vanascan)
        fprintf('Error: File %s for detector efficiency correction could not be loaded. You may change the settings in the options file.\n',vanafile); 
        det_eff = []; 
        return;
    end
    det_eff = sum(vanascan.MULTI, 1);  % do simple summation...
end

if vanacorr == 0
    % no correction
    det_eff = [];
else
    det_eff = det_eff / mean(det_eff); % Normalize so that average efficiency is 1
end


if nargin>0
    scan.MULTI      = scan.MULTI ./ repmat(det_eff, size(scan.MULTI,1), 1);
    scan.DATA.CNTS  = scan.DATA.CNTS / det_eff(getvar(scan,'chan'));
end

