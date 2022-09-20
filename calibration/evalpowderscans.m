function [results, oldzeros, filenames] = evalpowderscans(fileinput)

% automatically fit powder scans
% for use in TAS alignment
% tabledata : cell array {filename, a4theo, a4exp, diff, chi2, kfix, dm, a2, erra4}
% oldzeros  : za1,za2,za4 found in scans
% filenames : filenames ggf. with path

% P. Steffens 07/2014

if any(fileinput==',') && ~any(fileinput=='['), fileinput = ['[' fileinput ']']; end
files = multifilename(fileinput);
expdata = tasread(fileinput,'download'); % the choosen data is analyzed by 'tasread'
results=[];
filenames={};
filecount = 0;
a4errorstate = false;
za1=[]; za2=[]; za4=[];
for j=1:length(expdata) %loop over files
    Expdata = expdata(j);
    matchStr = regexpi(Expdata.COMND,'(?<=sc\s+a4\s+)\-?\d+.?\d*','match');
    if isempty(matchStr) % not found "sc a4" in COMND
        fprintf(['Warning: Cannot evaluate file ',Expdata.FILE,'. Data files must be of type "sc a4 ..."\n']); 
        continue;
    end
    a4center = str2double(matchStr{1});
    if isempty(za1), za1=Expdata.ZEROS.A1; za2=Expdata.ZEROS.A2; za4=Expdata.ZEROS.A4; end
    if any([za1,za2,za4] ~= [Expdata.ZEROS.A1, Expdata.ZEROS.A2, Expdata.ZEROS.A4])
        fprintf('Warning: In the chosen data files, not all the zeros of a1, a2, a4 are equal. Skip file %s.\n',Expdata.FILE);
        continue;
    end
    % Little consistency check for KI (not really necessary)
    if abs(2*asind(pi/Expdata.PARAM.DM/Expdata.PARAM.KFIX)*Expdata.PARAM.SM-Expdata.VARIA.A2)>.05
        fprintf('Warning: There is some inconsitency between the values of DM, KFIX and A2 in file %s. Please check.\n', Expdata.FILE);
    end
    % Consistency check for Bragg angle
    a4fromAS = 2 * asind(pi/Expdata.PARAM.AS/Expdata.PARAM.KFIX);
    dif = abs(abs(a4center) - a4fromAS);
    if dif > .1 && dif < 1 % for small deviations, use AS
        fprintf('Warning: Scan %s is not centered at the a4-value expected from the parameter PARAM.AS. Using the latter one. (Please check)\n', Expdata.FILE);
        a4theo = a4fromAS;
        a4errorstate = true;
    elseif dif>=1 % for large deviations, use scan center
        fprintf('Warning: Parameter PARAM.AS seems to be meaningless in scan %s (large deviation). Ignore and use scan center.\n', Expdata.FILE);
        a4theo = a4center;
        a4errorstate = true;
    else
        a4theo = a4center;
    end
        
    filecount = filecount+1;     
    filenames{filecount}  = files{j}; %#ok<*AGROW>
    
    fprintf('Analyzing file %s... ',Expdata.FILE);
    [~,S] = nplot(files{j}, 'fit', 'gauss', 'noplot', 'plotstyle or', 'nooutput');  % Do fit
    fprintf('A4 = %6.2f, Chi2 = %7.2f\n', S.optparam(2), S.chi2);

    % Create return structure
    results.file{filecount} = Expdata.FILE;
    results.scancenter(filecount) = a4theo;
    results.fita4(filecount) = S.optparam(2);
    results.erra4(filecount) = S.errors(2);
    results.chi2(filecount) = S.chi2;
    results.ki(filecount) = Expdata.PARAM.KFIX;
    results.dm(filecount) = Expdata.PARAM.DM;
    results.a2(filecount) = Expdata.VARIA.A2;
end

if a4errorstate
    fprintf(['Note: You have received a warning because for each scan a consistency check is performed if the lattice\n', ...
        'constant (AS) matches the center of the scan. Xbu-files created by powxbu automatically set this value.\n']);
end

oldzeros.a1 = za1;
oldzeros.a2 = za2;
oldzeros.a4 = za4;



