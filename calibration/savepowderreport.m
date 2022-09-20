function filename = savepowderreport(fitinfo, filename)

if nargin<2
    filename = ['powder_',datestr(now,'yyyy-mm-dd_HH'),'h',datestr(now,'MM')];
end

if filename==1, file=1; %screen output
else   file = fopen(filename,'w');
end

fprintf(file,'Powder peak calibration with powcal on %s\n',datestr(now));
fprintf(file,' Nr     File      Ki  Expected  Measured  Calculated Meas-calc\n');
for i=1:length(fitinfo.powderdata.file)
    fprintf(file,'%3d %8s  %6.3f  %8.2f  %8.3f    %8.3f  %8.3f ', ...
        i, fitinfo.powderdata.file{i}, fitinfo.powderdata.ki(i), fitinfo.powderdata.scancenter(i), ...
        fitinfo.powderdata.fita4(i), fitinfo.fitteda4(i), fitinfo.powderdata.fita4(i)-fitinfo.fitteda4(i));
    if ~ismember(i,fitinfo.fileindex), fprintf(file,'(excl.)'); end
    fprintf(file,'\n');
end
fprintf(file,'\nZEROS:\n');
fprintf(file,'a1 : %7.2f ', fitinfo.oldzeros.a1); 
if fitinfo.fitvarindex(1)==0, fprintf(file,'(fixed during fit)\n');
else fprintf(file,' (old) -->  %7.2f (new)   [+/- %6.3f]\n', fitinfo.newzeros.a1, fitinfo.fittedoffset.da2/2); end
fprintf(file,'a2 : %7.2f ', fitinfo.oldzeros.a2); 
if fitinfo.fitvarindex(1)==0, fprintf(file,'(fixed during fit)\n');
else fprintf(file,' (old) -->  %7.2f (new)   [+/- %6.3f]\n', fitinfo.newzeros.a2, fitinfo.fittedoffset.da2); end
fprintf(file,'a4 : %7.2f ', fitinfo.oldzeros.a4); 
if fitinfo.fitvarindex(2)==0, fprintf(file,'(fixed during fit)\n');
else fprintf(file,' (old) -->  %7.2f (new)   [+/- %6.3f]\n', fitinfo.newzeros.a4, fitinfo.fittedoffset.da4); end

if filename~=1, fclose(file); end
