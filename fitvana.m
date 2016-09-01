function [inten, fwhm, e0, di, df, de] = fitvana(scan)

% Fit Gaussian to all channels of Vanadium scan
%
% P. Steffens, 04/2008


%%
% If argument is filename, load scan first
if ~isstruct(scan)
    [scan, nscans] = tasread(scan);
end

energy = getvar(scan,'EN');

f = figure ;


for i = 1:31
    [y,pa,dpa] = funcfit(fsum(@const,@gaussA),energy,scan.MULTI(:,i),sqrt(scan.MULTI(:,i)),[0,18.6,max(scan.MULTI(:,i)),1.5],[1,1,1,1]);
    inten(i)= pa(3);
    fwhm(i) = pa(4);
    e0(i)   = pa(2);
    di(i)   = dpa(3);
    df(i)   = dpa(4);
    de(i)   = dpa(2);    
    errorbar(energy, scan.MULTI(:,i), sqrt(scan.MULTI(:,i)),'ob');
    hold;
    plot(energy, y, '-r');
    hold;
    title(['Channel ' num2str(i)]);
    pause;
end

close(f);
