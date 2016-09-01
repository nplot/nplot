function cscale(a1,a2)


try
    if findstr(a1,'lin')    
        [d,ps] = getfiguredata;
        ps.linlog = 'lin';
        doplot(ps,'keepview');
        return
    elseif strcmpi(a1,'log')
        [d,ps] = getfiguredata;
        ps.linlog = 'log';
        doplot(ps,'keepview');
        return
    end
catch
end

caxis([str2double(a1),str2double(a2)]);
