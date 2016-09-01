function colorlimits(figh)

% colorlimits(figh)
% Simple dialog window to set the limits of the colorscale
% figh is handle to the figure

% P. Steffens 08/2009 - 08/2014


% Find max and min of data in figure
try
    [figdat,windat] = getfiguredata(figh);
    allval = mergelist(figdat,'valuelist');
    dlims = [min(allval(:,1)), max(allval(:,1))];
    if strcmpi(windat.linlog,'LOG')
        m = min(allval(:,1));
        cd = log10(allval(:,1) + max(0,-m) + 1);
        dlims = [min(cd), max(cd)];
    end
    clims = caxis(windat.axeshandle);
catch   % if something is not present as it should, use [0,1]
    dlims = [0,1];
    clims = dlims;
end


smallstep =0;
bigstep  = 0;
calcstep;



% Some standard values for uicontrol positions
gr = 8;  ux = 50; uy= 100; uy2=200;

try
    scrnsz = get(0,'screensize');
    winpos = get(figh,'position');
    fp = [min(winpos(1)+winpos(3)+gr,scrnsz(3)-ux-80), min(winpos(2), scrnsz(4)-uy2-180), ux+80, uy2+160];
catch
    fp = get(0,'defaultfigureposition');
    fp = [fp(1) fp(2) ux+80, uy2+160];  
end

fig_props = {'name','Set color scale', 'color', get(0,'defaultUicontrolBackgroundColor'), 'resize', 'off', 'numbertitle', 'off', 'position', fp, ...
             'menubar', 'none', 'windowstyle', 'normal', 'closerequestfcn', 'delete(gcbf)' };
        
fig = figure(fig_props{:});




ok_btn      = uicontrol('style','pushbutton', 'string','Ok', 'position',[ux, gr, 50, 30], 'callback','delete(gcbf);');
reset_btn   = uicontrol('style','pushbutton', 'string','Reset', 'position',[ux, gr+40, 50, 30], 'UserData', [0,0], 'callback',@change);           

lab1  = uicontrol('style','text','string','Data range:','Fontsize',9,'position',[gr,fp(4)-20,fp(3)-gr,20]);               
lab1a = uicontrol('style','text','string',[num2str(dlims(1)) ' ... ' num2str(dlims(2))],'Fontsize',9,'position',[gr,fp(4)-40,fp(3)-gr,20]);               
lab2  = uicontrol('style','text','string','Color scale:','Fontsize',9,'position',[gr,fp(4)-70,fp(3)-gr,20]);
lab3  = uicontrol('style','text','string','Max','Fontsize',9,'position',[gr,uy2+65,ux-gr,20]);
lab4  = uicontrol('style','text','string','Min','Fontsize',9,'position',[gr,uy+65,ux-gr,20]);


loval  = uicontrol('style','edit','string',num2str(clims(1)),'Userdata',1,'position',[ux,uy+65,fp(3)-ux-gr, 20],'callback',@newval);               
upval  = uicontrol('style','edit','string',num2str(clims(2)),'Userdata',2,'position',[ux,uy2+65,fp(3)-ux-gr, 20],'callback',@newval);               

upm  = uicontrol('style','pushbutton', 'string','-',  'position',[ux,    uy2,    25, 25], 'UserData', [2,-1], 'callback',@change);
upp  = uicontrol('style','pushbutton', 'string','+',  'position',[ux,    uy2+30, 25, 25], 'UserData', [2,1], 'callback',@change);
upmm = uicontrol('style','pushbutton', 'string','--', 'position',[ux+30, uy2,    25, 25], 'UserData', [2,-2], 'callback',@change);
uppp = uicontrol('style','pushbutton', 'string','++', 'position',[ux+30, uy2+30, 25, 25], 'UserData', [2,2], 'callback',@change);

lom  = uicontrol('style','pushbutton', 'string','-',  'position',[ux,    uy,    25, 25], 'UserData', [1,-1], 'callback',@change);
lop  = uicontrol('style','pushbutton', 'string','+',  'position',[ux,    uy+30, 25, 25], 'UserData', [1,1], 'callback',@change);
lomm = uicontrol('style','pushbutton', 'string','--', 'position',[ux+30, uy,    25, 25], 'UserData', [1,-2], 'callback',@change);
lopp = uicontrol('style','pushbutton', 'string','++', 'position',[ux+30, uy+30, 25, 25], 'UserData', [1,2], 'callback',@change);



% movegui(fig)
% set(fig, 'visible','on'); drawnow;


%% Calcstep

    function calcstep
        diff=clims(2)-clims(1);
        smallstep = 10 ^ (floor(log10(diff))-1);
        bigstep = diff/4;
    end

%% Update
    function update
        try  caxis(windat.axeshandle,clims); catch end;
        set(upval,'string',num2str(clims(2)));
        set(loval,'string',num2str(clims(1)));
        calcstep;
    end
        
%% Newval
    function newval(hObject, evdata)
        n = get(hObject, 'Userdata');
        oldval = clims(n);
        try 
            clims(n) = str2double(get(hObject,'string'));
        catch end
        if isfinite(clims(n)) && (clims(2)>clims(1))
            update;
        else
            clims(n)=oldval;
            set(hObject,'string',num2str(oldval));
        end
    end

%% Change
    function change(hObject, evdata)
        ud = get(hObject,'Userdata');
        n = ud(1);
        c = ud(2);
        if abs(c)==1 % small change
            clims(n) = clims(n) + sign(c) * smallstep;
        elseif abs(c)==2 % big change
            prec = 10 ^ floor(log10(abs(clims(n)))-1);
            if clims(n)==0, prec = .1; end
            clims(n) = clims(n) + sign(c) * bigstep;
            clims(n) = round ( clims(n) / prec) * prec;
        else % Reset
            clims = dlims;
        end
        update;
    end

end
