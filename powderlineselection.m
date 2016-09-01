function [dvals,proc] = powderlineselection

% dvals = powderlineselection
% Simple dialog window to select powder lines
% returns an array of d-values
% proc stands for the process (ki / kf)

% P. Steffens 05/2010



% Standard values for uicontrol positions
gr = 8;     bw = 80;  bh = 30;  cbw = 70;   cbh = 20;
fp = get(0,'defaultfigureposition');
fp = [fp(1) fp(2) 5*(gr+cbw) bh+9*(gr+cbh)];  % keep upper left corner fixed

fig_props = {'name','Select powder lines', 'color', get(0,'defaultUicontrolBackgroundColor'), 'resize', 'off', 'numbertitle', 'off', 'position', fp, ...
             'menubar', 'none', 'windowstyle', 'modal', 'visible', 'off', 'createfcn', '', 'closerequestfcn', 'delete(gcbf)' };
        
fig = figure(fig_props{:});

ok_btn      = uicontrol('style','pushbutton', 'string','Ok', 'position',[gr, gr, bw, bh], 'callback',@doOK);
cancel_btn  = uicontrol('style','pushbutton', 'string','Cancel', 'position',[2*gr+bw, gr, bw, bh], 'callback',@doCancel);           

lab1  = uicontrol('style','text','string','Aluminium','Fontsize',11,'position',[gr,bh+7*(gr+cbh),cbw,cbh]);               
c.al111  = uicontrol('style','checkbox','string','1 1 1','position',[gr,          bh+6*(gr+cbh),cbw,cbh],'Value',1);               
c.al200  = uicontrol('style','checkbox','string','2 0 0','position',[gr+1*(gr+bw),bh+6*(gr+cbh),cbw,cbh],'Value',1);               
c.al220  = uicontrol('style','checkbox','string','2 2 0','position',[gr+2*(gr+bw),bh+6*(gr+cbh),cbw,cbh],'Value',1);               
c.al113  = uicontrol('style','checkbox','string','1 1 3','position',[gr+3*(gr+bw),bh+6*(gr+cbh),cbw,cbh],'Value',1);               

lab2  = uicontrol('style','text','string','Copper','Fontsize',11,'position',[gr,bh+5*(gr+cbh),cbw,cbh]);               
c.cu111  = uicontrol('style','checkbox','string','1 1 1','position',[gr,          bh+4*(gr+cbh),cbw,cbh]);               
c.cu200  = uicontrol('style','checkbox','string','2 0 0','position',[gr+1*(gr+bw),bh+4*(gr+cbh),cbw,cbh]);               
c.cu220  = uicontrol('style','checkbox','string','2 2 0','position',[gr+2*(gr+bw),bh+4*(gr+cbh),cbw,cbh]);               
c.cu113  = uicontrol('style','checkbox','string','1 1 3','position',[gr+3*(gr+bw),bh+4*(gr+cbh),cbw,cbh]);               

lab3  = uicontrol('style','text','string','Other','Fontsize',11,'position',[gr,bh+3*(gr+cbh),cbw,cbh]);
lab4  = uicontrol('style','text','string','(Give d-spacings in Angstrom)','Fontsize',8,'position',[gr+cbw-10,bh+3*(gr+cbh)-2,3*cbw,cbh]);
c.u1  = uicontrol('style','edit','string','','position',[gr,          bh+2*(gr+cbh),cbw,cbh]);               
c.u2  = uicontrol('style','edit','string','','position',[gr+1*(gr+bw),bh+2*(gr+cbh),cbw,cbh]);               
c.u3  = uicontrol('style','edit','string','','position',[gr+2*(gr+bw),bh+2*(gr+cbh),cbw,cbh]);               
c.u4  = uicontrol('style','edit','string','','position',[gr+3*(gr+bw),bh+2*(gr+cbh),cbw,cbh]);               

lab5 = uicontrol('style','text','string','Process:','position',[4*gr+3*bw-50, bh+1*(gr+cbh)-5, 50, cbh]);
c.p1 = uicontrol('style','checkbox','string','kf''=ki','position',[4*gr+3*bw, bh+1*(gr+cbh), cbw, cbh],'Value',1);
c.p2 = uicontrol('style','checkbox','string','ki''=kf','position',[4*gr+3*bw, bh+gr, cbw, cbh],'Value',0);



movegui(fig)
set(fig, 'visible','on'); drawnow;

proc=[1,1];

guidata(gcf,c);

try
    uiwait(fig);
catch
    if ishandle(fig)
        delete(fig)
    end
end

if isappdata(0,'DialogAppData__')
    ad = getappdata(0,'DialogAppData__');
    dvals = ad.dvals;
    proc = ad.proc;
    rmappdata(0,'DialogAppData__')
else
    % figure was deleted
    dvals = 0;
end


%% OK callback
    function doOK(ok_btn, evd, listbox) %#ok
        c = guidata(gcf);
        dvals = [];
        if get(c.al111,'value'), dvals = [dvals, 4.049 / sqrt(3)]; end
        if get(c.al200,'value'), dvals = [dvals, 4.049 / sqrt(4)]; end
        if get(c.al220,'value'), dvals = [dvals, 4.049 / sqrt(8)]; end
        if get(c.al113,'value'), dvals = [dvals, 4.049 / sqrt(11)]; end
        if get(c.cu111,'value'), dvals = [dvals, 3.61 / sqrt(3)]; end
        if get(c.cu200,'value'), dvals = [dvals, 3.61 / sqrt(4)]; end
        if get(c.cu220,'value'), dvals = [dvals, 3.61 / sqrt(8)]; end
        if get(c.cu113,'value'), dvals = [dvals, 3.61 / sqrt(11)]; end
        f = str2double(get(c.u1,'string')); if ~isnan(f) && f>0, dvals = [dvals, f]; end
        f = str2double(get(c.u2,'string')); if ~isnan(f) && f>0, dvals = [dvals, f]; end
        f = str2double(get(c.u3,'string')); if ~isnan(f) && f>0, dvals = [dvals, f]; end
        f = str2double(get(c.u4,'string')); if ~isnan(f) && f>0, dvals = [dvals, f]; end
        ad.dvals = dvals;
        if get(c.p1,'value'), ad.proc(1)=1; else ad.proc(1)=0; end
        if get(c.p2,'value'), ad.proc(2)=1; else ad.proc(2)=0; end    
        setappdata(0,'DialogAppData__',ad);
        delete(gcbf);
    end

%% Cancel callback
    function doCancel(cancel_btn, evd, listbox) %#ok
        ad.dvals = 0;
        ad.proc = [];
        setappdata(0,'ListDialogAppData__',ad)
        delete(gcbf);
    end


end
