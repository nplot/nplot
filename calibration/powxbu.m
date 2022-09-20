function varargout = powxbu(varargin)

% Type "powxbu" to open an interactive dialog window.
% *****************************
% Script version of "powxbu":
% (1)  powxbu create [powder] [ki] ([a4min]) ([a4max]) 
% (2)  powxbu set [reflections] ['ti'|'da4'|'np'] [value]
% (3)  powxbu writexbu ([filename]) ([reflections]) 
% (4)  powxbu transfer ([server]) ([directory]) 
% *****************************
% "powxbu create [powder] [ki] ([a4min]) ([a4max])" creates a list of reflections
%       [powder] : name of powder sample (see file powderdefinitions.m)
%       [ki]     : value of ki in A^-1
%       [a4min],[a4max] : limits of a4 range (optional)
%       Example: "powxbu Silicon 2.662"
% "powxbu set [reflections] ['ti'|'da4'|'np'] [value]" sets the scan parameters
%       [reflections] : list of reflections. "all" for all
%       Examples: "powxbu set all ti 1",  "powxbu set 2:5 da4 .15"
% "powxbu writexbu ([filename]) ([reflections])" creates the .xbu file
%       [filename]    : (optional) valid filename (including ".xbu")
%       [reflections] : (optional) which reflections to use
% "powxbu transfer" copy the file to the instrument computer (if possible)
%       [server],[directory] : (optional) server name and directory (use defaults if empty)

% P. Steffens, 07/2014




%% First, check for command line arguments (script version)
global passvar
if nargin > 1 && strcmpi(varargin{1},'create')
    try
        powdername = varargin{2};
        ki = str2double(varargin{3});
        if nargin>3, a4min = str2double(varargin{4}); else a4min = -120; end
        if nargin>4, a4max = str2double(varargin{4}); else a4max = 120; end
    catch
        fprintf('Error on reading input parameters of powxbu create. Type "help powxbu" for help.\n');
        return;
    end
    % determine powder
    materials = powderdefinitions; pind = 0;
    for ind=1:length(materials)
        if strcmpi(materials{ind}.name,powdername), pind = ind; end
    end
    if pind==0, fprintf('Error: Could not identify powder name. See file powderdefinitions.m.\n'); return; end
    [px.a4,px.d] = calcbraggangl(materials{pind}.lattice, materials{pind}.reflections, ki, a4min, a4max);
    px.name = materials{pind}.name;
    [px.ti,px.da4,px.np] = getoption('powparam.ti','powparam.da4','powparam.np');
    [~,px.a4,px.d,px.ki,px.ti,px.da4,px.np] = makesamesize(px.a4,px.d,ki,px.ti,px.da4,px.np);
    fprintf('Reflection list generated:\n');
    fprintf('Nb.     ki  a4(expected)    da4     np     ti\n');
    for n=1:numel(px.a4) fprintf('%3d %6.3f %13.2f %6.2f %6d %6d\n', n, px.ki(n), px.a4(n), px.da4(n), px.np(n), px.ti(n)); end
    fprintf('(Saved in workspace variable PXBU)\n');
    passvar = px; % pass variable pc as PXBU to caller's workspace
    evalin('caller','global passvar; PXBU=passvar; clear passvar;');
    clear passvar;
    return;
elseif nargin > 1 && strcmpi(varargin{1},'set')
    % Get PXBU
    if nargin<4, fprintf('Not enough parameters. See help powxbu.\n'); return; end
    try px = evalin('caller','PXBU');
    catch, fprintf('Error: Variable PXBU not found in workspace. (Run "powxbu create" first.)\n'); return; 
    end
    % Ref. list
    if strcmpi(varargin{2},'all'), reflections = 1:numel(px.a4); 
    else reflections = str2num(varargin{2});  %#ok<ST2NM>
        if isempty(reflections), fprintf('Cannot read reflection list. See help powxbu.\n'); return; end
    end
    % Value pair
    col = find(strcmpi(varargin{3},{'ti','np','da4'}));
    val = str2double(varargin{4});
    if isempty(col) || isnan(val), fprintf('Invalid value pair ti/da4/np. See help powxbu.\n'); return; end
    switch col
        case 1, px.ti(reflections) = val;
        case 2, px.np(reflections) = val;
        case 3, px.da4(reflections)= val;
    end
    fprintf('Updated reflection list:\n');
    fprintf('Nb.     ki  a4(expected)    da4     np     ti\n');
    for n=1:numel(px.a4) fprintf('%3d %6.3f %13.2f %6.2f %6d %6d\n', n, px.ki(n), px.a4(n), px.da4(n), px.np(n), px.ti(n)); end
    fprintf('(Saved in workspace variable PXBU)\n');
    passvar = px; % pass variable pc as PXBU to caller's workspace
    evalin('caller','global passvar; PXBU=passvar; clear passvar;');
    clear passvar;
    return;
elseif nargin >= 1 && strcmpi(varargin{1},'writexbu')
    % Get PXBU
    try px = evalin('caller','PXBU');
    catch, fprintf('Error: Variable PXBU not found in workspace. (Run "powxbu create" first.)\n'); return; 
    end
    %Filename
    if nargin > 1, filename = varargin{2};
    else filename = ['powder-', px.name, '-', lower(datestr(now,'ddmmmyy')),'.xbu']; end
    % Refl. list
    if nargin<3 || strcmpi(varargin{3},'all'), ind = 1:numel(px.a4); 
    else ind = str2num(varargin{3});  %#ok<ST2NM>
        if isempty(ind), fprintf('Cannot read reflection list. See help powxbu.\n'); return; end
    end
    fprintf('Save as %s\n',filename);
    powxbu_write(filename, px.ki(ind), px.d(ind), px.a4(ind), px.da4(ind), px.np(ind), px.ti(ind), [px.name,' powder'], 1);
    px.filename = filename; % save filename
    passvar = px; % pass variable pc as PXBU to caller's workspace
    evalin('caller','global passvar; PXBU=passvar; clear passvar;');
    clear passvar;
    return;  
elseif nargin >= 1 && strcmpi(varargin{1},'transfer')
    % Get PXBU (for filename)
    try px = evalin('caller','PXBU');
    catch, fprintf('Error: Variable PXBU not found in workspace. (Run "powxbu create" first.)\n'); return; 
    end
    if nargin>1, servername = varargin{2}; else servername = getoption('defaultserver'); end
    if nargin>2, serverpath = varargin{3}; else serverpath = ['/users/',servername,'/instr'];  end
    if serverpath(end) ~='/', serverpath = [serverpath,'/']; end
    [~, n, e]= fileparts(px.filename); destname=[n,e];
    fprintf('Attempt to copy to %s on %s...\n', [serverpath,destname], servername);
    putremote(px.filename, [serverpath,destname], servername, servername);
    return;
end







% POWXBU MATLAB code for powxbu.fig
%      POWXBU, by itself, creates a new POWXBU or raises the existing
%      singleton*.
%
%      H = POWXBU returns the handle to a new POWXBU or the handle to
%      the existing singleton*.
%
%      POWXBU('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in POWXBU.M with the given input arguments.
%
%      POWXBU('Property','Value',...) creates a new POWXBU or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before powxbu_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to powxbu_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help powxbu

% Last Modified by GUIDE v2.5 25-Jul-2014 11:03:09

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @powxbu_OpeningFcn, ...
                   'gui_OutputFcn',  @powxbu_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT




function putremote(source, target, user, server)
% Do remote copy
% This may depend on the OS. Under Windows, other command may be necessary if no scp installed.

% little piece of code like in "downloadfile.m"
[re,ou] = system(['scp "', source, '" ', user, '@', server, ':', target]);
if strfind(ou,'denied')
    fprintf('Access denied. Retry in new console window (enter password). Exit console manually.\nAvoid this complication by generating a key pair. Type help downloadfile for help.\n');
    [re,ou] = system(['scp ', source, ' ', user, '@', server, ':"', target, '" .']);  %retry same in new window (to allow for password typing). Note: in this case, [re,ou] always [0,''] ...
end
if re ~= 0
    if ~isempty(ou), fprintf('%s\n',ou); end
end


function handles = updatereflections(hObject,handles)
% update list to choose reflections
powderindex = get(handles.powderlist,'Value')-1;
if powderindex == 0
    handles.rlist = cell(0,4);
    handles.powdername = '';
else
    powders = powderdefinitions;
    powder = powders{powderindex};
    ki = str2double(get(handles.editKI,'String'));
    a4min = str2double(get(handles.edita4min,'String'));
    a4max = str2double(get(handles.edita4max,'String'));
    [~,~,~,handles.rlist] = calcbraggangl(powder.lattice, powder.reflections, ki, a4min, a4max);
    handles.powdername = powder.name;
end
set(handles.uitable1,'Data',handles.rlist(:,1:2));
guidata(hObject, handles);



% --- Executes just before powxbu is made visible.
function powxbu_OpeningFcn(hObject, eventdata, handles, varargin) %#ok<*INUSL,*DEFNU>
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to powxbu (see VARARGIN)

% Choose default command line output for powxbu
handles.output = hObject;

handles.xbufilename = [];
handles.scanlist = cell(0,7);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes powxbu wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = powxbu_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in powderlist.
function powderlist_Callback(hObject, eventdata, handles)

handles = updatereflections(hObject,handles);
button_reset_Callback(hObject,[],handles);




% --- Executes during object creation, after setting all properties.
function powderlist_CreateFcn(hObject, eventdata, handles) %#ok<*INUSD>
% hObject    handle to powderlist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
powders = powderdefinitions;
listdata = {'[Select a powder:]'};
for n=1:length(powders)
    listdata{n+1} = powders{n}.name;
end
set(hObject,'String',listdata);



function editKI_Callback(hObject, eventdata, handles) 

updatereflections(hObject,handles);


% --- Executes during object creation, after setting all properties.
function editKI_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editKI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_a4step_Callback(hObject, eventdata, handles)
% hObject    handle to edit_a4step (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_a4step as text
%        str2double(get(hObject,'String')) returns contents of edit_a4step as a double


% --- Executes during object creation, after setting all properties.
function edit_a4step_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_a4step (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_np_Callback(hObject, eventdata, handles)
% hObject    handle to edit_np (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_np as text
%        str2double(get(hObject,'String')) returns contents of edit_np as a double


% --- Executes during object creation, after setting all properties.
function edit_np_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_np (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_ti_Callback(hObject, eventdata, handles)
% hObject    handle to edit_ti (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_ti as text
%        str2double(get(hObject,'String')) returns contents of edit_ti as a double


% --- Executes during object creation, after setting all properties.
function edit_ti_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_ti (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in button_add.
function button_add_Callback(hObject, eventdata, handles)
% hObject    handle to button_add (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if ~isfield(handles,'listselection') || isempty(handles.listselection)
    ind = 1:size(get(handles.uitable1,'Data'),1); % if no selection, take all
else
    ind = unique(handles.listselection(:,1))';
end
for n = ind
    handles.scanlist{end+1,1} = str2double(get(handles.editKI,'String'));
    handles.scanlist{end,2}   = handles.rlist{n,4};
    handles.scanlist{end,3}   = str2double(get(handles.edit_a4step,'String'));
    handles.scanlist{end,4}   = str2double(get(handles.edit_np,'String'));
    handles.scanlist{end,5}   = str2double(get(handles.edit_ti,'String'));
    handles.scanlist{end,6}   = true;
    handles.scanlist{end,7}   = handles.rlist{n,3}; %d-value
end
set(handles.scantable,'Data',handles.scanlist(:,1:6));
guidata(hObject,handles);


% --- Executes on button press in button_reset.
function button_reset_Callback(hObject, eventdata, handles)
% hObject    handle to button_reset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.scantable,'Data',cell(0,6));
handles.scanlist = cell(0,7);
handles.xbufilename = [];
guidata(hObject,handles);



% --- Executes on button press in button_write.
function button_write_Callback(hObject, eventdata, handles)
% hObject    handle to button_write (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if size(get(handles.scantable,'Data'),1)==0 || ~any(cell2mat(handles.scanlist(:,6)))
    warndlg('No reflections selected.','Cannot create .xbu','modal');
    return;
end
if isfield(handles,'xbufilename') && ~isempty(handles.xbufilename), stdfilename = handles.xbufilename;
else stdfilename = ['powder-', handles.powdername, '-', lower(datestr(now,'ddmmmyy')),'.xbu']; end

[FileName,handles.pathname] = uiputfile('*','Save xbu file', fullfile(pwd,stdfilename));
if FileName==0, return; end
ki = cell2mat(handles.scanlist(:,1));
a4 = cell2mat(handles.scanlist(:,2));
da4= cell2mat(handles.scanlist(:,3));
np = cell2mat(handles.scanlist(:,4));
ti = cell2mat(handles.scanlist(:,5));
ind = cell2mat(handles.scanlist(:,6));
d  = cell2mat(handles.scanlist(:,7));
comment= [handles.powdername, ' powder'];
powxbu_write(fullfile(handles.pathname,FileName), ki(ind), d(ind), a4(ind), da4(ind), np(ind), ti(ind), comment, 1);
handles.xbufilename = FileName;
guidata(hObject,handles);


% --- Executes on button press in button_view.
function button_view_Callback(hObject, eventdata, handles)
% hObject    handle to button_view (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if ~isempty(handles.xbufilename)
    edit(fullfile(handles.pathname,handles.xbufilename));
end


function edita4min_Callback(hObject, eventdata, handles)
% hObject    handle to edita4min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edita4min as text
%        str2double(get(hObject,'String')) returns contents of edita4min as a double
updatereflections(hObject,handles);

% --- Executes during object creation, after setting all properties.
function edita4min_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edita4min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edita4max_Callback(hObject, eventdata, handles)
% hObject    handle to edita4max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edita4max as text
%        str2double(get(hObject,'String')) returns contents of edita4max as a double
updatereflections(hObject,handles);

% --- Executes during object creation, after setting all properties.
function edita4max_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edita4max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes when selected cell(s) is changed in uitable1.
function uitable1_CellSelectionCallback(hObject, eventdata, handles)
% hObject    handle to uitable1 (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) currently selecteds
% handles    structure with handles and user data (see GUIDATA)
handles.listselection = eventdata.Indices;
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function scantable_CreateFcn(hObject, eventdata, handles)
% hObject    handle to scantable (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
set(hObject,'Data',cell(0,7));


% --- Executes during object creation, after setting all properties.
function uitable1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to uitable1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
set(hObject,'Data',cell(0,4));


% --- Executes on button press in button_transfer.
function button_transfer_Callback(hObject, eventdata, handles)
% hObject    handle to button_transfer (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if isempty(handles.xbufilename), return; end
[username,servername] = getoption('defaultuser','defaultserver');
answ = inputdlg({'Server name','Destination path','Destination filename'}, ...
    'Try to tansfer .xbu file to instrument computer', 1, ...
    {servername, ['/users/',username,'/instr'], handles.xbufilename},'on');
if isempty(answ), return; end
servername = answ{1}; serverpath = answ{2}; serverfile = answ{3};
if ~strfind(servername,username)
    username = inputdlg(['User name on ',servername], 'Try to tansfer .xbu file to instrument computer',1); 
    if isempty(username), return; end
end
if serverpath(end) ~='/', serverpath = [serverpath,'/']; end
target = [serverpath,serverfile];
source = fullfile(handles.pathname,handles.xbufilename);
% remove Vol name and replace "\"->"/"
% source(source == '\') = '/';
if ~isempty(find(source==':',1,'last')), source = source((find(source==':',1,'last')+1):end); end
putremote(source, target, username, servername);


% --- Executes when entered data in editable cell(s) in scantable.
function scantable_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to scantable (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)
if ~isempty(eventdata.NewData)
    handles.scanlist{eventdata.Indices(1),eventdata.Indices(2)} = eventdata.NewData;
end
guidata(hObject,handles);
