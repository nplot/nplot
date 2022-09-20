function varargout = powcal(varargin)

% Type "powcal" to open an interactive dialog window.
% *****************************
% Script version of "powcal":
% (1)  powcal readfiles [filename] 
% (2)  powcal fitzeros
% (3)  powcal savereport 
% *****************************
% "powcal readfiles [filename]" loads the data files containing a4 scans.
%    [filename] can be given in format like 0123[45:67], 01234[0,2,4:7], etc.
%    The results are stored in the workspace variable "PCAL" and can be edited if desired
% "powcal fitzeros" fits the a2 and a4 offsets using the information stored in the "PCAL" variable
%    To use another workspace variable, type "powcal fitzeros [varname]".
%    To exclude scans, use PCAL.fileindex. To fix za2 or za4, set PCAL.fitvarindex to 0.
%    The results are again stored in the workspace variable "PCAL" 
% "powcal savereport [filename]" writes the report into the file [filename]
%    If [filename] is left empty, the default filename is ("powder_[date]_[time]" in the current folder)
%    To use the content of another workspace variable than "PCAL", type "powcal savereport [filename] [varname]"


%% First, check for command line arguments (script version)
global passvar
if nargin > 1 && strcmpi(varargin{1},'readfiles') % Evaluate a4 scans
    [pc.powderdata, pc.oldzeros] = evalpowderscans(varargin{2});
    if isempty(pc.powderdata), fprintf('No files read.\n'); return; end
    pc.fileindex = 1:numel(pc.powderdata.fita4); % By default, include all data
    pc.fitvarindex = [1,1]; % By default, fit a2 and a4
    passvar = pc;
    % pass variable pc as PCAL to caller's workspace
    evalin('caller','global passvar; PCAL=passvar; clear passvar;');
    clear passvar;
    fprintf('Results saved in workspace variable "PCAL".\n');
    fprintf('Change if necessary, and run "powcal fitzeros" for fitting.\n');
    return;
elseif nargin >= 1 && strcmpi(varargin{1},'fitzeros')   % Do fitting
    if nargin >=2, varname = varargin{2}; else varname = 'PCAL'; end
    try pc = evalin('caller',varname);
    catch, fprintf('Error: Variable name %s not found.\n',varname); return; 
    end
    [pc.newzeros,pc.fittedoffset,pc.fitteda4] = fita2a4zeros(pc.powderdata, pc.oldzeros, pc.fitvarindex, pc.fileindex);
    % pass variable pc as PCAL to caller's workspace (as above)
    passvar = pc;
    evalin('caller','global passvar; PCAL=passvar; clear passvar;');
    clear passvar;
    savepowderreport(pc,1); % Show report on screen
    fprintf('Results saved in workspace variable "PCAL".\n');
    fprintf('Run "powcal savereport [filename] ([pcal])" to save in a file.\n');
    return;
elseif nargin >= 1 && strcmpi(varargin{1},'savereport')   % Write Report file
    % obtain filename
    if nargin<2 || isempty(varargin{2})
        reportfilename = ['powder_',datestr(now,'yyyy-mm-dd_HH'),'h',datestr(now,'MM')];
    else reportfilename = varargin{2};
    end
    % obtain variable content from parent workspace
    if nargin >2, varname = varargin{3}; else varname = 'PCAL'; end
    try pc = evalin('caller',varname);
    catch, fprintf('Error: Variable name %s not found.\n',varname); return; 
    end
    % Save to disk
    savepowderreport(pc,reportfilename);
    fprintf('Report saved to %s.\n',reportfilename);
    return;
end
clear passvar
    


%% powcal MATLAB code for powcal.fig
%      powcal, by itself, creates a new powcal or raises the existing
%      singleton*.
%
%      H = powcal returns the handle to a new powcal or the handle to
%      the existing singleton*.
%
%      powcal('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in powcal.M with the given input arguments.
%
%      powcal('Property','Value',...) creates a new powcal or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before powcal_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to powcal_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help powcal

% Last Modified by GUIDE v2.5 22-Jul-2014 18:52:20

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @powcal_OpeningFcn, ...
                   'gui_OutputFcn',  @powcal_OutputFcn, ...
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


% --- Executes just before powcal is made visible.
function powcal_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to powcal (see VARARGIN)

handles.filepath = '';

% Choose default command line output for powcal
handles.output = hObject;

set(handles.fittable,'Data',[]);

% Update handles structure
guidata(hObject, handles);

% --- Outputs from this function are returned to the command line.
function varargout = powcal_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



% --- Evaluate Data Files.
function run_Callback(hObject, eventdata, handles) %#ok<*INUSL,DEFNU>

fileinput = get(handles.edit1, 'String');
[handles.powderdata, handles.oldzeros, handles.filenames] = evalpowderscans(fileinput);
guidata(hObject,handles);
updateTable(hObject,handles);
updateFitTable(hObject,handles);


function updateTable(hObject,handles)
% Fill table view
if ~isfield(handles,'powderdata') || isempty(handles.powderdata), tabledata = []; 
else
    for filecount=1:length(handles.powderdata.file)
        tabledata{filecount,1} = handles.powderdata.file{filecount}; %#ok<*AGROW>
        tabledata{filecount,2} = num2str(handles.powderdata.scancenter(filecount),'%7.2f');
        tabledata{filecount,3} = num2str(handles.powderdata.fita4(filecount),'%8.3f');
        tabledata{filecount,4} = num2str(handles.powderdata.fita4(filecount) - handles.powderdata.scancenter(filecount),'%8.3f');;
        if handles.powderdata.chi2(filecount)==0
            tabledata{filecount,5} = '--';
        else
            tabledata{filecount,5} = num2str(handles.powderdata.chi2(filecount),'%8.3f');
        end
    end
end
set(handles.uitable1,'Data',tabledata);
set(handles.oza1,'String',num2str(handles.oldzeros.a1,'%7.2f')); 
set(handles.oza2,'String',num2str(handles.oldzeros.a2,'%7.2f')); 
set(handles.oza4,'String',num2str(handles.oldzeros.a4,'%7.2f')); 





% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function updateFitTable(hObject,handles)
% fills the "Fit" part of the window
isfit = hasfield(handles,'fitteda4'); % a valid fit result exists
oldtable = get(handles.fittable,'Data');
try
    for i=1:length(handles.powderdata.file)
        fittable{i,1} = handles.powderdata.file{i};
        fittable{i,2} = num2str(handles.powderdata.fita4(i),'%7.2f');
        if isfit
            fittable{i,3} = num2str(handles.fitteda4(i),'%7.2f');
            fittable{i,4} = num2str(handles.powderdata.fita4(i) - handles.fitteda4(i),'%7.2f');
        else
            fittable{i,3} = '--';
            fittable{i,4} = '--';
        end
        if all(size(oldtable) >= [i,5]), fittable{i,5} = oldtable{i,5}; else fittable{i,5}=true; end 
    end
    set(handles.fittable,'Data',fittable);
    if hasfield(handles,'newzeros')
        set(handles.nza1,'String',num2str(handles.newzeros.a1,'%7.2f'));
        set(handles.nza2,'String',num2str(handles.newzeros.a2,'%7.2f'));
        set(handles.nza4,'String',num2str(handles.newzeros.a4,'%7.2f'));
    end   
catch, set(handles.fittable,'Data',{});
end

% --- Executes on button press in fitbutton.
function fitbutton_Callback(hObject, eventdata, handles) %#ok<*DEFNU>
% fitbutton is the function called by the button 'Fit'.
% the function performs the fit of the a2 and a4 offsets
a2ind = 1 - get(handles.fixa2,'Value');
fittabledata = get(handles.fittable,'Data');
if isempty(fittabledata) || size(fittabledata,1)<length(handles.powderdata.file)
    try fileindex = 1:length(handles.powderdata.file); catch, return; end
else
    for i= 1:size(fittabledata,1)
        fileindex(i) = fittabledata{i,5};
    end
    fileindex = find(fileindex);
end
%fit
[handles.newzeros,handles.fittedoffset,handles.fitteda4] = fita2a4zeros(handles.powderdata, handles.oldzeros, [a2ind, 1], fileindex);
handles.fileindex = fileindex;
handles.fitvarindex = [a2ind, 1]; %(save for later)
handles.reportfilename = savepowderreport(handles,fullfile(handles.filepath,'powcal_report')); % write report file
guidata(hObject,handles);
updateFitTable(hObject, handles);
savepowderreport(handles);




% --- Executes on mouse motion over figure - except title and menu.
function figure1_WindowButtonMotionFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% % generates the Help text for hovering mouse
pos = get(hObject, 'currentpoint'); % get mouse location on figure
x = pos(1); 
y = pos(2); % assign locations to x and y
if x>66 && x<562 && y>505 && y<527 % position of file input line
        set (handles.text10, 'string', ['Enter filenames. Enter multiple files like 0123[44:55], 0123[40,42,44:47,50], etc. (With or without path.) '...
            'Will try to download from server, if not found locally. Alternatively, click "Choose files" to open dialog window.']); 
elseif x>60 && x<206 && y>420 && y<456 % 'Evaluate data files' button
        set (handles.text10, 'string', 'Evaluate and fit selected files.'); 
elseif x>226 && x<628 && y>292 && y<460 % Scan fit results table
        set (handles.text10, 'string', ['Expected a4 value is deduced from data file. (If scan center does not match the d-value in PARAM.AS, '...
            'see message in matlab main window). "Measured" is the fitted center (Gauss+const.), '...
            'Chi2 refers to this fit. Values can be edited manually.']); 
elseif x>69 && x<391 && y>87 && y<243 % Fitting table for zeros
        set (handles.text10, 'string', ['Data used for za2/za4 fitting. "Measured" as above, "fitted" is the a4-value using the fitted offsets.' ...
            ' Untick reflections that are not to be used in za2/za4 determination.']);
elseif x>399 && x<474 && y>220 && y<243 % Check box fix a2
        set (handles.text10, 'string', 'Tick to keep za2 fixed and fit only za4.'); 
elseif x>397 && x<473 && y>90 && y<122   % Save report button
        set (handles.text10, 'string', 'Save the report (as in "view report") to a file.'); 

else
    set(handles.text10, 'string', ' ');
end
%set(handles.lbl_x, 'string', ['x loc:' num2str(x)]); % update text for x loc
%set(handles.lbl_y, 'string', ['y loc:' num2str(y)]); % update text for y loc 


% --- Executes when selected cell(s) is changed in uitable1.
function uitable1_CellSelectionCallback(hObject, eventdata, handles)
% hObject    handle to uitable1 (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) currently selecteds
% handles    structure with handles and user data (see GUIDATA)
if ~isempty(eventdata.Indices) && eventdata.Indices(2)==1
%     table = get(hObject,'Data');
    nplot(handles.filenames{eventdata.Indices(1)},'fit gauss','nolegend','nooutput','showfit');
end


% --- Executes on button press in filedlg.
function filedlg_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to filedlg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[FileName,PathName] = uigetfile('*','Select Data Files',handles.filepath,'Multiselect','on');
if isnumeric(FileName) && FileName==0, return; end
if ~iscell(FileName) 
    files = FileName; 
else  
    files = FileName{1};
    for n=2:length(FileName)
        files = [files,',',FileName{n}];
    end
end
if ~strcmpi([pwd,'\'],PathName) && ~strcmpi([pwd,'/'],PathName) && ~strcmpi(pwd,PathName)
    % files in other than the current folder -> add path name
    files=[fullfile(PathName,'['),files,']'];
end
set(handles.edit1,'String',files);
handles.filepath = PathName;
if isfield(handles,'fitteda4'), handles = rmfield(handles,'fitteda4'); end
guidata(hObject, handles);
run_Callback(hObject, eventdata, handles);


% --- Executes when entered data in editable cell(s) in uitable1.
function uitable1_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to uitable1 (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)
newtabledata = get(hObject,'Data');
if ~isempty(eventdata.Error)
    fprintf('%s\n',eventdata.Error);
    return;
% elseif isnan(tabledata(eventdata.Indices(1),eventdata.Indices(2)))
end
if eventdata.Indices(2)==2
    handles.powderdata.scancenter(eventdata.Indices(1)) = str2double(newtabledata{eventdata.Indices(1),2});
elseif eventdata.Indices(2)==3
    handles.powderdata.fita4(eventdata.Indices(1)) = str2double(newtabledata{eventdata.Indices(1),3});
else
    return;
end
handles.powderdata.chi2(eventdata.Indices(1)) = 0;
handles.powderdata.erra4(eventdata.Indices(1)) = 0.01;
guidata(hObject, handles);
updateTable(hObject, handles);    
updateFitTable(hObject, handles);



% --- Executes on button press in reportbutton.
function reportbutton_Callback(hObject, eventdata, handles)
% hObject    handle to reportbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if isfield(handles,'reportfilename'), edit(handles.reportfilename); end


% --- Executes on button press in savereportbutton.
function savereportbutton_Callback(hObject, eventdata, handles)
% hObject    handle to savereportbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[FileName,PathName] = uiputfile('*','Save Calibration Report',  ...
    fullfile(handles.filepath,['powder_',datestr(now,'yyyy-mm-dd_HH'),'h',datestr(now,'MM')]));
if FileName==0, return; end
savepowderreport(handles,fullfile(PathName,FileName));
handles.reportfilename = fullfile(PathName,FileName);
guidata(hObject, handles);


function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double
if isfield(handles,'fitteda4'), handles = rmfield(handles,'fitteda4'); end
guidata(hObject, handles);
