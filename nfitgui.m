function varargout = nfitgui(varargin)
% NFITGUI MATLAB code for nfitgui.fig
%      NFITGUI, by itself, creates a new NFITGUI or raises the existing
%      singleton*.
%
%      H = NFITGUI returns the handle to a new NFITGUI or the handle to
%      the existing singleton*.
%
%      NFITGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in NFITGUI.M with the given input arguments.
%
%      NFITGUI('Property','Value',...) creates a new NFITGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before nfitgui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to nfitgui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help nfitgui

% Last Modified by GUIDE v2.5 06-Feb-2015 19:15:13

% if call without arguments (open new nfitgui), capture actual axes handle to pass to OpeningFcn
if isempty(varargin), varargin = { 'actualaxes', gca }; end

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @nfitgui_OpeningFcn, ...
                   'gui_OutputFcn',  @nfitgui_OutputFcn, ...
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




%%

% --- Executes just before nfitgui is made visible.
function nfitgui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to nfitgui (see VARARGIN)

% Choose default command line output for nfitgui
handles.output = hObject;
guidata(hObject, handles);

% Treat different cases of input (in varargin)
% - first argument is nfit object
% - first argument is axes handle
% - nothing provided, but gca handle inserted (in main function)
newobject = false;
if nargin>3 && isa(varargin{1},'nfit')      % is nfit object
    handles.nfitobj = varargin{1};
else
    axhandle = readinput('actualaxes',varargin);
    if isempty(axhandle) && nargin>3 && isgraphics(varargin{ 1 },'axes')  
        axhandle = varargin{1};
    end
    % check if axes have an nfit object attached (in guidata)
    if ishandle(axhandle), figdata = guidata(axhandle); else figdata = []; end
    if ~isempty(figdata) && isfield(figdata,'nfitobj')
        if iscell(figdata.nfitobj)
            for ii=1:length(figdata.nfitobj)
                if figdata.nfitobj{ii}.graphhandle == axhandle, handles.nfitobj = figdata.nfitobj{ii}; end
            end
        elseif figdata.nfitobj.graphhandle == axhandle, handles.nfitobj=figdata.nfitobj; 
        end
    else
        % try to use current axes to create new nfit object
        handles.nfitobj = nfit(axhandle);  newobject = true;
    end
end  
if isfield(handles,'nfitobj') && isempty(handles.nfitobj.xdata) %|| isempty(handles.nfitobj.graphhandle))
    fprintf('Error on opening nfitgui: No data found.\n'); 
    if newobject, delete(handles.nfitobj); end 
    delete(hObject); return;
elseif ~isfield(handles,'nfitobj')
	fprintf('Error on opening nfitgui: Could not set up fitting.\n'); delete(hObject); return;
end
% Initialize...
handles.selectedfunctionindex = [];
handles.selectedconstraintindex = [];
% Add listener
handles.listener = event.listener(handles.nfitobj,'nfitupdate',@(src,evnt)nfitgui_listenerfcn(hObject,src,evnt));
% Fill table
tableupdate(hObject,handles);
guidata(hObject, handles);


function tableupdate(hObject,handles)
% Fill entries in functions table
tabledata = {}; ptabledata = {};
for fct = 1:length(handles.nfitobj.fitfunction.components)
    try % function call 
        func = handles.nfitobj.fitfunction.components{fct};
        [~, ~, ~, tabledata{fct,1}] = func([],[]); 
    catch
        fprintf('Error: problem on calling function\n'); return; 
    end
    tabledata{fct,2} = handles.nfitobj.fitfunction.tag{fct}; %#ok<*AGROW>
end
set(handles.functable,'Data',tabledata);
% Fill parameter table
for prm = 1:numel(handles.nfitobj.parameters.values)
    ptabledata{prm,1} = ['[',num2str(handles.nfitobj.parameters.belongtocomponent(prm),'%2d'),']'];
    ptabledata{prm,2} = handles.nfitobj.parameters.names{prm};
    ptabledata{prm,3} = handles.nfitobj.parameters.values(prm);
    ptabledata{prm,4} = handles.nfitobj.parameters.errors(prm);
    ptabledata{prm,5} = logical(handles.nfitobj.parameters.fixed(prm));
end
set(handles.paramtable,'Data',ptabledata)
% Fill constraints table
set(handles.constrainttable,'Data',handles.nfitobj.constraints(:));


function nfitgui_listenerfcn(hObject,source,evnt) %#ok<*INUSD>
% listener function for notifications from nfit (--> react to changes)
tableupdate(hObject,guidata(hObject));



% --- Outputs from this function are returned to the command line.
function varargout = nfitgui_OutputFcn(hObject, eventdata, handles)  %#ok<*INUSL>
% Get default command line output from handles structure
if nargout>0, if isfield(handles,'output'), varargout{1} = handles.output; else varargout{1}=[]; end; end



% --- Button: remove Function
function pb_fcrem_Callback(hObject, eventdata, handles) %#ok<*DEFNU>
for ind = handles.selectedfunctionindex(:)'
    handles.nfitobj.rmfunc(ind);
end
handles.selectedfunctionindex = [];
% tableupdate(hObject,handles);
handles.nfitobj.plot;

% --- Button: add Function
function pb_fcadd_Callback(hObject, eventdata, handles)
% obtain list of available functions
flist = fitfunctions; 
flist = [{'','------- Standard functions -------'; 'gaussA.m','Gauss'; 'lorentzA.m','Lorentz'; ... 
                    'const.m','Constant'; 'linear.m','Linear'; 'poly5.m','Polynomial'; 'expo.m','Exponential'; ...
          '','---------- All functions ----------'}; flist];
% Dialog window    
[Selection,ok] = listdlg('ListString',flist(:,2), 'Selectionmode','single','Initialvalue',2, ...
                    'Name','Add function', 'Promptstring',{'Select a function:','(to define user functions, see template.m)'});
if ~ok || isempty(flist{Selection,1}), return; end
fname = flist{Selection,1};
handles.nfitobj.addfunc(str2func(fname(1:end-2)));
% tableupdate(hObject,handles);
handles.nfitobj.plot;


% --- Executes when selected cell(s) is changed in functable.
function functable_CellSelectionCallback(hObject, eventdata, handles)
% simply store the selected rows in the handles structure
handles.selectedfunctionindex = eventdata.Indices(:,1);
guidata(hObject,handles);


% --- Executes when entered data in editable cell(s) in paramtable.
function paramtable_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to paramtable (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)
if eventdata.Indices(2)==5  % check box
    handles.nfitobj.parameters.fixed(eventdata.Indices(1)) = ~handles.nfitobj.parameters.fixed(eventdata.Indices(1));
elseif eventdata.Indices(2)==3 && ~isempty(eventdata.NewData) 
    handles.nfitobj.setparam(eventdata.Indices(1),eventdata.NewData); 
%     handles.nfitobj.plot;
end


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% handles.nfitobj.guihandle = [];
if isfield(handles,'listener'), delete(handles.listener); end
delete(hObject); % closes the figure


% --- Executes on button press in fitbutton.
function fitbutton_Callback(hObject, eventdata, handles)
handles.nfitobj.fit;


% --- Executes on button press in pb_consrem.
function pb_consrem_Callback(hObject, eventdata, handles)
% remove selected constraints from list
handles.nfitobj.rmconstraint(handles.selectedconstraintindex);
handles.selectedconstraintindex = [];
tableupdate(hObject,handles);


% --- Executes on button press in pb_consadd.
function pb_consadd_Callback(hObject, eventdata, handles)
% add constraint to list
str = inputdlg('Enter a linear constraint (e.g. p3=2*p2)','Add constraint',1);
if ~isempty(str{1}), handles.nfitobj.addconstraint(str{1}); end
tableupdate(hObject,handles);


% --- Executes on button press in pb_consedit.
function pb_consedit_Callback(hObject, eventdata, handles)  
% Edit the selected constraint
ind = handles.selectedconstraintindex;
if numel(ind)~=1, msgbox('Select one constraint from the list','modal'); return; end
str = inputdlg('Enter a linear constraint (e.g. p3=2*p2)','Edit constraint',1,handles.nfitobj.constraints(ind));
if ~isempty(str{1}), handles.nfitobj.constraints{ind} = str{1}; end
tableupdate(hObject,handles);


% --- Executes when selected cell(s) is changed in constrainttable.
function constrainttable_CellSelectionCallback(hObject, eventdata, handles)
% simply store the selected rows in the handles structure
handles.selectedconstraintindex = eventdata.Indices(:,1);
guidata(hObject,handles);
