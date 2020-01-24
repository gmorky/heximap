function varargout = georef(varargin)
% GEOREF MATLAB code for georef.fig
%      GEOREF, by itself, creates a new GEOREF or raises the existing
%      singleton*.
%
%      H = GEOREF returns the handle to a new GEOREF or the handle to
%      the existing singleton*.
%
%      GEOREF('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GEOREF.M with the given input arguments.
%
%      GEOREF('Property','Value',...) creates a new GEOREF or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before georef_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to georef_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help georef

% Last Modified by GUIDE v2.5 26-Aug-2019 12:55:20

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @georef_OpeningFcn, ...
                   'gui_OutputFcn',  @georef_OutputFcn, ...
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


% --- Executes just before georef is made visible.
function georef_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to georef (see VARARGIN)

% Choose default command line output for georef
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes georef wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% Clear workspace
evalin('base','clear,clc');

% Initialize checkboxes
lVis = false;
assignin('base','lVis',lVis);
lPoly = false;
assignin('base','lPoly',lPoly);

% Make sure all folders are added to search path
evalin('base','addpath(genpath(cd))'); 


% --- Outputs from this function are returned to the command line.
function varargout = georef_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% User selects reference DEM .tif file
[strFile,strPath] = uigetfile('*.tif', ...
    'Select the reference DEM .tif file:');
assignin('base','strRef',[strPath strFile]);


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% User selects folder containing unstable terrain shapefiles
strPath = uigetdir(cd,['Select top folder containing the unstable ' ...
    'terrain .shp files:']);
assignin('base','strShpPath',[strPath '\']);


% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% User selects folder containing extracted Hexagon DEM files
strPath = uigetdir(cd,['Select top folder containing the extracted ' ...
    'Hexagon DEM files:']);
assignin('base','strWinPath',[strPath '\']);


% --- Executes on button press in pushbutton2.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Save control point selection method
cRes = get(handles.popupmenu1,'String');
strCP = cRes{get(handles.popupmenu1,'Value')};
assignin('base','strCP',strCP);

% Close the GUI
close(handles.figure1)

% Run main script
evalin('base','mainGeoref');


% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1


% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% Populate the menu
set(hObject,'String',{'automatic';'manual'});


% --- Executes on button press in checkbox1.
function checkbox1_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox1

% Get toggle state of checkbox
lVis = get(hObject,'Value');
assignin('base','lVis',lVis);


% --- Executes on button press in checkbox2.
function checkbox2_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox2

% Get toggle state of checkbox
lPoly = get(hObject,'Value');
assignin('base','lPoly',lPoly);
