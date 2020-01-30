function varargout = extract(varargin)
% EXTRACT MATLAB code for extract.fig
%      EXTRACT, by itself, creates a new EXTRACT or raises the existing
%      singleton*.
%
%      H = EXTRACT returns the handle to a new EXTRACT or the handle to
%      the existing singleton*.
%
%      EXTRACT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in EXTRACT.M with the given input arguments.
%
%      EXTRACT('Property','Value',...) creates a new EXTRACT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before extract_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to extract_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help extract

% Last Modified by GUIDE v2.5 29-Jan-2020 22:55:19

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @extract_OpeningFcn, ...
                   'gui_OutputFcn',  @extract_OutputFcn, ...
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


% --- Executes just before extract is made visible.
function extract_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to extract (see VARARGIN)

% Choose default command line output for extract
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes extract wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% Clear workspace
evalin('base','clear,clc');

% Make sure all folders are added to search path
evalin('base','addpath(genpath(cd))'); 

% Set default selection for spatial resolution of output
set(handles.popupmenu1,'Value',2);

% Initialize parameters
iBlkSz = 3;
assignin('base','iBlkSz',iBlkSz);


% --- Outputs from this function are returned to the command line.
function varargout = extract_OutputFcn(hObject, eventdata, handles) 
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

% User selects the stitched Hexagon .mat files
[cHexFile,strHexPath] = uigetfile('*.mat', strcat(['Select the' ...
    ' stitched Hexagon .mat files:']),'MultiSelect','on');
if ~iscell(cHexFile)
    cHexFile = {cHexFile};
end
assignin('base','cHexFile',cHexFile);
assignin('base','strHexPath',strHexPath);


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Save spatial resolution of output
cRes = get(handles.popupmenu1,'String');
strRes = cRes{get(handles.popupmenu1,'Value')};
assignin('base','strRes',strRes);

% Close the GUI
close(handles.figure1)

% Run main script
evalin('base','mainExtract');


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
set(hObject,'String',{'Full';'1/2';'1/4';'1/8'});



function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double

% Get value
iBlkSz = str2double(get(hObject,'String'));
assignin('base','iBlkSz',iBlkSz);

% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
