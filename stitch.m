function varargout = stitch(varargin)
%STITCH M-file for stitch.fig
%      STITCH, by itself, creates a new STITCH or raises the existing
%      singleton*.
%
%      H = STITCH returns the handle to a new STITCH or the handle to
%      the existing singleton*.
%
%      STITCH('Property','Value',...) creates a new STITCH using the
%      given property value pairs. Unrecognized properties are passed via
%      varargin to stitch_OpeningFcn.  This calling syntax produces a
%      warning when there is an existing singleton*.
%
%      STITCH('CALLBACK') and STITCH('CALLBACK',hObject,...) call the
%      local function named CALLBACK in STITCH.M with the given input
%      arguments.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help stitch

% Last Modified by GUIDE v2.5 29-Aug-2016 16:42:52

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @stitch_OpeningFcn, ...
                   'gui_OutputFcn',  @stitch_OutputFcn, ...
                   'gui_LayoutFcn',  [], ...
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


% --- Executes just before stitch is made visible.
function stitch_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   unrecognized PropertyName/PropertyValue pairs from the
%            command line (see VARARGIN)

% Choose default command line output for stitch
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes stitch wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% Clear workspace
evalin('base','clear,clc');

% Make sure all folders are added to search path
evalin('base','addpath(genpath(cd))'); 


% --- Outputs from this function are returned to the command line.
function varargout = stitch_OutputFcn(hObject, eventdata, handles)
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

% User selects Hexagon images
[cFile,strPath] = uigetfile('*.tif', strcat(['Select the Hexagon' ...
    ' .tif images:']),'MultiSelect','on');
cFile = sort(cFile);
assignin('base','cFile',cFile);
assignin('base','strPath',strPath);


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Close the GUI
close(handles.figure1)

% Run main script
evalin('base','mainStitch');
