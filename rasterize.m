function varargout = rasterize(varargin)
% RASTERIZE MATLAB code for rasterize.fig
%      RASTERIZE, by itself, creates a new RASTERIZE or raises the existing
%      singleton*.
%
%      H = RASTERIZE returns the handle to a new RASTERIZE or the handle to
%      the existing singleton*.
%
%      RASTERIZE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in RASTERIZE.M with the given input arguments.
%
%      RASTERIZE('Property','Value',...) creates a new RASTERIZE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before rasterize_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to rasterize_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help rasterize

% Last Modified by GUIDE v2.5 11-Oct-2016 13:43:30

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @rasterize_OpeningFcn, ...
                   'gui_OutputFcn',  @rasterize_OutputFcn, ...
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


% --- Executes just before rasterize is made visible.
function rasterize_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to rasterize (see VARARGIN)

% Choose default command line output for rasterize
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes rasterize wait for user response (see UIRESUME)
% uiwait(handles.figure1);

%Clear workspace
evalin('base','clear,clc');

%Make sure all folders are added to search path
evalin('base','addpath(genpath(cd))'); 


% --- Outputs from this function are returned to the command line.
function varargout = rasterize_OutputFcn(hObject, eventdata, handles) 
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

% User selects folder containing extracted Hexagon DEM files
strPath = uigetdir(cd,['Select top folder containing the extracted ' ...
    'Hexagon DEM files:']);
assignin('base','strWinPath',[strPath '\']);


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Close the GUI
close(handles.figure1)

% Run main script
evalin('base','mainRasterize');
