function varargout = stats(varargin)
% STATS MATLAB code for stats.fig
%      STATS, by itself, creates a new STATS or raises the existing
%      singleton*.
%
%      H = STATS returns the handle to a new STATS or the handle to
%      the existing singleton*.
%
%      STATS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in STATS.M with the given input arguments.
%
%      STATS('Property','Value',...) creates a new STATS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before stats_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to stats_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help stats

% Last Modified by GUIDE v2.5 04-Sep-2019 12:01:53

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @stats_OpeningFcn, ...
                   'gui_OutputFcn',  @stats_OutputFcn, ...
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


% --- Executes just before stats is made visible.
function stats_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to stats (see VARARGIN)

% Choose default command line output for stats
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

updateAxes(handles)





% UIWAIT makes stats wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = stats_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function updateAxes(handles)
global im0
data = table2array(im0.data);
sgmnts = data(:,1)';
dias = im0.segDiam(sgmnts);

l01 = str2double(get(handles.lwrD_hist01,'String'));
u01 = str2double(get(handles.uprD_hist01,'String'));
l02 = str2double(get(handles.lwrD_hist02,'String'));
u02 = str2double(get(handles.uprD_hist02,'String'));
lmt_01 = [l01 u01];
lmt_02 = [l02 u02];

%%-- histogram 01
nbins = str2double(get(handles.bin_hist01,'String'));
flwV = data((dias >= lmt_01(1) & dias <= lmt_01(2)),3);
absChk = get(handles.absChk_hist01,'Value');
if absChk
    flwV = abs(flwV);
end
histogram(handles.hist01, flwV, nbins)
xlabel(handles.hist01, 'mm/s')
handles.hist01.Toolbar.Visible = 'on';
enableDefaultInteractivity(handles.hist01)

%%-- histogram 02
nbins = str2double(get(handles.bin_hist02,'String'));
flwV = data((dias >= lmt_02(1) & dias <= lmt_02(2)),3);
absChk = get(handles.absChk_hist02,'Value');
if absChk
    flwV = abs(flwV);
end
histogram(handles.hist02, flwV, nbins)
xlabel(handles.hist02, 'mm/s')
handles.hist02.Toolbar.Visible = 'on';
enableDefaultInteractivity(handles.hist02)






function lwrD_hist01_Callback(hObject, eventdata, handles)
% hObject    handle to lwrD_hist01 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
updateAxes(handles)
% Hints: get(hObject,'String') returns contents of lwrD_hist01 as text
%        str2double(get(hObject,'String')) returns contents of lwrD_hist01 as a double


% --- Executes during object creation, after setting all properties.
function lwrD_hist01_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lwrD_hist01 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function lwrD_hist02_Callback(hObject, eventdata, handles)
% hObject    handle to lwrD_hist02 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
updateAxes(handles)
% Hints: get(hObject,'String') returns contents of lwrD_hist02 as text
%        str2double(get(hObject,'String')) returns contents of lwrD_hist02 as a double


% --- Executes during object creation, after setting all properties.
function lwrD_hist02_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lwrD_hist02 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function uprD_hist01_Callback(hObject, eventdata, handles)
% hObject    handle to uprD_hist01 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
updateAxes(handles)
% Hints: get(hObject,'String') returns contents of uprD_hist01 as text
%        str2double(get(hObject,'String')) returns contents of uprD_hist01 as a double


% --- Executes during object creation, after setting all properties.
function uprD_hist01_CreateFcn(hObject, eventdata, handles)
% hObject    handle to uprD_hist01 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function uprD_hist02_Callback(hObject, eventdata, handles)
% hObject    handle to uprD_hist02 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
updateAxes(handles)
% Hints: get(hObject,'String') returns contents of uprD_hist02 as text
%        str2double(get(hObject,'String')) returns contents of uprD_hist02 as a double


% --- Executes during object creation, after setting all properties.
function uprD_hist02_CreateFcn(hObject, eventdata, handles)
% hObject    handle to uprD_hist02 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function bin_hist01_Callback(hObject, eventdata, handles)
% hObject    handle to bin_hist01 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
updateAxes(handles)
% Hints: get(hObject,'String') returns contents of bin_hist01 as text
%        str2double(get(hObject,'String')) returns contents of bin_hist01 as a double


% --- Executes during object creation, after setting all properties.
function bin_hist01_CreateFcn(hObject, eventdata, handles)
% hObject    handle to bin_hist01 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function bin_hist02_Callback(hObject, eventdata, handles)
% hObject    handle to bin_hist02 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
updateAxes(handles)
% Hints: get(hObject,'String') returns contents of bin_hist02 as text
%        str2double(get(hObject,'String')) returns contents of bin_hist02 as a double


% --- Executes during object creation, after setting all properties.
function bin_hist02_CreateFcn(hObject, eventdata, handles)
% hObject    handle to bin_hist02 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in absChk_hist01.
function absChk_hist01_Callback(hObject, eventdata, handles)
% hObject    handle to absChk_hist01 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
updateAxes(handles)
% Hint: get(hObject,'Value') returns toggle state of absChk_hist01


% --- Executes on button press in absChk_hist02.
function absChk_hist02_Callback(hObject, eventdata, handles)
% hObject    handle to absChk_hist02 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
updateAxes(handles)
% Hint: get(hObject,'Value') returns toggle state of absChk_hist02
