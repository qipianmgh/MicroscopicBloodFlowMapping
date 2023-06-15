function varargout = cutoffMethodDialog(varargin)
% CUTOFFMETHODDIALOG MATLAB code for cutoffMethodDialog.fig
%      CUTOFFMETHODDIALOG, by itself, creates a new CUTOFFMETHODDIALOG or raises the existing
%      singleton*.
%
%      H = CUTOFFMETHODDIALOG returns the handle to a new CUTOFFMETHODDIALOG or the handle to
%      the existing singleton*.
%
%      CUTOFFMETHODDIALOG('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CUTOFFMETHODDIALOG.M with the given input arguments.
%
%      CUTOFFMETHODDIALOG('Property','Value',...) creates a new CUTOFFMETHODDIALOG or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before cutoffMethodDialog_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to cutoffMethodDialog_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help cutoffMethodDialog

% Last Modified by GUIDE v2.5 03-Dec-2019 17:54:33

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @cutoffMethodDialog_OpeningFcn, ...
                   'gui_OutputFcn',  @cutoffMethodDialog_OutputFcn, ...
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


% --- Executes just before cutoffMethodDialog is made visible.
function cutoffMethodDialog_OpeningFcn(hObject, eventdata, handles, varargin)
handles.output = hObject;
guidata(hObject, handles);

global im0

if isfield(im0, 'cuoffMethd')
    if isfield(im0.cuoffMethd, 'out')
        ctM = im0.cuoffMethd.out;
        l = size(ctM, 1);
        for itr = 1:l
            set(eval(['handles.segIdMethod_0' num2str(itr)]), 'Value', (ctM(itr,1) + 2))
            set(eval(['handles.diamLrLim_0' num2str(itr)]), 'String',num2str(ctM(itr,2)))
            set(eval(['handles.diamUprLim_0' num2str(itr)]), 'String',num2str(ctM(itr,3)))
            set(eval(['handles.cutoffMethod_0' num2str(itr)]), 'Value', (ctM(itr,4) + 2))
            set(eval(['handles.fixedCutoff_0' num2str(itr)]), 'String',num2str(ctM(itr,5)))
        end
    end
end

im0.cuoffMethd.canceled = true;



% --- Outputs from this function are returned to the command line.
function varargout = cutoffMethodDialog_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function diamLrLim_01_Callback(hObject, eventdata, handles)
% hObject    handle to diamLrLim_01 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of diamLrLim_01 as text
%        str2double(get(hObject,'String')) returns contents of diamLrLim_01 as a double


% --- Executes during object creation, after setting all properties.
function diamLrLim_01_CreateFcn(hObject, eventdata, handles)
% hObject    handle to diamLrLim_01 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function diamUprLim_01_Callback(hObject, eventdata, handles)
% hObject    handle to diamUprLim_01 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of diamUprLim_01 as text
%        str2double(get(hObject,'String')) returns contents of diamUprLim_01 as a double


% --- Executes during object creation, after setting all properties.
function diamUprLim_01_CreateFcn(hObject, eventdata, handles)
% hObject    handle to diamUprLim_01 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in cutoffMethod_01.
function cutoffMethod_01_Callback(hObject, eventdata, handles)
mthd = get(hObject,'Value');
if mthd == 5
    set(handles.fixedCutoff_01, 'enable', 'on')
else
    set(handles.fixedCutoff_01, 'enable', 'off')
end

% hObject    handle to cutoffMethod_01 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns cutoffMethod_01 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from cutoffMethod_01


% --- Executes during object creation, after setting all properties.
function cutoffMethod_01_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cutoffMethod_01 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function fixedCutoff_01_Callback(hObject, eventdata, handles)
% hObject    handle to fixedCutoff_01 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fixedCutoff_01 as text
%        str2double(get(hObject,'String')) returns contents of fixedCutoff_01 as a double


% --- Executes during object creation, after setting all properties.
function fixedCutoff_01_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fixedCutoff_01 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function diamLrLim_02_Callback(hObject, eventdata, handles)
% hObject    handle to diamLrLim_02 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of diamLrLim_02 as text
%        str2double(get(hObject,'String')) returns contents of diamLrLim_02 as a double


% --- Executes during object creation, after setting all properties.
function diamLrLim_02_CreateFcn(hObject, eventdata, handles)
% hObject    handle to diamLrLim_02 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function diamUprLim_02_Callback(hObject, eventdata, handles)
% hObject    handle to diamUprLim_02 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of diamUprLim_02 as text
%        str2double(get(hObject,'String')) returns contents of diamUprLim_02 as a double


% --- Executes during object creation, after setting all properties.
function diamUprLim_02_CreateFcn(hObject, eventdata, handles)
% hObject    handle to diamUprLim_02 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in cutoffMethod_02.
function cutoffMethod_02_Callback(hObject, eventdata, handles)
mthd = get(hObject,'Value');
if mthd == 5
    set(handles.fixedCutoff_02, 'enable', 'on')
else
    set(handles.fixedCutoff_02, 'enable', 'off')
end
% hObject    handle to cutoffMethod_02 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns cutoffMethod_02 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from cutoffMethod_02


% --- Executes during object creation, after setting all properties.
function cutoffMethod_02_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cutoffMethod_02 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function fixedCutoff_02_Callback(hObject, eventdata, handles)
% hObject    handle to fixedCutoff_02 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fixedCutoff_02 as text
%        str2double(get(hObject,'String')) returns contents of fixedCutoff_02 as a double


% --- Executes during object creation, after setting all properties.
function fixedCutoff_02_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fixedCutoff_02 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in segIdMethod_01.
function segIdMethod_01_Callback(hObject, eventdata, handles)

mthd = get(hObject,'Value');
if mthd == 3
    set(handles.diamLrLim_01, 'enable', 'on')
    set(handles.diamUprLim_01, 'enable', 'on')
else
    set(handles.diamLrLim_01, 'enable', 'off')
    set(handles.diamUprLim_01, 'enable', 'off')
end


% hObject    handle to segIdMethod_01 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns segIdMethod_01 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from segIdMethod_01


% --- Executes during object creation, after setting all properties.
function segIdMethod_01_CreateFcn(hObject, eventdata, handles)
% hObject    handle to segIdMethod_01 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in segIdMethod_02.
function segIdMethod_02_Callback(hObject, eventdata, handles)
mthd = get(hObject,'Value');
if mthd == 3
    set(handles.diamLrLim_02, 'enable', 'on')
    set(handles.diamUprLim_02, 'enable', 'on')
else
    set(handles.diamLrLim_02, 'enable', 'off')
    set(handles.diamUprLim_02, 'enable', 'off')
end
% hObject    handle to segIdMethod_02 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns segIdMethod_02 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from segIdMethod_02


% --- Executes during object creation, after setting all properties.
function segIdMethod_02_CreateFcn(hObject, eventdata, handles)
% hObject    handle to segIdMethod_02 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function diamLrLim_03_Callback(hObject, eventdata, handles)
% hObject    handle to diamLrLim_03 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of diamLrLim_03 as text
%        str2double(get(hObject,'String')) returns contents of diamLrLim_03 as a double


% --- Executes during object creation, after setting all properties.
function diamLrLim_03_CreateFcn(hObject, eventdata, handles)
% hObject    handle to diamLrLim_03 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function diamUprLim_03_Callback(hObject, eventdata, handles)
% hObject    handle to diamUprLim_03 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of diamUprLim_03 as text
%        str2double(get(hObject,'String')) returns contents of diamUprLim_03 as a double


% --- Executes during object creation, after setting all properties.
function diamUprLim_03_CreateFcn(hObject, eventdata, handles)
% hObject    handle to diamUprLim_03 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in cutoffMethod_03.
function cutoffMethod_03_Callback(hObject, eventdata, handles)
mthd = get(hObject,'Value');
if mthd == 5
    set(handles.fixedCutoff_03, 'enable', 'on')
else
    set(handles.fixedCutoff_03, 'enable', 'off')
end
% hObject    handle to cutoffMethod_03 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns cutoffMethod_03 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from cutoffMethod_03


% --- Executes during object creation, after setting all properties.
function cutoffMethod_03_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cutoffMethod_03 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function fixedCutoff_03_Callback(hObject, eventdata, handles)
% hObject    handle to fixedCutoff_03 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fixedCutoff_03 as text
%        str2double(get(hObject,'String')) returns contents of fixedCutoff_03 as a double


% --- Executes during object creation, after setting all properties.
function fixedCutoff_03_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fixedCutoff_03 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in segIdMethod_03.
function segIdMethod_03_Callback(hObject, eventdata, handles)
mthd = get(hObject,'Value');
if mthd == 3
    set(handles.diamLrLim_03, 'enable', 'on')
    set(handles.diamUprLim_03, 'enable', 'on')
else
    set(handles.diamLrLim_03, 'enable', 'off')
    set(handles.diamUprLim_03, 'enable', 'off')
end
% hObject    handle to segIdMethod_03 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns segIdMethod_03 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from segIdMethod_03


% --- Executes during object creation, after setting all properties.
function segIdMethod_03_CreateFcn(hObject, eventdata, handles)
% hObject    handle to segIdMethod_03 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in cancelButton.
function cancelButton_Callback(hObject, eventdata, handles)
global im0
answer = questdlg('are you sure you want to cancel this operation?', ...
    'confirm cancelation', 'Yes','No','No');
switch answer
    case 'Yes'
        im0.cuoffMethd.canceled = true;
        delete(handles.figure1)
    case 'No'
        return
end

% hObject    handle to cancelButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in OKButton.
function OKButton_Callback(hObject, eventdata, handles)
global im0
numOfOp = 3;

out = zeros(numOfOp, 5);

for itr = 1:numOfOp
    i = get(eval(['handles.segIdMethod_0' num2str(itr)]),'Value') - 2;
    ii = str2double(get(eval(['handles.diamLrLim_0' num2str(itr)]),'String'));
    iii = str2double(get(eval(['handles.diamUprLim_0' num2str(itr)]),'String'));
    iv = get(eval(['handles.cutoffMethod_0' num2str(itr)]),'Value') - 2;
    v = str2double(get(eval(['handles.fixedCutoff_0' num2str(itr)]),'String'));
    
    out(itr, :) = [i, ii, iii, iv, v];
end

lgc = (out(:,1) < 1 | out(:,4) < 1)';
if all(lgc)
    msgbox({'Invalid inputI',...
        'Please complete all input elements for atleast one row'},...
        'Error','error');
    return
end

out(lgc, :) = [];
im0.cuoffMethd.out = out;
im0.cuoffMethd.canceled = false;

delete(handles.figure1)


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
cancelButton_Callback(hObject, eventdata, handles)
% Hint: delete(hObject) closes the figure
% delete(hObject);
