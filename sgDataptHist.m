%%
function varargout = sgDataptHist(varargin)
% SGDATAPTHIST MATLAB code for sgDataptHist.fig
%      SGDATAPTHIST, by itself, creates a new SGDATAPTHIST or raises the existing
%      singleton*.
%
%      H = SGDATAPTHIST returns the handle to a new SGDATAPTHIST or the handle to
%      the existing singleton*.
%
%      SGDATAPTHIST('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SGDATAPTHIST.M with the given input arguments.
%
%      SGDATAPTHIST('Property','Value',...) creates a new SGDATAPTHIST or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before sgDataptHist_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to sgDataptHist_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help sgDataptHist

% Last Modified by GUIDE v2.5 04-Dec-2019 16:45:49

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @sgDataptHist_OpeningFcn, ...
                   'gui_OutputFcn',  @sgDataptHist_OutputFcn, ...
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


% --- Executes just before sgDataptHist is made visible.
function sgDataptHist_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to sgDataptHist (see VARARGIN)

% Choose default command line output for sgDataptHist
handles.output = hObject;
% Update handles structure
guidata(hObject, handles);
% UIWAIT makes sgDataptHist wait for user response (see UIRESUME)
% uiwait(handles.figure1);
global im0 G

refreshLst(handles)
val = find(G.sgLst(:,1) == im0.gui.segmentID);
set(handles.sgmntLst, 'value', val)
updateHistogram(handles, im0.gui.segmentID)




% --- Outputs from this function are returned to the command line.
function varargout = sgDataptHist_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function refreshLst(handles)
global im0 G
sgmnts = unique(im0.selectedEdgesList(1,:));
numOfEdge = ones(1, length(sgmnts));
for i = 1:length(sgmnts)
    sg = sgmnts(i);
    numOfEdge(i) = length(find(im0.selectedEdgesList(1,:) == sg));
end

sgLst = sortrows([sgmnts; numOfEdge]');
G.sgLst = sgLst;
set(handles.sgmntLst, 'string',num2str(sgLst, '%8d'))


function dataPt = sgDatapt(sgmntID)
global im0
edges = im0.selectedEdgesList(2, im0.selectedEdgesList(1,:) == sgmntID);
dataPt = [];
for e = edges
    dataPt = [dataPt, getDataPt(e, 1)];
end


function updateHistogram(handles, sgmntID)
dataPt = sgDatapt(sgmntID);
histogram(handles.hist, dataPt, 50, 'Orientation', 'horizontal')
enableDefaultInteractivity(handles.hist)




% --- Executes on selection change in sgmntLst.
function sgmntLst_Callback(hObject, eventdata, handles)
% hObject    handle to sgmntLst (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global G
sgmntID = G.sgLst(get(hObject,'Value'), 1);
updateHistogram(handles, sgmntID)
% Hints: contents = cellstr(get(hObject,'String')) returns sgmntLst contents as cell array
%        contents{get(hObject,'Value')} returns selected item from sgmntLst


% --- Executes during object creation, after setting all properties.
function sgmntLst_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sgmntLst (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in closeButton.
function closeButton_Callback(hObject, eventdata, handles)
delete(handles.figure1)
% hObject    handle to closeButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- 
function dataPt = getDataPt(edgeID, method)

%%- method can be 'all' which would return all datapoints in the cylinder
%%- or can be 'zeroCut' which applies a threshold at zero.

%%
global im0;

r = im0.cylDia(edgeID)/2;

pp0 = im0.cylNodes(edgeID, 1:3);
pp1 = im0.cylNodes(edgeID, 4:6);

p0 = [pp0(2), pp0(1), pp0(3)];
p1 = [pp1(2), pp1(1), pp1(3)];

v = p1-p0;       %edge vector
vScld = (v.*im0.Hvox);  %scaled edge vector
vScldMag = norm(vScld);

x0 = min(p0(1),p1(1))-r;
x1 = max(p0(1),p1(1))+r;
y0 = min(p0(2),p1(2))-r;
y1 = max(p0(2),p1(2))+r;
z0 = min(p0(3),p1(3))-r;
z1 = max(p0(3),p1(3))+r;
c0 = floor([x0,y0,z0]);
c1 = ceil([x1,y1,z1]);
Vsize = c1-c0;          %Volume Size
Vsize = max(Vsize(:));
Vsize = [Vsize Vsize Vsize];
c1 = c0 + Vsize;

p0 = p0 - c0;
dataPt = [];

if c1(3) <= size(im0.OCT_Flow,3) && c0(3) > 0
    flowVol = im0.OCT_Flow(c0(1):c1(1), c0(2):c1(2), c0(3):c1(3)); %flow Volume
    
    [X, Y, Z] = meshgrid(0:Vsize(1), 0:Vsize(2), 0:Vsize(3));
    X = (X - p0(1))* im0.Hvox(1);
    Y = (Y - p0(2))* im0.Hvox(2);
    Z = (Z - p0(3))* im0.Hvox(3);
    
    d = sqrt((Y*vScld(3) - Z*vScld(2)).^2 + (Z*vScld(1) - X*vScld(3)).^2 +...
        (X*vScld(2) - Y*vScld(1)).^2) ./ vScldMag;
    v0Prj = (X*vScld(1) + Y*vScld(2) + Z*vScld(3))/vScldMag;
    
    dataPt = flowVol(d<=r & v0Prj<=vScldMag & v0Prj>=0);
    dataPt = dataPt(dataPt ~= 0);
    
    if method == 0  %- make a cuoff at zero
        dataMean = mean(dataPt);
        if dataMean > 0
            dataPt = dataPt(dataPt > 0);
        else
            dataPt = dataPt(dataPt < 0);
        end
    end
    dataPt = dataPt';
    vScld = (im0.nodePos_um(im0.nodeEdges(edgeID, 2), :) - im0.nodePos_um(im0.nodeEdges(edgeID, 1), :));
    vScldMag = norm(vScld);
    cosPhZ = dot([0,0,1],vScld)/vScldMag;
    dataPt = dataPt/cosPhZ;
end


% --- Executes on button press in refreshButton.
function refreshButton_Callback(hObject, eventdata, handles)
refreshLst(handles)
% hObject    handle to refreshButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
delete(hObject);
