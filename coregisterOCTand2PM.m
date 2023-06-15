function varargout = coregisterOCTand2PM_CarolineResult(varargin)
%      COREGISTEROCTAND2PM M-file for coregisterOCTand2PM.fig
%      COREGISTEROCTAND2PM, by itself, creates a new COREGISTEROCTAND2PM or raises the existing
%      singleton*.
%
%      H = COREGISTEROCTAND2PM returns the handle to a new COREGISTEROCTAND2PM or the handle to
%      the existing singleton*.
%
%      COREGISTEROCTAND2PM('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in COREGISTEROCTAND2PM.M with the given input arguments.
%
%      COREGISTEROCTAND2PM('Property','Value',...) creates a new COREGISTEROCTAND2PM or raises the
%      existing singleton*. Starting from the left, property value pairs are
%      applied to the GUI before coregisterOCTand2PM_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to coregisterOCTand2PM_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help coregisterOCTand2PM

% Last Modified by GUIDE v2.5 07-Aug-2018 13:42:39

% Begin initialization code - DO NOT EDIT

gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @coregisterOCTand2PM_OpeningFcn, ...
                   'gui_OutputFcn',  @coregisterOCTand2PM_OutputFcn, ...
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


% --- Executes just before coregisterOCTand2PM is made visible.
function coregisterOCTand2PM_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to coregisterOCTand2PM (see VARARGIN)
% Choose default command line output for coregisterOCTand2PM
handles.output = hObject;

% Update handles structure

guidata(hObject, handles);

clear global im;
global im

quest = 'Loading Co-registered Data?';
answer = questdlg(quest,'Data type','Yes','No','No');
    
 switch answer
 case 'Yes'
        answer_flag=1;
 case 'No'
        answer_flag=0;
 end
    
  if answer_flag
    
    [fooFname, fooPath] = uigetfile('*.mat','Select MAT file of co-registered data');
    infname = [fooPath fooFname];
    load(infname);
    if ~isfield(im,'coreg_flag')
        im.coreg_flag=1;  %%% coregistration performed
    end
    
    if ~isfield(im,'flow_flag')
        im.flow_flag=0;   %%% flow data not loaded
    end
        
  else
    
    % load 'It' - target, 2PM stack
    im.coreg_flag=0;
    im.flow_flag=0;
    
    [fooFname, fooPath] = uigetfile('*.mat','Select MAT 3D file of TPM data');
    infname = [fooPath fooFname];
    load(infname);
    x=whos('-file',infname);
    eval(['im.It = ' x.name ';' 'clear ' x.name ';']); % it is Z, Y, X matrix
    % create Y,X,Z matrix...
    im.It = shiftdim(im.It,1);
    
    if max(im.It(:))<65535
       im.It=im.It/max(im.It(:))*65535; 
    else  
    end
    
    [fooFname, fooPath] = uigetfile('*.mat','Select MAT 3D file of OCT data');
    infname = [fooPath fooFname];
    foo = load(infname);
    
    if isfield(foo,'OCT_DLS_angio')
           [Nz Nx Ny] = size(foo.OCT_DLS_angio);
           im.Im=single(zeros(Nz,Nx,Ny));
           im.Im=foo.OCT_DLS_angio;
                  
           for ii = 1:Nz,
              foo1 = squeeze(im.Im(ii,:,:));
              im.Im(ii,:,:) = rot90(foo1,2);
           end
                  
           im.Im = shiftdim(im.Im,1);
           clear foo1;
           clear foo; 
        
    elseif isfield(foo,'OCT_Doppler_angio')
           [Nz Nx Ny] = size(foo.OCT_Doppler_angio);
           im.Im=single(zeros(Nz,Nx,Ny));
           im.Im=foo.OCT_Doppler_angio;
           
             for ii = 1:Nz,
              foo1 = squeeze(im.Im(ii,:,:));
              im.Im(ii,:,:) = rot90(foo1,2);
             end
                  
           im.Im = shiftdim(im.Im,1);
           clear foo1;
           clear foo; 
            
     elseif isfield(foo,'OCT_DLS_flow')
               [Nz Nx Ny] = size(foo.OCT_DLS_flow);
               im.Im=single(zeros(Nz,Nx,Ny));
               im.Im=foo.OCT_DLS_flow;
                for ii = 1:Nz,
                    foo1 = squeeze(im.Im(ii,:,:));
                    im.Im(ii,:,:) = rot90(foo1,2);
                end
                  
               im.Im = shiftdim(im.Im,1);
               clear foo1;
               clear foo; 
                   
      elseif isfield(foo,'OCT_Doppler_flow')
          
      elseif isfield(foo,'OCT_Coreg_flow')
             
             im.flow_flag = 1;
             [Nz Nx Ny] = size(foo.OCT_Coreg_flow);
             im.Im = single(zeros(Nz,Nx,Ny));
             im.Im = foo.OCT_Coreg_flow;
             clear foo; 
                  
    end
  
    im.Im0=im.Im;   %  save the original OCT stack in case log operation is needed     
    if min(im.Im(:))>=0
       im.Im=log(im.Im+1);
    else    
    end

    im.CrangeM = [min(im.Im(:)) max(im.Im(:))];
    im.CrangeT = [min(im.It(:)) max(im.It(:))];
    
    set(handles.editAxes2Crange,'string',sprintf('%d ',im.CrangeM) );
    set(handles.editAxes1Crange,'string',sprintf('%d ',im.CrangeT) );
    im.h = 10;
    im.h_moveable = 10;
    im.ht = [1 1 1];
    im.hm = [1 1 1];
    set(handles.editVoxelSizeTarget,'string',sprintf('%.1f ',im.ht));
    set(handles.editVoxelSizeMoveable,'string',sprintf('%.1f ',im.hm));
    set(handles.editWid,'string',sprintf('%d',im.h));
    set(handles.editWid_moveable,'string',sprintf('%d',im.h_moveable));
    
    im.flagOverlay = 0;
    set(handles.checkboxOverlay,'value',im.flagOverlay);
 
    im.pos = [1 1 1];
    im.pos_moveable = [1 1 1];
    
    im.Z2fit(1) = im.ht(3);
    im.Z2fit(2) = 0;
    
    im.Z2fit_moveable(1) = im.hm(3);
    im.Z2fit_moveable(2) = 0;
    
    [nyt nxt nzt] = size(im.It);
    [nym nxm nzm] = size(im.Im);
    
    im.posMax = [nxt nyt nzt];
    im.posMax_moveable = [nxm nym nzm];

    im.MarkerNum = 0;
    im.MarkerState = 0;  
    
    % 0 - no active marker
    % 1 - marker placed in axes 1, waiting for
    %     corresponding marker in axes 2
    % 2 - marker placed in axes 2, waiting for
    %     corresponding marker in axes 1
    
    im.MarkerPos1 = [];
    im.MarkerPos2 = [];
    im.T = eye(4,4);
    
    im.flagButton = 1;
    im.flagView = 3;
    im.MarkerPoint = 0;
    
    im.ROI_MarkerState = 0;
    im.ROI_MarkerPos1 = [];
    im.ROI_MarkerPos2 = [];
    im.ROI_MarkerNum = 0;
end
updateImages( handles )



function varargout = coregisterOCTand2PM_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Get default command line output from handles structure
varargout{1} = handles.output;


function pushbuttonUp_Callback(hObject, eventdata, handles)
global im
fv = im.flagView; 

im.pos(fv)=max(im.pos(fv)-1,1);

if im.flagOverlay==1 & isfield(im,'reg_method')
im.pos_moveable= im.pos;
im.h_moveable= im.h;
im.hm = im.ht;
end

updateImages( handles );



function pushbuttonUp10_Callback(hObject, eventdata, handles)
global im
fv = im.flagView; 


im.pos(fv)=max(im.pos(fv)-10,1);

if im.flagOverlay==1 & isfield(im,'reg_method')
im.pos_moveable= im.pos;
im.h_moveable= im.h;
im.hm = im.ht;
end

updateImages( handles );


function pushbuttonDown_Callback(hObject, eventdata, handles)
global im
fv = im.flagView; 
im.pos(fv)=min(im.pos(fv)+1,im.posMax(fv));

if im.flagOverlay==1 & isfield(im,'reg_method')
im.pos_moveable=im.pos;
im.h_moveable=im.h;
im.hm = im.ht;
end

updateImages( handles );


function pushbuttonDown10_Callback(hObject, eventdata, handles)
global im
fv = im.flagView; 
im.pos(fv)= min(im.pos(fv)+10,im.posMax(fv));

if im.flagOverlay==1 & isfield(im,'reg_method')
im.pos_moveable=im.pos;
im.h_moveable=im.h;
im.hm = im.ht;
end

updateImages( handles );

function editWid_Callback(hObject, eventdata, handles)
global im

foo = str2num( get(handles.editWid,'string') );
if isempty(foo)
     set(handles.editWid,'string',num2str(im.h));
else if foo==0
     set(handles.editWid,'string',num2str(im.h));
     else
     im.h = foo;
     end
end
    
if im.flagOverlay==1 & isfield(im,'reg_method')
im.pos_moveable=im.pos;
im.h_moveable=im.h;
im.hm = im.ht;
end
    
updateImages( handles );

function pushbuttonUp_moveable_Callback(hObject, eventdata, handles)
global im
fv = im.flagView; 

im.pos_moveable(fv)=max(im.pos_moveable(fv)-1,1);

if im.flagOverlay==1 & isfield(im,'reg_method')
im.pos=im.pos_moveable;
im.h=im.h_moveable;
im.ht = im.hm;
end

updateImages( handles );

function pushbuttonUp10_moveable_Callback(hObject, eventdata, handles)
global im
fv = im.flagView; 

im.pos_moveable(fv)=max(im.pos_moveable(fv)-10,1);

if im.flagOverlay==1 & isfield(im,'reg_method')
im.pos = im.pos_moveable;
im.h = im.h_moveable;
im.ht = im.hm;
end
updateImages( handles );

function pushbuttonDown_moveable_Callback(hObject, eventdata, handles)
global im
fv = im.flagView; 

im.pos_moveable(fv)=min(im.pos_moveable(fv)+1,im.posMax_moveable(fv));

if im.flagOverlay==1 & isfield(im,'reg_method')
im.pos = im.pos_moveable;
im.h = im.h_moveable;
im.ht = im.hm;
end

updateImages( handles );

function pushbuttonDown10_moveable_Callback(hObject, eventdata, handles)
global im
fv = im.flagView; 

im.pos_moveable(fv)=min(im.pos_moveable(fv)+10,im.posMax_moveable(fv));

if im.flagOverlay==1 & isfield(im,'reg_method')
im.pos=im.pos_moveable;
im.h=im.h_moveable;
im.ht = im.hm;
end

updateImages( handles );


function editWid_moveable_Callback(hObject, eventdata, handles)
global im
foo = str2num( get(handles.editWid_moveable,'string') );
if isempty(foo)
    set(handles.editWid_moveable,'string',num2str(im.h_moveable));
else if foo==0
     set(handles.editWid_moveable,'string',num2str(im.h_moveable));
     else
     im.h_moveable = foo;
     end
end
    
if im.flagOverlay==1 & isfield(im,'reg_method')
im.pos=im.pos_moveable;
im.h=im.h_moveable;
im.ht = im.hm;
end
    
updateImages(handles);

function editVoxelSizeTarget_Callback(hObject, eventdata, handles)
global im

if im.flagOverlay==1
    foo = str2num( get(handles.editVoxelSizeTarget,'string') );
    if isempty(foo)
         set(handles.editWid,'string',num2str(im.ht));
    else if length(foo)~=3 || min(foo)<=0
         set(handles.editWid,'string',num2str(im.ht));
         else
         im.ht = foo;
         end
    end
    
    [nyt nxt nzt] = size(im.It);
    im.posMax = [nxt nyt nzt];
    im.posMax_moveable = im.posMax;
    im.hm= im.ht;    
else
    foo = str2num( get(handles.editVoxelSizeTarget,'string') );
    if isempty(foo)
         set(handles.editWid,'string',num2str(im.h));
    else if length(foo)~=3 || min(foo)<=0
         set(handles.editWid,'string',num2str(im.h));
         else
         im.ht = foo;
         end
    end
    [nyt nxt nzt] = size(im.It);
    im.posMax = [nxt nyt nzt];   
end
updateImages(handles);

function editVoxelSizeMoveable_Callback(hObject, eventdata, handles)
global im

if  im.flagOverlay == 1
    im.posMax_moveable = im.posMax;
    im.hm = im.ht;     
else
    foo = str2num( get(handles.editVoxelSizeMoveable,'string') );
    if   isempty(foo)
         set(handles.editWid_moveable,'string',num2str(im.h_moveable));
    else if length(foo)~=3 || min(foo)<=0
         set(handles.editWid_moveable,'string',num2str(im.h_moveable));
         else
         im.hm = foo;
        end
    end     
    [nym nxm nzm] = size(im.Im);
    im.posMax_moveable = [nxm nym nzm];
end
updateImages( handles );

% --- Executes when selected object is changed in uipanelProjection1.
function uipanelProjection1_SelectionChangeFcn(hObject, eventdata, handles)
updateImages( handles )


% --- Executes when selected object is changed in uipanelProjection2.
function uipanelProjection2_SelectionChangeFcn(hObject, eventdata, handles)
updateImages( handles )


% --- Executes on mouse press over axes background.
function axes1_ButtonDownFcn(~, eventdata, handles)

global im

pos = get(gca,'CurrentPoint');
pos = pos(1,:);
   
fv = im.flagView;
if fv==1      %%% xz
    fv0 = 1;
    fv1 = 2;
    fv2 = 3;
    pos = [im.pos(1) pos(1)/im.ht(fv1)  pos(2)/im.ht(fv2)];  %%% Y X Z
elseif fv==2  %%% yz
    fv0 = 2;
    fv1 = 1;
    fv2 = 3;
    pos = [pos(1)/im.ht(fv1)   im.pos(2)  pos(2)/im.ht(fv2)];  %%% Y X Z
elseif fv==3  %%% xy
    fv0 = 3;
    fv1 = 2;
    fv2 = 1;
    pos = [pos(2)/im.ht(fv2)   pos(1)/im.ht(fv1)   im.pos(3)];  %%% Y X Z
end

% % im.flagButton = 1; % Marker add
% % im.flagButton = 2; % move
% % im.flagButton = 3; % Zoom
% % im.flagButton = 4; % point
% % im.flagButton = 5; % ROI add

if im.flagButton==2  % move to Marker
    if im.MarkerState==0 | im.MarkerState==2
        pos1 = im.MarkerPos1;
        r = sum( (ones(size(pos1,1),1)*pos-pos1).^2,2).^0.5;
        [foo,idx] = min(r);
        pos = pos1(idx,:);
        im.MarkerPoint = idx;  % The sequence number of Marker to be moved to
        im.MarkerState = 1;
        set(handles.textMarkerStatus,'string', sprintf('Move Marker %d in Target Axes',im.MarkerPoint))
        set(handles.pushbuttonMarkerCancel,'visible','on')
        set(handles.pushbuttonMarkerDelete,'visible','on')
        
    elseif im.MarkerState==1
        im.MarkerPos1(im.MarkerPoint,:) = pos;
        im.MarkerState = 0;
        im.MarkerPoint = 0;
        set(handles.textMarkerStatus,'string','');
        set(handles.pushbuttonSave,'visible','on')
        set(handles.pushbuttonMarkerCancel,'visible','off')
        set(handles.pushbuttonMarkerDelete,'visible','off')
    end

elseif im.MarkerState==0 & im.flagButton~=4 & im.flagButton~=5,
        % Add marker and wait for corresponding marker in axes 2
        im.MarkerState = 1;
        im.MarkerPos1(end+1,1:3) = pos;
        set(handles.textMarkerStatus,'string', sprintf('PlaceS Marker %d in Moveable Axes',im.MarkerNum+1))
        set(handles.pushbuttonMarkerDelete,'visible','on')
    
elseif im.MarkerState==2 & im.flagButton~=4 & im.flagButton~=5,
        % Marker pair completed
        im.MarkerState = 0;
        im.MarkerNum = im.MarkerNum + 1;
        im.MarkerPos1(end+1,1:3) = pos;
        set(handles.textMarkerStatus,'string','');
        set(handles.pushbuttonMarkerDelete,'visible','off');
        set(handles.pushbuttonSave,'visible','on');
        
elseif im.flagButton == 5 & fv == 3       % ROI add
    if  im.ROI_MarkerState == 0, % no active ROI
        im.ROI_MarkerPos1(end+1,1:3) = pos;  % this is not pixel coordinate, but taking into account im.ht voxel sizes
        im.ROI_MarkerNum = im.ROI_MarkerNum + 1;
        im.ROI_MarkerState = 1;
        set(handles.textMarkerStatus,'string', sprintf('Place ROI Marker %d in Moveable (Flow) Axes',im.ROI_MarkerNum));
    elseif im.ROI_MarkerState == 2, % ROI added in FLOW but not in FITC axes
        im.ROI_MarkerPos1(end+1,1:3) = pos;
        im.ROI_MarkerState = 0;
    elseif im.ROI_MarkerState == 1,
        set(handles.textMarkerStatus,'string', sprintf('Error: ROI %d already selected in FITC but not in FLOW Axes',im.ROI_MarkerNum));
    end
end

im.pos = pos;          %%%   Y X Z
updateImages( handles )



% --- Executes on mouse press over axes background.
function axes2_ButtonDownFcn(hObject, eventdata, handles)

global im

pos = get(gca,'CurrentPoint');
pos = pos(1,:);
posC = im.pos_moveable;

fv = im.flagView;

if fv==1  %%%  xz
    fv0 = 1;
    fv1 = 2;
    fv2 = 3;
    pos = [im.T*[posC(1) pos(1)/im.hm(fv1)   pos(2)/im.hm(fv2)  1]']';
elseif fv==2  %%% yz
    fv0 = 2;
    fv1 = 1;
    fv2 = 3;
    pos = [im.T*[pos(1)/im.hm(fv1) posC(2) pos(2)/im.hm(fv2) 1]']';
elseif fv==3  %%% xy
    fv0 = 3;
    fv1 = 2;
    fv2 = 1;
    pos = [im.T*[pos(2)/im.hm(fv2) pos(1)/im.hm(fv1) posC(3) 1]']';
end

pos=pos(1:3);

% % im.flagButton = 1; % Marker add
% % im.flagButton = 2; % move
% % im.flagButton = 3; % Zoom
% % im.flagButton = 4; % point
% % im.flagButton = 5; % ROI add

if im.flagButton==2 % move
    if im.MarkerState==0 | im.MarkerState==1
        pos2 = im.MarkerPos2;
        r = sum( (ones(size(pos2,1),1)*pos-pos2).^2,2).^0.5;
        [foo,idx] = min(r);
        pos = pos2(idx,:);
        im.MarkerPoint = -idx;
        im.MarkerState = 2;
        set(handles.textMarkerStatus,'string', sprintf('Move Marker %d in Moveable Axes',im.MarkerPoint))
        set(handles.pushbuttonMarkerCancel,'visible','on')
        set(handles.pushbuttonMarkerDelete,'visible','on')
        pos = [im.T\[pos 1]']';
        pos = pos(1:3);
        
    elseif im.MarkerState==2
        im.MarkerPos2(-im.MarkerPoint,:) = pos;
        im.MarkerState = 0;
        im.MarkerPoint = 0;
        set(handles.textMarkerStatus,'string','');
        set(handles.pushbuttonSave,'visible','on')
        set(handles.pushbuttonMarkerCancel,'visible','off')
        set(handles.pushbuttonMarkerDelete,'visible','off')
        pos = [im.T\[pos 1]']';
        pos = pos(1:3);
    end
    
elseif im.MarkerState==0 & im.flagButton~=4 & im.flagButton~=5,
      % Add marker and wait for corresponding marker in axes 2
        im.MarkerState = 2;
        im.MarkerPos2(end+1,1:3) = pos;
        set(handles.textMarkerStatus,'string', sprintf('Place Marker %d in Target Axes',im.MarkerNum+1))
        set(handles.pushbuttonMarkerDelete,'visible','on')
        pos = [im.T\[pos 1]']';
        pos = pos(1:3);

elseif im.MarkerState==1 & im.flagButton~=4 & im.flagButton~=5,
    % Marker pair completed
       im.MarkerState = 0;
       im.MarkerNum = im.MarkerNum + 1;
       im.MarkerPos2(end+1,1:3) = pos;

       set(handles.textMarkerStatus,'string','');
       set(handles.pushbuttonMarkerDelete,'visible','off')
       set(handles.pushbuttonSave,'visible','on')
       pos = [im.T\[pos 1]']';
       pos = pos(1:3);
    
elseif  im.flagButton == 5 & fv == 3,  

    %   select a circle...
    xc = pos(1); % this is not pixel coordinate, but taking into account im.hm voxel sizes
    yc = pos(2);
    hhh = msgbox('Select one point on circle circumference');
    uiwait(hhh);
    [xr,yr] = ginput(1); % selected point on circumference, % this is not pixel coordinate, but taking into account im.hm voxel sizes
    r = sqrt((xc-xr)^2 + (yc-yr)^2); % pdist([xc xr],[yc yr]); % radius of circle (taking into account im.hm voxel sizes
        
    % create BW masks and circle lines in transformed and original coordinates...
    dangle = 1/(4*r);
    angle = 0:dangle:2*pi;
    
    %xCircl0 = xc0+r.*cos(angle);
    %yCircl0 = yc0+r.*sin(angle);
    %BW0 = roipoly(im.enface1,xCircl0,yCircl0); % I will not remember masks. They can always be generated with this command
    
    xCircl = xc+r.*cos(angle); % coordinates taking into account im.hm voxel sizes
    yCircl = yc+r.*sin(angle);
    xCircl = xCircl'; yCircl = yCircl';
    CirclCircumf = [xCircl yCircl];
    %BW = roipoly(im.newEnface,xCircl,yCircl);
    fooFlag =1;
    if im.ROI_MarkerState == 0, % no active marker
        im.ROI_MarkerNum = im.ROI_MarkerNum + 1;
        im.ROI_MarkerState = 2;
        set(handles.textMarkerStatus,'string', sprintf('Place ROI Marker %d in Target (FITC) Axes',im.ROI_MarkerNum));
    elseif im.ROI_MarkerState == 1, % marker placed on axes 1 and waiting for axes 2 here
        im.ROI_MarkerState = 0;
        set(handles.textMarkerStatus,'string','');
    elseif im.ROI_MarkerState == 2, % marker placed on axes 2 already => ERROR 
        set(handles.textMarkerStatus,'string','ERROR: Marker placed alredy on Axes 2, but not 1');
        fooFlag =0 ;
    end;
    if fooFlag,
        im.ROI_MarkerPos2.circCenter(im.ROI_MarkerNum,1:3) = pos;
        im.ROI_MarkerPos2.circR(im.ROI_MarkerNum) = r;
        im.ROI_MarkerPos2.circCircumf{im.ROI_MarkerNum} = CirclCircumf;
    end;
        
elseif im.flagButton ~= 5,
    pos = [im.T\[pos 1]']';
    pos = pos(1:3);
    
end

im.pos_moveable = pos(1:3);
updateImages( handles )


% --- Executes on button press in pushbuttonMarkerDelete.
function pushbuttonMarkerDelete_Callback(hObject, eventdata, handles)
global im

if im.MarkerPoint == 0
    if im.MarkerState==1
        im.MarkerPos1(end,:) = [];
    elseif im.MarkerState==2
        im.MarkerPos2(end,:) = [];
    end
else
    im.MarkerPos1(abs(im.MarkerPoint),:) = [];
    im.MarkerPos2(abs(im.MarkerPoint),:) = [];
    im.MarkerPoint = 0;
    set(handles.pushbuttonSave,'visible','on')
end


im.MarkerState=0;
set(handles.textMarkerStatus,'string','');
set(handles.pushbuttonMarkerDelete,'visible','off')
set(handles.pushbuttonMarkerCancel,'visible','off')
updateImages( handles )







% --- Executes on button press in pushbuttonAxes2Linc.
function pushbuttonAxes2Linc_Callback(hObject, eventdata, handles)
global im
im.CrangeM(1) = min(im.CrangeM(1)*2,im.CrangeM(2));
if im.CrangeM(1)==0
    im.CrangeM(1) = 1;
end
set(handles.editAxes2Crange,'string',sprintf('%d ',im.CrangeM) );
updateImages( handles )

% --- Executes on button press in pushbuttonAxes2Ldec.
function pushbuttonAxes2Ldec_Callback(hObject, eventdata, handles)
global im
im.CrangeM(1) = im.CrangeM(1)/2; %max(im.CrangeM(1)/2,0);
% if im.CrangeM(1)<=1
%     im.CrangeM(1) = 0;
% end
set(handles.editAxes2Crange,'string',sprintf('%d ',im.CrangeM) );
updateImages( handles )


% --- Executes on button press in pushbuttonAxes2Udec.
function pushbuttonAxes2Udec_Callback(hObject, eventdata, handles)
global im
im.CrangeM(2) = max(im.CrangeM(2)/2,im.CrangeM(1));
set(handles.editAxes2Crange,'string',sprintf('%d ',im.CrangeM) );
updateImages( handles )


% --- Executes on button press in pushbuttonAxes2Uinc.
function pushbuttonAxes2Uinc_Callback(hObject, eventdata, handles)
global im
im.CrangeM(2) = im.CrangeM(2)*2;
set(handles.editAxes2Crange,'string',sprintf('%d ',im.CrangeM) );
updateImages( handles )


function editAxes2Crange_Callback(hObject, eventdata, handles)
global im

crange = str2num(get(handles.editAxes2Crange,'string'));
if length(crange)~=2
    set(handles.editAxes2Crange,'string',sprintf('%d ',im.CrangeM) );
elseif crange(1)>crange(2)
    set(handles.editAxes2Crange,'string',sprintf('%d ',im.CrangeM) );
else
    im.CrangeM = crange;
    updateImages( handles );
end



function editAxes1Crange_Callback(hObject, eventdata, handles)
global im

crange = str2num(get(handles.editAxes1Crange,'string'));
if length(crange)~=2
    set(handles.editAxes1Crange,'string',sprintf('%d ',im.CrangeT) );
elseif crange(1)<0 | crange(1)>crange(2)
    set(handles.editAxes1Crange,'string',sprintf('%d ',im.CrangeT) );
else
    im.CrangeT = crange;
    updateImage( handles );
end


% --- Executes on button press in pushbuttonAxes1Ldec.
function pushbuttonAxes1Ldec_Callback(hObject, eventdata, handles)
global im
if im.CrangeT(1)<=1
    im.CrangeT(1) = 0;
end
im.CrangeT(1) = max(im.CrangeT(1)/2,0);
set(handles.editAxes1Crange,'string',sprintf('%d ',im.CrangeT) );
updateImages( handles )


% --- Executes on button press in pushbuttonAxes1Uinc.
function pushbuttonAxes1Uinc_Callback(hObject, eventdata, handles)
global im
im.CrangeT(2) = im.CrangeT(2)*2;
set(handles.editAxes1Crange,'string',sprintf('%d ',im.CrangeT) );
updateImages( handles )



% --- Executes on button press in pushbuttonAxes1Udec.
function pushbuttonAxes1Udec_Callback(hObject, eventdata, handles)
global im
im.CrangeT(2) = max(im.CrangeT(2)/2,im.CrangeT(1));
set(handles.editAxes1Crange,'string',sprintf('%d ',im.CrangeT) );
updateImages( handles )



% --- Executes on button press in pushbuttonAxes1Linc.
function pushbuttonAxes1Linc_Callback(hObject, eventdata, handles)
global im
im.CrangeT(1) = min(im.CrangeT(1)*2,im.CrangeT(2));
if im.CrangeT(1)==0
    im.CrangeT(1) = 1;
end
set(handles.editAxes1Crange,'string',sprintf('%d ',im.CrangeT) );
updateImages( handles )




% --- Executes on button press in pushbuttonRegisterXY.
function pushbuttonRegisterXY_Callback(hObject, eventdata, handles)

global im
im.reg_method=1;
pT = im.MarkerPos1;
pM = im.MarkerPos2;


% Get XY transformation
A = [pT(:,1:2) ones(size(pT,1),1)];
Ttmp = [A\pM(:,1) A\pM(:,2) [0 0 1]']';
T = eye(4,4);
T(1:2,1:2) = Ttmp(1:2,1:2);
T(1:2,4) = Ttmp(1:2,3);

im.Im_transform = vol2vol(im.Im,im.It,T);
im.T = T;
im.coreg_flag=1;
updateImages( handles )
set(handles.pushbuttonSave,'visible','on');


% --- Executes on button press in pushbuttonRegisterXYZ.
function pushbuttonRegisterXYZ_Callback(hObject, eventdata, handles)
global im
im.reg_method=1;

pT = im.MarkerPos1;
pM = im.MarkerPos2;


% Get XYZ transformation
A = [pT(:,1:3) ones(size(pT,1),1)];
Ttmp = [A\pM(:,1) A\pM(:,2) A\pM(:,3)]';
T = eye(4,4);
T(1:3,1:4) = Ttmp(1:3,1:4);

% transform the moveable image
im.Im_transform = vol2vol(im.Im,im.It,T);

im.T = T;
im.coreg_flag=1;
updateImages( handles )
set(handles.pushbuttonSave,'visible','on');

%%%%%%%%%%%%%%%%%%%%%%%
% UPDATE IMAGES
%%%%%%%%%%%%%%%%%%%%%%%

function updateImages( handles )

global im


fv = im.flagView;
if fv==1
    fv0 = 1;
    fv1 = 2;
    fv2 = 3;
    set(handles.textRange,'string',sprintf('Y [%.1f %.1f]',(im.pos(fv0)-im.h/2)*im.ht(fv0),(im.pos(fv0)+im.h/2)*im.ht(fv0)));
    set(handles.textRange_moveable,'string',sprintf('Y [%.1f %.1f]',(im.pos_moveable(fv0)-im.h_moveable/2)*im.hm(fv0),(im.pos_moveable(fv0)+im.h_moveable/2)*im.hm(fv0)));

elseif fv==2
    fv0 = 2;
    fv1 = 1;
    fv2 = 3;
    set(handles.textRange,'string',sprintf('X [%.1f %.1f]',(im.pos(fv0)-im.h/2)*im.ht(fv0),(im.pos(fv0)+im.h/2)*im.ht(fv0)));
    set(handles.textRange_moveable,'string',sprintf('X [%.1f %.1f]',(im.pos_moveable(fv0)-im.h_moveable/2)*im.hm(fv0),(im.pos_moveable(fv0)+im.h_moveable/2)*im.hm(fv0)));

elseif fv==3
    fv0 = 3;
    fv1 = 2;
    fv2 = 1;
    set(handles.textRange,'string',sprintf('Z [%.1f %.1f]',(im.pos(fv0)-im.h/2)*im.ht(fv0),(im.pos(fv0)+im.h/2)*im.ht(fv0)));
    set(handles.textRange_moveable,'string',sprintf('Z [%.1f %.1f]',(im.pos_moveable(fv0)-im.h_moveable/2)*im.hm(fv0),(im.pos_moveable(fv0)+im.h_moveable/2)*im.hm(fv0)));

end

set(handles.editWid,'string',num2str(im.h));
set(handles.editWid_moveable,'string',num2str(im.h_moveable));

% get overlay options
im.flagOverlay = get(handles.checkboxOverlay,'value');
overlayRight = 0;
if strcmpi(get(handles.togglebuttonLeftRight,'string'),'Right')
    overlayRight = 1;
end
overlayInvertRight = get(handles.checkboxInvertRight,'value');
overlayInvertLeft = get(handles.checkboxInvertLeft,'value');
overlayChColor = get(handles.popupmenuChColor,'value');

% Axes 1
axes(handles.axes1)
pos = [1:size(im.It,fv)];
lst = find(pos>=(im.pos(fv0)-im.h/2) & pos<=(im.pos(fv0)+im.h/2));
if get(handles.radiobuttonMin1,'value')
    if fv==1
        foo = squeeze(min(im.It(lst,:,:),[],1))';
    elseif fv==2
        foo = squeeze(min(im.It(:,lst,:),[],2))';
    elseif fv==3
        foo = min(im.It(:,:,lst),[],3);
    end
else
    if fv==1
        foo = squeeze(max(im.It(lst,:,:),[],1))';
    elseif fv==2
        foo = squeeze(max(im.It(:,lst,:),[],2))';
    elseif fv==3
        foo = max(im.It(:,:,lst),[],3);
    end
end
imagesc([1:size(foo,2)]*im.ht(fv1),[1:size(foo,1)]*im.ht(fv2),foo, im.CrangeT);

axis fill
colormap gray

hold on
for ii=1:size(im.MarkerPos1,1)
    if im.MarkerPos1(ii,fv0)>=(im.pos(fv0)-im.h/2) & im.MarkerPos1(ii,fv0)<=(im.pos(fv0)+im.h/2)
        ht = text((im.MarkerPos1(ii,fv1)+8)*im.ht(fv1),im.MarkerPos1(ii,fv2)*im.ht(fv2),sprintf('%d',ii),'fontSize',12);
        plot(im.MarkerPos1(ii,fv1)*im.ht(fv1),im.MarkerPos1(ii,fv2)*im.ht(fv2),'o','MarkerSize',5,'MarkerEdgeColor','r','MarkerFaceColor','r');
        set(ht,'color','r');
        if im.MarkerPoint==ii
            set(ht,'color','g');
        end
    end
end

if size(im.MarkerPos1,1)<size(im.MarkerPos2,1)
    pT = [im.T\[im.MarkerPos2(end,:) 1]']';
    ht = text(pT(fv1)*im.ht(fv1), pT(fv2)*im.ht(fv2),sprintf('%d',size(im.MarkerPos2,1)),'fontsize',12);
    set(ht,'color',[1 1 0]);
end

% % im.flagButton = 1; % Marker add
% % im.flagButton = 2; % move
% % im.flagButton = 3; % Zoom
% % im.flagButton = 4; % point
% % im.flagButton = 5; % ROI add

if im.flagButton==4 | im.flagButton==3
    pT = im.pos;
    ht = text(pT(fv1)*im.ht(fv1),pT(fv2)*im.ht(fv2),'X');
    set(ht,'color',[1 0 1]);
end
hold off

freezeColors;



% Axes 2
axes(handles.axes2);
    
if im.flow_flag==1
    
    image_src_axes2=im.Im; 
   
elseif im.coreg_flag==0
        
    image_src_axes2=im.Im; 
    
elseif im.coreg_flag==1
    
    image_src_axes2=im.Im_transform;

end

pos = [1:size(image_src_axes2,fv)];

lst = find(pos>=(im.pos_moveable(fv0)-im.h_moveable/2) & pos<=(im.pos_moveable(fv0)+im.h_moveable/2));

if im.flow_flag==1
       
    temp_stack=im.Im(:,:,lst);
    foo2=max(abs(temp_stack),[],3);
    max_c=max(foo2(:));
    [nx ny]=size(foo2);
    for k=1:nx
        for l=1:ny
            
        sign_stack(k,l)=sign(temp_stack(k,l,find(abs(temp_stack(k,l,:))==foo2(k,l),1)));
   
        end
    end
    
    foo2=foo2.*sign_stack;
    
else
    
if get(handles.radiobuttonMin2,'value')
    if fv==1
        foo2 = squeeze(min(image_src_axes2(lst,:,:),[],1))';
    elseif fv==2
        foo2 = squeeze(min(image_src_axes2(:,lst,:),[],2))';
    elseif fv==3
        foo2 = min(image_src_axes2(:,:,lst),[],3);
    end
elseif get(handles.radiobuttonMax2,'value')
    if fv==1
        foo2 = squeeze(max(image_src_axes2(lst,:,:),[],1))';
    elseif fv==2
        foo2 = squeeze(max(image_src_axes2(:,lst,:),[],2))';
    elseif fv==3
        foo2 = max(image_src_axes2(:,:,lst),[],3);
    end
else
    if fv==1
        foo2 = squeeze(mean(image_src_axes2(lst,:,:),1))';
    elseif fv==2
        foo2 = squeeze(mean(image_src_axes2(:,lst,:),2))';
    elseif fv==3
        foo2 = mean(image_src_axes2(:,:,lst),3);
    end
end


end


if (im.flagOverlay==0 | overlayRight==0) & ~isfield(im,'Im_transform')    %%% Display non-transformed OCT image

  if im.flow_flag==0
    imagesc([1:size(foo2,2)]*im.hm(fv1),[1:size(foo2,1)]*im.hm(fv2),foo2,im.CrangeM);
    axis image
    colormap gray

    if ~isempty(im.MarkerPos2)
        
        pM = [im.T\[im.MarkerPos2 ones(size(im.MarkerPos2,1),1)]']';
        for ii=1:size(im.MarkerPos2,1)
            if pM(ii,fv0)>=(im.pos_moveable(fv0)-im.h_moveable/2) & pM(ii,fv0)<=(im.pos_moveable(fv0)+im.h_moveable/2)
                
                hold on
                ht = text((pM(ii,fv1)+8)*im.hm(fv1),pM(ii,fv2)*im.hm(fv2),sprintf('%d',ii),'fontsize',12);
                set(ht,'color','m');
                plot(pM(ii,fv1)*im.hm(fv1),pM(ii,fv2)*im.hm(fv2),'o','MarkerSize',5,'MarkerEdgeColor','m','MarkerFaceColor','m');
                                           
                if im.MarkerPoint==-ii
                    set(ht,'color','g');
                end                            
                
            end
        end
        
        if size(im.MarkerPos2,1)<size(im.MarkerPos1,1)
            pT = im.MarkerPos1(end,:);
            ht = text(pT(fv1)*im.hm(fv1),pT(fv2)*im.hm(fv2),sprintf('%d',size(im.MarkerPos1,1)),'fontsize',12);
            set(ht,'color',[0 1 1]);
        end
        if im.flagButton==4 | im.flagButton==3
            pM = im.pos;
            ht = text(pM(fv1)*im.hm(fv1),pM(fv2)*im.hm(fv2),'X');
            set(ht,'color',[1 0 1]);
        end
        hold off
    end
    
  else
    imagesc([1:size(foo2,2)]*im.hm(fv1),[1:size(foo2,1)]*im.hm(fv2),foo2,[-max_c max_c]);
    axis image
    cm1a = [32:-1:1]'*[0 1 0]/32; 
    cm1a(:,3) = 1;
    cm1b = [32:-1:1]'*[0 0 1]/32;
    cm2 = hot(64);
    colormap([cm1a; cm1b; cm2]);
    colorbar
    
    if ~isempty(im.MarkerPos2)
        
        pM = [im.T\[im.MarkerPos2 ones(size(im.MarkerPos2,1),1)]']';
        for ii=1:size(im.MarkerPos2,1)
            if pM(ii,fv0)>=(im.pos_moveable(fv0)-im.h_moveable/2) & pM(ii,fv0)<=(im.pos_moveable(fv0)+im.h_moveable/2)
                
                hold on
                ht = text((pM(ii,fv1)+8)*im.hm(fv1),pM(ii,fv2)*im.hm(fv2),sprintf('%d',ii),'fontsize',12);
                set(ht,'color','m');
                plot(pM(ii,fv1)*im.hm(fv1),pM(ii,fv2)*im.hm(fv2),'o','MarkerSize',5,'MarkerEdgeColor','m','MarkerFaceColor','m');
                                           
                if im.MarkerPoint==-ii
                    set(ht,'color','g');
                end                            
                
            end
        end
        if size(im.MarkerPos2,1)<size(im.MarkerPos1,1)
            pT = im.MarkerPos1(end,:);
            ht = text(pT(fv1)*im.hm(fv1),pT(fv2)*im.hm(fv2),sprintf('%d',size(im.MarkerPos1,1)),'fontsize',12);
            set(ht,'color',[0 1 1]);
        end
        if im.flagButton==4 | im.flagButton==3
            pM = im.pos;
            ht = text(pM(fv1)*im.hm(fv1),pM(fv2)*im.hm(fv2),'X');
            set(ht,'color',[1 0 1]);
        end
        hold off
      
      end
   end
end

if (im.flagOverlay==0 | overlayRight==0) & isfield(im,'Im_transform')
    
    imagesc([1:size(foo2,2)]*im.ht(fv1),[1:size(foo2,1)]*im.ht(fv2),foo2,im.CrangeM);
    axis image
    colormap gray

        if ~isempty(im.MarkerPos2)

        pT = [im.T\[im.MarkerPos2 ones(size(im.MarkerPos2,1),1)]']';
        pT_0= im.MarkerPos1;
        for ii=1:size(im.MarkerPos2,1)
            if pT(ii,fv0)>=(im.pos(fv0)-im.h/2) & pT(ii,fv0)<=(im.pos(fv0)+im.h/2)
                
                hold on
                ht = text((pT(ii,fv1)+8)*im.ht(fv1),pT(ii,fv2)*im.ht(fv2),sprintf('%d',ii),'fontsize',12);
                set(ht,'color','b');
                plot(pT(ii,fv1)*im.ht(fv1),pT(ii,fv2)*im.ht(fv2),'o','MarkerSize',5,'MarkerEdgeColor','b','MarkerFaceColor','b');

                hold on
                ht_0 = text((pT_0(ii,fv1)+8)*im.ht(fv1),pT_0(ii,fv2)*im.ht(fv2),sprintf('%d',ii),'fontsize',12);
                set(ht_0,'color','r');
                plot(pT_0(ii,fv1)*im.ht(fv1),pT_0(ii,fv2)*im.ht(fv2),'o','MarkerSize',5,'MarkerEdgeColor','r','MarkerFaceColor','r');

                
                if im.MarkerPoint==-ii
                    set(ht,'color','g');
                end
            end 
        end
        if size(im.MarkerPos2,1)<size(im.MarkerPos1,1)
            pT = im.MarkerPos1(end,:);
            ht = text(pT(fv1)*im.ht(fv1),pT(fv2)*im.ht(fv2),sprintf('%d',size(im.MarkerPos1,1)),'fontsize',12);
            set(ht,'color',[0 1 1]);
        end
        if im.flagButton==4 | im.flagButton==3
            pM = im.pos;
            ht = text(pM(fv1)*im.ht(fv1),pM(fv2)*im.ht(fv2),'X');
            set(ht,'color',[1 0 1]);
        end
        hold off
    end
end

freezeColors;

% check for overlay
if im.flagOverlay==1 
      
    % dispaly method 2
    temp1=1-exp(-4*foo/max(foo(:)));
    green=cat(3,zeros(size(temp1)),ones(size(temp1)),zeros(size(temp1)));

         
    [ny,nx] = size(foo);
    dx = 1/double(nx-1); nxa = 0:dx:1;
    dy = 1/double(ny-1); nya = 0:dy:1;
    % interpolate foo2 (OCT image) on foo size (FITC image size)
    [ny2,nx2] = size(foo2);
    dx2 = 1/double(nx2-1); nx2a = 0:dx2:1;
    dy2 = 1/double(ny2-1); ny2a = 0:dy2:1;    
    foo2 = interp2(nx2a', ny2a, double(foo2), nxa', nya);
    
    if overlayRight==1
        axes(handles.axes2)
    else
        axes(handles.axes1)
    end
    
    if im.flow_flag==1
    %imagesc([1:size(foo2,2)]*im.hm(fv1),[1:size(foo2,1)]*im.hm(fv2),foo2,[-max_c max_c]);
    imagesc([1:size(foo2,2)],[1:size(foo2,1)],[-max_c max_c]);
    axis off
    axis fill  
    cm1a = [32:-1:1]'*[0 1 0]/32; 
    cm1a(:,3) = 1;
    cm1b = [32:-1:1]'*[0 0 1]/32;
    cm2 = hot(64);
    colormap([cm1a; cm1b; cm2]);
    colorbar
    else
    %imagesc([1:size(foo2,2)]*im.hm(fv1),[1:size(foo2,1)]*im.hm(fv2),foo2,im.CrangeM);
    imagesc([1:size(foo2,2)],[1:size(foo2,1)],foo2,im.CrangeM);
    axis off
    axis fill    
    colormap gray 
    end
    
    hold on
    h_overlay=imshow(green);
    axis off
    axis fill
    hold off
    set(h_overlay,'AlphaData',temp1);

end

 fhandle_mousemove_figure='coregisterOCTand2PM(''figure_WindowButtonMotionFcn'',gcbo,[],guidata(gcbo))';
 set(handles.figure1,'WindowButtonMotionFcn',fhandle_mousemove_figure);

% set up button down callbacks
axes(handles.axes1)   
   
if  im.flagButton==1 | im.flagButton==2 | im.flagButton==4 | im.flagButton ==5
    fhandle1='coregisterOCTand2PM(''axes1_ButtonDownFcn'',gcbo,[],guidata(gcbo))';
    set(handles.axes1,'ButtonDownFcn',fhandle1);
    set(get(handles.axes1,'children'), 'ButtonDownFcn',fhandle1);

end

axes(handles.axes2)

if  im.flagButton==1 | im.flagButton==2 | im.flagButton==4 | im.flagButton ==5
    zoom off
    fhandle2='coregisterOCTand2PM(''axes2_ButtonDownFcn'',gcbo,[],guidata(gcbo))';
    set(handles.axes2,'ButtonDownFcn',fhandle2);
    set(get(handles.axes2,'children'), 'ButtonDownFcn',fhandle2);
elseif im.flagButton == 3
    h = zoom;
    set(h,'enable','on');
    linkaxes([handles.axes1 handles.axes2],'xy')
end

function axes1_WindowButtonMotionFcn(hObject, eventdata, handles)
mouse_pos_axes1 = get(gca,'currentpoint');
xlim = get(gca,'xlim');
ylim = get(gca,'ylim');
outX = ~any(diff([xlim(1) mouse_pos_axes1(1,1) xlim(2)])<0);
outY = ~any(diff([ylim(1) mouse_pos_axes1(1,2) ylim(2)])<0);
if outX&outY
    set(handles.coordinate_TPM,'string', sprintf('[%.1f %.1f]',mouse_pos_axes1(1,1), mouse_pos_axes1(1,2)));
else
    set(handles.coordinate_TPM,'string',['']);
end

function figure_WindowButtonMotionFcn(hObject, eventdata, handles)

mouse_pos= get(gcf,'currentpoint');

panel_pos_axes1 = handles.axes1.Position;
panel_pos_axes2 = handles.axes2.Position;

if mouse_pos(1)>=panel_pos_axes1(1) & mouse_pos(1)<=(panel_pos_axes1(1)+panel_pos_axes1(3)) & mouse_pos(2)>=panel_pos_axes1(2) & mouse_pos(2)<=(panel_pos_axes1(2)+panel_pos_axes1(4)) 

    axes(handles.axes1)
    mouse_pos_ax=get(gca,'currentpoint');
    set(handles.coordinate_TPM,'string', sprintf('[%.1f %.1f]',mouse_pos_ax(1,1), mouse_pos_ax(1,2)),'fontSize',12);
    set(handles.coordinate_OCT,'string', ['']);

else if mouse_pos(1)>=panel_pos_axes2(1) & mouse_pos(1)<=(panel_pos_axes2(1)+panel_pos_axes2(3)) & mouse_pos(2)>=panel_pos_axes2(2) & mouse_pos(2)<=(panel_pos_axes2(2)+panel_pos_axes2(4)) 
 
    axes(handles.axes2)
    mouse_pos_ax=get(gca,'currentpoint');
    set(handles.coordinate_OCT,'string', sprintf('[%.1f %.1f]',mouse_pos_ax(1,1), mouse_pos_ax(1,2)),'fontSize',12);
    set(handles.coordinate_TPM,'string', ['']);

    else 
        
            set(handles.coordinate_TPM,'string', ['']);
            set(handles.coordinate_OCT,'string', ['']);
    end
    
end
    

% --- Executes on button press in pushbuttonSave.
function pushbuttonSave_Callback(hObject, eventdata, handles)
global im

[fn,pn] = uiputfile( '*.mat', 'Save file');
if fn~=0
    hwait=waitbar(0,'Saving');
    save([pn fn],'im','-v7.3');
    close(hwait)
    %set(handles.pushbuttonSave,'visible','off')
end


% --- Executes when selected object is changed in uipanelAxesButton.
function uipanelAxesButton_SelectionChangeFcn(hObject, eventdata, handles)
global im

if get(handles.radiobuttonMarker,'value')
    im.flagButton = 1; % Marker add
elseif get(handles.radiobuttonMarkerMove,'value')
    im.flagButton = 2; % move
elseif get(handles.radiobuttonZoom,'value')
    im.flagButton = 3; % Zoom
elseif get(handles.radiobuttonPoint,'value')
    im.flagButton = 4; % point
elseif get(handles.radiobutton_ROI_Add,'value')
    im.flagButton = 5; % ROI add
end

im.MarkerPoint = 0;
    
if im.MarkerState==1
    im.MarkerPos1(end,:) = [];
elseif im.MarkerState==2
    im.MarkerPos2(end,:) = [];
end

im.MarkerState=0;
set(handles.textMarkerStatus,'string','');
set(handles.pushbuttonMarkerDelete,'visible','off')
set(handles.pushbuttonMarkerCancel,'visible','off')
updateImages( handles )


% --- Executes when selected object is changed in uipanelView.
function uipanelView_SelectionChangeFcn(hObject, eventdata, handles)
global im

if get(handles.radiobuttonXY,'value')
    im.flagView=3;
elseif get(handles.radiobuttonXZ,'value')
    im.flagView=1;
elseif get(handles.radiobuttonYZ,'value')
    im.flagView=2;
end    
updateImages( handles )

function radiobuttonXY_Callback(hObject, eventdata, handles)
global im
if get(handles.radiobuttonXY,'value')
    im.flagView=3;
end
updateImages( handles )

function radiobuttonXZ_Callback(hObject, eventdata, handles)
global im
if get(handles.radiobuttonXY,'value')
    im.flagView=1;
end
updateImages( handles )

function radiobuttonYZ_Callback(hObject, eventdata, handles)
global im
if get(handles.radiobuttonYZ,'value')
    im.flagView=2;
end
updateImages( handles )

% --- Executes on button press in pushbuttonMarkerCancel.
function pushbuttonMarkerCancel_Callback(hObject, eventdata, handles)
global im

im.MarkerPoint = 0;
im.MarkerState=0;
set(handles.textMarkerStatus,'string','');
set(handles.pushbuttonMarkerCancel,'visible','off')
updateImages( handles )


% --- Executes on button press in checkboxOverlay.
function checkboxOverlay_Callback(hObject, eventdata, handles)
global im
if im.flagOverlay==1 & isfield(im,'reg_method')

im.pos_moveable=im.pos;
im.h_moveable=im.h;
im.hm= im.ht;

else
    
    
    
end

updateImages( handles )


% --- Executes on button press in togglebuttonLeftRight.
function togglebuttonLeftRight_Callback(hObject, eventdata, handles)
if strcmpi(get(hObject,'string'),'right')
    set(hObject,'string','Left');
else
    set(hObject,'string','Right');
end
    

% --- Executes on button press in checkboxInvertLeft.
function checkboxInvertLeft_Callback(hObject, eventdata, handles)
updateImages( handles )


% --- Executes on button press in checkboxInvertRight.
function checkboxInvertRight_Callback(hObject, eventdata, handles)
updateImages( handles )


% --- Executes on selection change in popupmenuChColor.
function popupmenuChColor_Callback(hObject, eventdata, handles)
updateImages( handles )





% --- Executes on button press in pushbutton_ExportFlow.
function pushbutton_ExportFlow_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_ExportFlow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global im;

% constants for getting scale coeff between phase and velocity
lambda = 1310e-6; %center wavelength in mm
n = 1.35; %index
exp_time = 21.377e-6; %in s
phase_to_vel_scale_factor = (lambda/(2*n))/(2*pi)/exp_time; % lambda / 2 because double path

% How to get velocities and flow
% doublemat = pi/(2^7-1)*double(handles.totaldopplermat(dz,dy,dx,:)); %127 phase steps convert to [rad]
% velocity_mmPerSec = phase_to_vel_scale_factor * unwrap(angle(exp(j*doublemat(i_z,i_y,i_x,i_t))));
% flow_L_per_min = velocity_mmPerSec * handles.pixel_area_mm2 / litersPer_mm3 * 60;

hwait = waitbar(0,'Processing ROIs...');


im.ROI_maxVel_mmPerSec=[];
im.ROI_velocityMap_mmPerSec=[];
im.ROI_averFlowLitPerMin=[];
    
for idx = 1:im.ROI_MarkerNum,
    waitbar(idx/im.ROI_MarkerNum,hwait,['Processing ROI number ' num2str(idx)]);
    xyzC_um = im.ROI_MarkerPos2.circCenter(idx,:);
    circleXY_um = im.ROI_MarkerPos2.circCircumf{idx};
    % find center location in pixels
    pos = polyval(im.Z2fit,[1:size(im.Im,3)]);
    lstZ = find(pos>=(xyzC_um(3)-2*im.hm(3)) & pos<=(xyzC_um(3)+im.hm(3)));
    BWcirc = double(roipoly( (1:size(im.Im,1)).*im.hm(1), (1:size(im.Im,2)).*im.hm(2), im.newEnface, circleXY_um(:,1), circleXY_um(:,2) ) );
    
    xmax = max(circleXY_um(:,1));
    xmin = min(circleXY_um(:,1));
    ymax = max(circleXY_um(:,2));
    ymin = min(circleXY_um(:,2));
    
    BWsq = double(roipoly( (1:size(im.Im,1)).*im.hm(1), (1:size(im.Im,2)).*im.hm(2), im.newEnface, [xmin xmax xmax xmin xmin], [ymin ymin ymax ymax ymin]) ); 
    
    %doublemat = pi/(2^7-1)*double(im.Im(:,:,lstZ));
    lst = find(BWcirc(:) > 0);
    %fooVel3D = doublemat.*0.0;
    
    velFrame = squeeze(mean(im.Im(:,:,lstZ),3));
    
%     velFrame = zeros(size(im.Im,1), size(im.Im,2) );
%     for jj = 1:length(lstZ),
%         fooPlane = squeeze(doublemat(:,:,jj));
%         for ii = 1:length(lst),
%             velFrame(lst(ii)) = velFrame(lst(ii)) + phase_to_vel_scale_factor * unwrap(angle( exp(j*fooPlane(lst(ii)))));
%         end;
%     end;
%     velFrame = velFrame ./ length(lstZ);  % this is now mean velocity (mm/sec) in 2D frame equal to newEnface image size


    averFlow_L_per_min = sum(sum(velFrame.*BWcirc)) * (im.hm(1)*im.hm(2)*1e-6) * 1e-6 * 60;
    %extract square ROI around vessel from velFrame...
%     xmin1 = round(xmin/im.hm(1));
%     xmax1 = round(xmax/im.hm(1));
%     ymin1 = round(ymin/im.hm(2));
%     ymax1 = round(ymax/im.hm(2));
%     extractVelSq = velFrame(ymin1:ymax1,xmin1:xmax1);
    % better extraction of ROI:
    BWfoo = bwboundaries(BWsq,8,'noholes');
    poo = BWfoo{1};
    xfoomin= min(poo(:,1));
    xfoomax= max(poo(:,1));
    yfoomin= min(poo(:,2));
    yfoomax= max(poo(:,2));
    extractVelSq = velFrame(xfoomin:xfoomax, yfoomin:yfoomax);
    
    im.ROI_averFlowLitPerMin(idx) = averFlow_L_per_min;
    im.ROI_velocityMap_mmPerSec{idx} = extractVelSq;
    
    velMax = max(extractVelSq(:));
    velMin = min(extractVelSq(:));
    
    if averFlow_L_per_min < 0, % if integral is negative
        im.ROI_maxVel_mmPerSec(idx) = velMin; % look for negative max Velocity
    else
        im.ROI_maxVel_mmPerSec(idx) = velMax; % else, look for positive maxVelocity
    end;
       
    
    
end;
close(hwait);
[fname, pname] = uiputfile('*.mat','Select File Name for imView3d Flow Export Data','_ExportedOCTflowFor_imView3d.mat');
    OCTflow.ROI_maxVel_mmPerSec = im.ROI_maxVel_mmPerSec;
    OCTflow.ROI_velocityMap_mmPerSec = im.ROI_velocityMap_mmPerSec;
    OCTflow.ROI_averFlowLitPerMin = im.ROI_averFlowLitPerMin;
    OCTflow.ROI_MarkerNum = im.ROI_MarkerNum;
    OCTflow.ROI_MarkerPosOCT = im.ROI_MarkerPos2;
    OCTflow.ROI_MarkerPosFITC = im.ROI_MarkerPos1;
    OCTflow.voxelSizes_FITC = im.ht;
    OCTflow.voxelSizes_OCT = im.hm;
    save([pname fname],'OCTflow','-mat');



% --- Executes on button press in pushbutton_DeleteROI.
function pushbutton_DeleteROI_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_DeleteROI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global im;

answer = inputdlg('ROI number to delete (cancel to exit)','Delete ROI',1,{'-1'});
if ~isempty(answer),
    foo = str2num(answer{1});
    if foo>0 & foo<=im.ROI_MarkerNum,
        im.ROI_MarkerNum = im.ROI_MarkerNum-1;
        im.ROI_MarkerPos1(foo,:)=[];
        im.ROI_MarkerPos2.circCenter(foo,:) = [];
        im.ROI_MarkerPos2.circCircumf(foo) = [];
        im.ROI_MarkerPos2.circR(foo) = [];
    end;
end;
