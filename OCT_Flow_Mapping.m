function varargout = OCT_Flow_Mapping(varargin)
% OCT_FLOW_MAPPING MATLAB code for OCT_Flow_Mapping.fig
%      OCT_FLOW_MAPPING, by itself, creates a new OCT_FLOW_MAPPING or raises the existing
%      singleton*.
%
%      H = OCT_FLOW_MAPPING returns the handle to a new OCT_FLOW_MAPPING or the handle to
%      the existing singleton*.
%
%      OCT_FLOW_MAPPING('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in OCT_FLOW_MAPPING.M with the given input arguments.
%
%      OCT_FLOW_MAPPING('Property','Value',...) creates a new OCT_FLOW_MAPPING or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before OCT_Flow_Mapping_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to OCT_Flow_Mapping_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help OCT_Flow_Mapping

% Last Modified by GUIDE v2.5 16-Jun-2020 10:39:14

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @OCT_Flow_Mapping_OpeningFcn, ...
                   'gui_OutputFcn',  @OCT_Flow_Mapping_OutputFcn, ...
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
end
% End initialization code - DO NOT EDIT

% --- Executes just before OCT_Flow_Mapping is made visible.
function OCT_Flow_Mapping_OpeningFcn(hObject, ~, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to OCT_Flow_Mapping (see VARARGIN)

% Choose default command line output for OCT_Flow_Mapping
handles.output = hObject;
% Update handles structure
guidata(hObject, handles);

global hndl
hndl = handles;

openingDlg(handles)

end


% UIWAIT makes OCT_Flow_Mapping wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = OCT_Flow_Mapping_OutputFcn(~, ~, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
end


%-- Update TP Angiogram Axis 
function update_TPMA_image(handles)
global im0 V hndl %#ok<*GVMIS> 

hndl = handles;
xLimit = xlim(handles.axes_TP_Ang);
yLimit = ylim(handles.axes_TP_Ang);

ll = round(get(handles.slider, 'Value'));

switch im0.gui.viewPlane
    case 'XY'
        if ll <= size(V.TP_XY, 3)
            imagesc(handles.axes_TP_Ang, V.TP_XY(:,:,ll), [im0.gui.TPthreshold*32 32])
        else
            im = zeros(size(V.TP_XY, 1), size(V.TP_XY, 2));
            imagesc(handles.axes_TP_Ang, im, [im0.gui.TPthreshold*32 32])
        end
    case 'YZ'
        if ll <= size(V.TP_YZ, 3)
            imagesc(handles.axes_TP_Ang, V.TP_YZ(:,:,ll), [im0.gui.TPthreshold*32 32])
        else
            im = zeros(size(V.TP_YZ, 1), size(V.TP_YZ, 2));
            imagesc(handles.axes_TP_Ang, im, [im0.gui.TPthreshold*32 32])
        end
    case 'XZ'
        if ll <= size(V.TP_XZ, 3)
            imagesc(handles.axes_TP_Ang, V.TP_XZ(:,:,ll), [im0.gui.TPthreshold*32 32])
        else
            im = zeros(size(V.TP_XZ, 1), size(V.TP_XZ, 2));
            imagesc(handles.axes_TP_Ang, im, [im0.gui.TPthreshold*32 32])
        end
end


xlim(handles.axes_TP_Ang, xLimit)
ylim(handles.axes_TP_Ang, yLimit)
colormap(handles.axes_TP_Ang, gray);
% colormap(handles.axes_TP_Ang, [0,0,0; 0,1,1; 1,1,0; 1,0,0]);
% caxis(handles.axes_TP_Ang, [0, 3])
handles.axes_TP_Ang.Color = 'k';
enableDefaultInteractivity(handles.axes_TP_Ang)
end


%-- Update OCT Flow Axis
function update_OCTV_image(handles)
global im0 V
xLimit = xlim(handles.axes_OCT_Flow);
yLimit = ylim(handles.axes_OCT_Flow);

ll = round(get(handles.slider, 'Value'));

switch im0.gui.viewPlane
    case 'XY'
        if ll <= size(V.OCT_Flow_XY, 3)
            imagesc(handles.axes_OCT_Flow, V.OCT_Flow_XY(:,:,ll), [im0.gui.colorbalLowerLim im0.gui.colorbarUpperLim])
        else
            im = zeros(size(V.OCT_Flow_XY, 1), size(V.OCT_Flow_XY, 2));
            imagesc(handles.axes_OCT_Flow, im, [im0.gui.colorbalLowerLim im0.gui.colorbarUpperLim])
        end
    case 'YZ'
        if ll <= size(V.OCT_Flow_YZ, 3)
            imagesc(handles.axes_OCT_Flow, V.OCT_Flow_YZ(:,:,ll), [im0.gui.colorbalLowerLim im0.gui.colorbarUpperLim])
        else
            im = zeros(size(V.OCT_Flow_YZ, 1), size(V.OCT_Flow_YZ, 2));
            imagesc(handles.axes_OCT_Flow, im, [im0.gui.colorbalLowerLim im0.gui.colorbarUpperLim])
        end
    case 'XZ'
        if ll <= size(V.OCT_Flow_XZ, 3)
            imagesc(handles.axes_OCT_Flow, V.OCT_Flow_XZ(:,:,ll), [im0.gui.colorbalLowerLim im0.gui.colorbarUpperLim])
        else
            im = zeros(size(V.OCT_Flow_XZ, 1), size(V.OCT_Flow_XZ, 2));
            imagesc(handles.axes_OCT_Flow, im, [im0.gui.colorbalLowerLim im0.gui.colorbarUpperLim])
        end
end

xlim(handles.axes_OCT_Flow, xLimit)
ylim(handles.axes_OCT_Flow, yLimit)
yticks(handles.axes_OCT_Flow, [])
colormap(handles.axes_OCT_Flow, im0.gui.flowClrMp);
handles.axes_OCT_Flow.Color = 'k';
enableDefaultInteractivity(handles.axes_OCT_Flow)
end


function update_Colorbar(handles)
global im0
colormap(handles.axes_Colorbar, im0.gui.flowClrMp);
caxis(handles.axes_Colorbar, [im0.gui.colorbalLowerLim im0.gui.colorbarUpperLim])
c = colorbar(handles.axes_Colorbar);
c.Label.String = 'mm/sec';
end

function overlayGraph(handles)
global im0 V
if ~(im0.gui.graphDispState(1)||im0.gui.graphDispState(2))
    %--- Link axes limits
    linkaxes([handles.axes_TP_Ang, handles.axes_OCT_Flow],'xy')
    %--- Set axes_TP_Ang button down function
    fhandle1='OCT_Flow_Mapping(''axes_TP_Ang_ButtonDownFcn'',gcbo,[],guidata(gcbo))';
    set(handles.axes_TP_Ang,'ButtonDownFcn',fhandle1);
    set(get(handles.axes_TP_Ang,'children'), 'ButtonDownFcn',fhandle1);
    
    handles.axes_TP_Ang.Color = 'k';
    enableDefaultInteractivity(handles.axes_TP_Ang)
    
    %--- Set axes_OCT_Flow button down function
    fhandle1='OCT_Flow_Mapping(''axes_OCT_Flow_ButtonDownFcn'',gcbo,[],guidata(gcbo))';
    set(handles.axes_OCT_Flow,'ButtonDownFcn',fhandle1);
    set(get(handles.axes_OCT_Flow,'children'), 'ButtonDownFcn',fhandle1);
    
    handles.axes_OCT_Flow.Color = 'k';
    enableDefaultInteractivity(handles.axes_OCT_Flow)
    return
end

f0 = im0.gui.frameNum - floor((im0.gui.numOfFrames - 1)/2);
f1 = im0.gui.frameNum + ceil((im0.gui.numOfFrames - 1)/2);




switch im0.gui.viewPlane
    case 'XY'
        i = V.Z(1,:) >= f0 & V.Z(1,:) <= f1 & V.Z(2,:) >= f0 & V.Z(2,:) <= f1;
        flags = im0.edgeFlag(i);
        plotX = V.X(:,i); plotY = V.Y(:,i);
        
    case 'YZ'
        i = V.X(1,:) >= f0 & V.X(1,:) <= f1 & V.X(2,:) >= f0 & V.X(2,:) <= f1;
        flags = im0.edgeFlag(i);
        plotX = ((im0.szOCT(1) + 1) - V.Y(:,i)); plotY = V.Z(:,i);
        
    case 'XZ'
        i = V.Y(1,:) >= f0 & V.Y(1,:) <= f1 & V.Y(2,:) >= f0 & V.Y(2,:) <= f1;
        flags = im0.edgeFlag(i);
        plotX = ((im0.szOCT(1) + 1) - V.X(:,i)); plotY = V.Z(:,i);
end

V.projVctrs = V.vctrs(i',:);





flagUp = flags == 1;
flagNeut = flags == 0;
flagDown = flags == -1;
flagTemp = flags == 2;

Xup = plotX(:,flagUp);  Yup = plotY(:,flagUp);
Xneut = plotX(:,flagNeut);  Yneut = plotY(:,flagNeut);
Xdown = plotX(:,flagDown);  Ydown = plotY(:,flagDown);
Xtemp = plotX(:,flagTemp);  Ytemp = plotY(:,flagTemp);

%%%%%%%%%%
if im0.gui.graphDispState(1)
    hold(handles.axes_TP_Ang, 'on')
    plot(handles.axes_TP_Ang, Xneut(:), Yneut(:), 'LineWidth', im0.gui.edgeLinewidth_axes,...
        'color', 'g','Marker','.', 'MarkerSize', im0.gui.nodePointWidth_axes,...
        'MarkerEdgeColor', 'r');
    plot(handles.axes_TP_Ang, Xup(:), Yup(:), 'LineWidth', im0.gui.edgeLinewidth_axes,...
        'color', 'r','Marker','.', 'MarkerSize', im0.gui.nodePointWidth_axes,...
        'MarkerEdgeColor', 'r');
    plot(handles.axes_TP_Ang, Xdown(:), Ydown(:), 'LineWidth', im0.gui.edgeLinewidth_axes,...
        'color', 'b','Marker','.', 'MarkerSize', im0.gui.nodePointWidth_axes,...
        'MarkerEdgeColor', 'r');
    plot(handles.axes_TP_Ang, Xtemp(:), Ytemp(:), 'LineWidth', im0.gui.edgeLinewidth_axes,...
        'color', 'm','Marker','.', 'MarkerSize', im0.gui.nodePointWidth_axes,...
        'MarkerEdgeColor', 'r');
    hold(handles.axes_TP_Ang, 'off')
end

%%-------------------------------------------------------------------------
if im0.gui.graphDispState(2)
    hold(handles.axes_OCT_Flow, 'on')
    plot(handles.axes_OCT_Flow, Xneut(:), Yneut(:), 'LineWidth', im0.gui.edgeLinewidth_axes,...
        'color', 'g','Marker','.', 'MarkerSize', im0.gui.nodePointWidth_axes,...
        'MarkerEdgeColor', 'r');
    plot(handles.axes_OCT_Flow, Xup(:), Yup(:), 'LineWidth', im0.gui.edgeLinewidth_axes,...
        'color', 'r','Marker','.', 'MarkerSize', im0.gui.nodePointWidth_axes,...
        'MarkerEdgeColor', 'r');
    plot(handles.axes_OCT_Flow, Xdown(:), Ydown(:), 'LineWidth', im0.gui.edgeLinewidth_axes,...
        'color', 'b','Marker','.', 'MarkerSize', im0.gui.nodePointWidth_axes,...
        'MarkerEdgeColor', 'r');
    plot(handles.axes_OCT_Flow, Xtemp(:), Ytemp(:), 'LineWidth', im0.gui.edgeLinewidth_axes,...
        'color', 'w','Marker','.', 'MarkerSize', im0.gui.nodePointWidth_axes,...
        'MarkerEdgeColor', 'r');
    hold(handles.axes_OCT_Flow, 'off')
    yticks(handles.axes_OCT_Flow, [])
end


%--- Link axes limits
linkaxes([handles.axes_TP_Ang, handles.axes_OCT_Flow], 'xy')
%--- Set axes_TP_Ang button down function
fhandle1= 'OCT_Flow_Mapping(''axes_TP_Ang_ButtonDownFcn'',gcbo,[],guidata(gcbo))';
set(handles.axes_TP_Ang,'ButtonDownFcn',fhandle1);
set(get(handles.axes_TP_Ang,'children'), 'ButtonDownFcn',fhandle1);

handles.axes_TP_Ang.Color = 'k';
enableDefaultInteractivity(handles.axes_TP_Ang)

%--- Set axes_OCT_Flow button down function
fhandle1= 'OCT_Flow_Mapping(''axes_OCT_Flow_ButtonDownFcn'',gcbo,[],guidata(gcbo))';
set(handles.axes_OCT_Flow,'ButtonDownFcn',fhandle1);
set(get(handles.axes_OCT_Flow,'children'), 'ButtonDownFcn',fhandle1);

handles.axes_OCT_Flow.Color = 'k';
enableDefaultInteractivity(handles.axes_OCT_Flow)
end



function update_Monitors(handles, edge, lgc_update)
%%
global im0 V
cylMode = get(handles.cylMode, 'value');
inspectMode = get(handles.inspectMode, 'value');
if (~cylMode && ~inspectMode) || ~lgc_update; return; end

if cylMode && ~isempty(edge)
    
    p0 = im0.nodeEdges(edge, 1);    p0 = im0.nodePos(p0,:);
    p1 = im0.nodeEdges(edge, 2);    p1 = im0.nodePos(p1,:);
    midP = round((p0 + p1) / 2);
    c0 = midP - floor((im0.gui.winSize - 1) / 2);
    c1 = midP + ceil((im0.gui.winSize - 1) / 2);
    
    sz = size(im0.OCT_Flow);
    
    p0_cyl = im0.cylNodes(edge, 1:3); p1_cyl = im0.cylNodes(edge, 4:6);
    midP = round((p0_cyl + p1_cyl)/2);
    midP =  midP - floor(im0.gui.numOfFrames/2);
    midP(midP < 1) = 1;
    r = im0.cylDia(edge)/2;
    
    [c0_cyl, c1_cyl, mat_cyl] = cylVis(p0_cyl, p1_cyl, r);
    global s
    s.c0_cyl = c0_cyl;
    s.c1_cyl = c1_cyl;
    s.mat_cyl = mat_cyl;


    XY_cyl = squeeze(max(mat_cyl,[],3));
    YZ_cyl = squeeze(max(mat_cyl,[],2));
    XZ_cyl = squeeze(max(mat_cyl,[],1));
    
    im_cyl_XY = zeros(sz([1,2]));
    im_cyl_XY(c0_cyl(2):c1_cyl(2), c0_cyl(1):c1_cyl(1)) = XY_cyl;
    
    im_cyl_YZ = zeros(sz([3, 2]));
    im_cyl_YZ(c0_cyl(3):c1_cyl(3), ((im0.szOCT(1) + 1)-c1_cyl(2)):((im0.szOCT(1) + 1)-c0_cyl(2))) = fliplr(YZ_cyl');
    
    im_cyl_XZ = zeros(sz([3, 1]));
    im_cyl_XZ(c0_cyl(3):c1_cyl(3), ((im0.szOCT(1) + 1)-c1_cyl(1)):((im0.szOCT(1) + 1)-c0_cyl(1))) = fliplr(XZ_cyl');
    
    
    %%% window_01 --XY
    nS = midP(3);
    if nS > V.prjVolSz_OCT(3); nS = V.prjVolSz_OCT(3); end
    imagesc(handles.axes_OCT_Win_01, V.OCT_Flow_XY(:,:,nS),...
        [im0.gui.colorbalLowerLim im0.gui.colorbarUpperLim])
    
    hold(handles.axes_OCT_Win_01, 'on')
    im = imagesc(handles.axes_OCT_Win_01, im_cyl_XY);
    im.AlphaData = (max(im_cyl_XY,[],3)/999) * im0.gui.alpha;
    
    xlim(handles.axes_OCT_Win_01, [c0(1) c1(1)])
    ylim(handles.axes_OCT_Win_01, [c0(2) c1(2)])
    colormap(handles.axes_OCT_Win_01, im0.gui.flowClrMp);
    xticks(handles.axes_OCT_Win_01, [])
    yticks(handles.axes_OCT_Win_01, [])
    handles.axes_OCT_Win_01.Color = 'k';
    
    plot(handles.axes_OCT_Win_01, [p0(1), p1(1)], [p0(2), p1(2)], 'LineWidth', im0.gui.edgeLinewidth_monitors,...
        'color', 'g','Marker','.', 'MarkerSize', im0.gui.nodePointWidth_monitors,...
        'MarkerEdgeColor', 'r');
    hold(handles.axes_OCT_Win_01, 'off')
    
    
    
    %%% window_02 --YZ
    nS = midP(1);
    if nS > V.prjVolSz_OCT(1); nS = V.prjVolSz_OCT(1); end
    
    imagesc(handles.axes_OCT_Win_02, V.OCT_Flow_YZ(:,:,nS),...
        [im0.gui.colorbalLowerLim im0.gui.colorbarUpperLim])
    
    hold(handles.axes_OCT_Win_02, 'on')
    im = imagesc(handles.axes_OCT_Win_02, im_cyl_YZ);
    im.AlphaData = (max(im_cyl_YZ,[],3)/999) * im0.gui.alpha;
    
    xlim(handles.axes_OCT_Win_02, [((im0.szOCT(1) + 1) - c1(2)) ((im0.szOCT(1) + 1) - c0(2))])
    ylim(handles.axes_OCT_Win_02, [c0(3) c1(3)])
    colormap(handles.axes_OCT_Win_02, im0.gui.flowClrMp);
    xticks(handles.axes_OCT_Win_02, [])
    yticks(handles.axes_OCT_Win_02, [])
    handles.axes_OCT_Win_02.Color = 'k';
    
    plot(handles.axes_OCT_Win_02, [((im0.szOCT(1) + 1)-p0(2)), ((im0.szOCT(1) + 1)-p1(2))], [p0(3), p1(3)], 'LineWidth', im0.gui.edgeLinewidth_monitors,...
        'color', 'g','Marker','.', 'MarkerSize', im0.gui.nodePointWidth_monitors,...
        'MarkerEdgeColor', 'r');
    hold(handles.axes_OCT_Win_02, 'off')
    
    %%% window_03 --XZ
    nS = midP(2);
    if nS > V.prjVolSz_OCT(2); nS = V.prjVolSz_OCT(2); end
    
    imagesc(handles.axes_OCT_Win_03, V.OCT_Flow_XZ(:,:,nS),...
        [im0.gui.colorbalLowerLim im0.gui.colorbarUpperLim])
    
    hold(handles.axes_OCT_Win_03, 'on')
    im = imagesc(handles.axes_OCT_Win_03, im_cyl_XZ);
    im.AlphaData = (max(im_cyl_XZ,[],3)/999) * im0.gui.alpha;
    
    xlim(handles.axes_OCT_Win_03, [((im0.szOCT(1) + 1) - c1(1)) ((im0.szOCT(1) + 1) - c0(1))])
    ylim(handles.axes_OCT_Win_03, [c0(3) c1(3)])
    colormap(handles.axes_OCT_Win_03, im0.gui.flowClrMp);
    xticks(handles.axes_OCT_Win_03, [])
    yticks(handles.axes_OCT_Win_03, [])
    handles.axes_OCT_Win_03.Color = 'k';
    
    plot(handles.axes_OCT_Win_03, [((im0.szOCT(1) + 1)-p0(1)), ((im0.szOCT(1) + 1)-p1(1))], [p0(3), p1(3)], 'LineWidth', im0.gui.edgeLinewidth_monitors,...
        'color', 'g','Marker','.', 'MarkerSize', im0.gui.nodePointWidth_monitors,...
        'MarkerEdgeColor', 'r');
    hold(handles.axes_OCT_Win_03, 'off')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %%% wind_TPAng
    switch im0.gui.angWind
        case 'XY'
            nS = midP(3);
            if nS > V.prjVolSz_TPM(3); nS = V.prjVolSz_TPM(3); end
            imagesc(handles.axes_TPAng_win, V.TP_XY(:,:,nS), [im0.gui.TPthreshold*32 32])
            xlim(handles.axes_TPAng_win, [c0(1) c1(1)])
            ylim(handles.axes_TPAng_win, [c0(2) c1(2)])
            xticks(handles.axes_TPAng_win, [])
            yticks(handles.axes_TPAng_win, [])
            colormap(handles.axes_TPAng_win, gray);
            
            hold(handles.axes_TPAng_win, 'on')
            plot(handles.axes_TPAng_win, [p0(1), p1(1)], [p0(2), p1(2)], 'LineWidth', im0.gui.edgeLinewidth_monitors,...
                'color', 'g','Marker','.', 'MarkerSize', im0.gui.nodePointWidth_monitors,...
                'MarkerEdgeColor', 'r');
            hold(handles.axes_TPAng_win, 'off')
            
            
            im = imagesc(handles.axes_ang_cyl, im_cyl_XY);
            im.AlphaData = (max(im_cyl_XY,[],3)/999) * im0.gui.alpha;
            colormap(handles.axes_ang_cyl, [1 0 0])
            handles.axes_ang_cyl.Visible = 'off';
            xlim(handles.axes_ang_cyl, [c0(1) c1(1)])
            ylim(handles.axes_ang_cyl, [c0(2) c1(2)])
            linkaxes([handles.axes_ang_cyl, handles.axes_TPAng_win],'xy')
            xticks(handles.axes_ang_cyl, [])
            yticks(handles.axes_ang_cyl, [])
            handles.axes_TPAng_win.Color = 'k';
            
        case 'YZ'
            nS = midP(1);
            if nS > V.prjVolSz_TPM(1); nS = V.prjVolSz_TPM(1); end
            imagesc(handles.axes_TPAng_win, V.TP_YZ(:,:,nS), [im0.gui.TPthreshold*32 32])
            xlim(handles.axes_TPAng_win, [((im0.szOCT(1) + 1) - c1(2)) ((im0.szOCT(1) + 1) - c0(2))])
            ylim(handles.axes_TPAng_win, [c0(3) c1(3)])
            xticks(handles.axes_TPAng_win, [])
            yticks(handles.axes_TPAng_win, [])
            colormap(handles.axes_TPAng_win, gray);
            
            hold(handles.axes_TPAng_win, 'on')
            plot(handles.axes_TPAng_win, [((im0.szOCT(1) + 1)-p0(2)), ((im0.szOCT(1) + 1)-p1(2))], [p0(3), p1(3)], 'LineWidth', im0.gui.edgeLinewidth_monitors,...
                'color', 'g','Marker','.', 'MarkerSize', im0.gui.nodePointWidth_monitors,...
                'MarkerEdgeColor', 'r');
            hold(handles.axes_TPAng_win, 'off')
            
            
            im = imagesc(handles.axes_ang_cyl, im_cyl_YZ);
            im.AlphaData = (max(im_cyl_YZ,[],3)/999) * im0.gui.alpha;
            colormap(handles.axes_ang_cyl, [1 0 0])
            handles.axes_ang_cyl.Visible = 'off';
            xlim(handles.axes_ang_cyl, [((im0.szOCT(1) + 1) - c1(2)) ((im0.szOCT(1) + 1) - c0(2))])
            ylim(handles.axes_ang_cyl, [c0(3) c1(3)])
            linkaxes([handles.axes_ang_cyl, handles.axes_TPAng_win],'xy')
            xticks(handles.axes_ang_cyl, [])
            yticks(handles.axes_ang_cyl, [])
            handles.axes_TPAng_win.Color = 'k';
            
        case 'XZ'
            nS = midP(2);
            if nS > V.prjVolSz_TPM(2); nS = V.prjVolSz_TPM(2); end
            imagesc(handles.axes_TPAng_win, V.TP_XZ(:,:,nS), [im0.gui.TPthreshold*32 32])
            xlim(handles.axes_TPAng_win, [((im0.szOCT(1) + 1) - c1(1)) ((im0.szOCT(1) + 1) - c0(1))])
            ylim(handles.axes_TPAng_win, [c0(3) c1(3)])
            xticks(handles.axes_TPAng_win, [])
            yticks(handles.axes_TPAng_win, [])
            colormap(handles.axes_TPAng_win, gray);
            
            hold(handles.axes_TPAng_win, 'on')
            plot(handles.axes_TPAng_win, [((im0.szOCT(1) + 1)-p0(1)), ((im0.szOCT(1) + 1)-p1(1))], [p0(3), p1(3)], 'LineWidth', im0.gui.edgeLinewidth_monitors,...
                'color', 'g','Marker','.', 'MarkerSize', im0.gui.nodePointWidth_monitors,...
                'MarkerEdgeColor', 'r');
            hold(handles.axes_TPAng_win, 'off')
            
            
            im = imagesc(handles.axes_ang_cyl, im_cyl_XZ);
            im.AlphaData = (max(im_cyl_XZ,[],3)/999) * im0.gui.alpha;
            colormap(handles.axes_ang_cyl, [1 0 0])
            handles.axes_ang_cyl.Visible = 'off';
            xlim(handles.axes_ang_cyl, [((im0.szOCT(1) + 1) - c1(1)) ((im0.szOCT(1) + 1) - c0(1))])
            ylim(handles.axes_ang_cyl, [c0(3) c1(3)])
            linkaxes([handles.axes_ang_cyl, handles.axes_TPAng_win],'xy')
            xticks(handles.axes_ang_cyl, [])
            yticks(handles.axes_ang_cyl, [])
            handles.axes_TPAng_win.Color = 'k';
    end
    
    
elseif inspectMode && ~isempty(im0.gui.mntrPt)
    midP = im0.gui.mntrPt - ceil(im0.gui.numOfFrames/2) + 1;
    c0 = midP - floor((im0.gui.winSize - 1) / 2);
    c1 = midP + ceil((im0.gui.winSize - 1) / 2);
    midP(midP < 1) = 1;
    
    %%% window_01 --XY
    if midP(3) <= size(V.OCT_Flow_XY, 3)
        imagesc(handles.axes_OCT_Win_01, V.OCT_Flow_XY(:,:,midP(3)),...
            [im0.gui.colorbalLowerLim im0.gui.colorbarUpperLim])
    else
        imagesc(handles.axes_OCT_Win_01, zeros(im0.szOCT(1), im0.szOCT(2)),...
            [im0.gui.colorbalLowerLim im0.gui.colorbarUpperLim])
    end
    
    xlim(handles.axes_OCT_Win_01, [c0(1) c1(1)])
    ylim(handles.axes_OCT_Win_01, [c0(2) c1(2)])
    colormap(handles.axes_OCT_Win_01, im0.gui.flowClrMp);
    xticks(handles.axes_OCT_Win_01, [])
    yticks(handles.axes_OCT_Win_01, [])
    handles.axes_OCT_Win_01.Color = 'k';
    
    %%% window_02 --YZ
    if midP(1) <= size(V.OCT_Flow_YZ, 3)
        imagesc(handles.axes_OCT_Win_02, V.OCT_Flow_YZ(:,:,midP(1)),...
            [im0.gui.colorbalLowerLim im0.gui.colorbarUpperLim])
    else
        imagesc(handles.axes_OCT_Win_02, zeros(im0.szOCT(2), im0.szOCT(3)),...
            [im0.gui.colorbalLowerLim im0.gui.colorbarUpperLim])
    end
    xlim(handles.axes_OCT_Win_02, [((im0.szOCT(1) + 1) - c1(2)) ((im0.szOCT(1) + 1) - c0(2))])
    ylim(handles.axes_OCT_Win_02, [c0(3) c1(3)])
    colormap(handles.axes_OCT_Win_02, im0.gui.flowClrMp);
    xticks(handles.axes_OCT_Win_02, [])
    yticks(handles.axes_OCT_Win_02, [])
    handles.axes_OCT_Win_02.Color = 'k';
    
    %%% window_03 --XZ
    if midP(2) <= size(V.OCT_Flow_XZ, 3)
        imagesc(handles.axes_OCT_Win_03, V.OCT_Flow_XZ(:,:,midP(2)),...
            [im0.gui.colorbalLowerLim im0.gui.colorbarUpperLim])
    else
        imagesc(handles.axes_OCT_Win_03, zeros(im0.szOCT(1), im0.szOCT(3)),...
            [im0.gui.colorbalLowerLim im0.gui.colorbarUpperLim])
    end
    
    xlim(handles.axes_OCT_Win_03, [((im0.szOCT(1) + 1) - c1(1)) ((im0.szOCT(1) + 1) - c0(1))])
    ylim(handles.axes_OCT_Win_03, [c0(3) c1(3)])
    colormap(handles.axes_OCT_Win_03, im0.gui.flowClrMp);
    xticks(handles.axes_OCT_Win_03, [])
    yticks(handles.axes_OCT_Win_03, [])
    handles.axes_OCT_Win_03.Color = 'k';
    
    %%% wind_TPAng
    switch im0.gui.angWind
        case 'XY'
            imagesc(handles.axes_TPAng_win, V.TP_XY(:,:,midP(3)), [im0.gui.TPthreshold*32 32])
            xlim(handles.axes_TPAng_win, [c0(1) c1(1)])
            ylim(handles.axes_TPAng_win, [c0(2) c1(2)])
            xticks(handles.axes_TPAng_win, [])
            yticks(handles.axes_TPAng_win, [])
            colormap(handles.axes_TPAng_win, gray);
            handles.axes_TPAng_win.Color = 'k';
            
            imagesc(handles.axes_ang_cyl, []);
            handles.axes_ang_cyl.Visible = 'off';
            xlim(handles.axes_ang_cyl, [c0(1) c1(1)])
            ylim(handles.axes_ang_cyl, [c0(2) c1(2)])
            linkaxes([handles.axes_ang_cyl, handles.axes_TPAng_win],'xy')
            xticks(handles.axes_ang_cyl, [])
            yticks(handles.axes_ang_cyl, [])
            
        case 'YZ'
            imagesc(handles.axes_TPAng_win, V.TP_YZ(:,:,midP(1)), [im0.gui.TPthreshold*32 32])
            xlim(handles.axes_TPAng_win, [((im0.szOCT(1) + 1) - c1(2)) ((im0.szOCT(1) + 1) - c0(2))])
            ylim(handles.axes_TPAng_win, [c0(3) c1(3)])
            xticks(handles.axes_TPAng_win, [])
            yticks(handles.axes_TPAng_win, [])
            colormap(handles.axes_TPAng_win, gray);
            
            imagesc(handles.axes_ang_cyl, []);
            handles.axes_ang_cyl.Visible = 'off';
            xlim(handles.axes_ang_cyl, [((im0.szOCT(1) + 1) - c1(2)) ((im0.szOCT(1) + 1) - c0(2))])
            ylim(handles.axes_ang_cyl, [c0(3) c1(3)])
            linkaxes([handles.axes_ang_cyl, handles.axes_TPAng_win],'xy')
            xticks(handles.axes_ang_cyl, [])
            yticks(handles.axes_ang_cyl, [])
            handles.axes_TPAng_win.Color = 'k';
            
            
        case 'XZ'
            imagesc(handles.axes_TPAng_win, V.TP_XZ(:,:,midP(2)), [im0.gui.TPthreshold*32 32])
            xlim(handles.axes_TPAng_win, [((im0.szOCT(1) + 1) - c1(1)) ((im0.szOCT(1) + 1) - c0(1))])
            ylim(handles.axes_TPAng_win, [c0(3) c1(3)])
            xticks(handles.axes_TPAng_win, [])
            yticks(handles.axes_TPAng_win, [])
            colormap(handles.axes_TPAng_win, gray);
            
            imagesc(handles.axes_ang_cyl, []);
            handles.axes_ang_cyl.Visible = 'off';
            xlim(handles.axes_ang_cyl, [((im0.szOCT(1) + 1) - c1(1)) ((im0.szOCT(1) + 1) - c0(1))])
            ylim(handles.axes_ang_cyl, [c0(3) c1(3)])
            linkaxes([handles.axes_ang_cyl, handles.axes_TPAng_win],'xy')
            xticks(handles.axes_ang_cyl, [])
            yticks(handles.axes_ang_cyl, [])
            handles.axes_TPAng_win.Color = 'k';
    end
    
else
    %%-- window_01
    imagesc(handles.axes_OCT_Win_01, zeros(im0.gui.winSize),...
        [im0.gui.colorbalLowerLim im0.gui.colorbarUpperLim])
    colormap(handles.axes_OCT_Win_01, im0.gui.flowClrMp);
    xticks(handles.axes_OCT_Win_01, [])
    yticks(handles.axes_OCT_Win_01, [])
    handles.axes_OCT_Win_01.Color = 'k';
    
    %%-- window_02
    imagesc(handles.axes_OCT_Win_02, zeros(im0.gui.winSize),...
        [im0.gui.colorbalLowerLim im0.gui.colorbarUpperLim])
    colormap(handles.axes_OCT_Win_02, im0.gui.flowClrMp);
    xticks(handles.axes_OCT_Win_02, [])
    yticks(handles.axes_OCT_Win_02, [])
    handles.axes_OCT_Win_02.Color = 'k';
    
    %%-- window_03
    imagesc(handles.axes_OCT_Win_03, zeros(im0.gui.winSize),...
        [im0.gui.colorbalLowerLim im0.gui.colorbarUpperLim])
    colormap(handles.axes_OCT_Win_03, im0.gui.flowClrMp);
    xticks(handles.axes_OCT_Win_03, [])
    yticks(handles.axes_OCT_Win_03, [])
    handles.axes_OCT_Win_03.Color = 'k';
    
    %%-- wind_TPAng
    imagesc(handles.axes_TPAng_win, zeros(im0.gui.winSize),...
        [im0.gui.TPthreshold*32 32])
    colormap(handles.axes_TPAng_win, gray);
    xticks(handles.axes_TPAng_win, [])
    yticks(handles.axes_TPAng_win, [])
    handles.axes_TPAng_win.Color = 'k';
end


fHandle1 = 'OCT_Flow_Mapping(''monitorsBtnDownFnc'',gcbo,[],guidata(gcbo))';
set(handles.axes_OCT_Win_01,'ButtonDownFcn', fHandle1);
set(get(handles.axes_OCT_Win_01, 'children'), 'ButtonDownFcn', fHandle1);
set(handles.axes_OCT_Win_02,'ButtonDownFcn', fHandle1);
set(get(handles.axes_OCT_Win_02, 'children'), 'ButtonDownFcn', fHandle1);
set(handles.axes_OCT_Win_03,'ButtonDownFcn', fHandle1);
set(get(handles.axes_OCT_Win_03, 'children'), 'ButtonDownFcn', fHandle1);
handles.axes_OCT_Win_01.Tag = 'axes_OCT_Win_01.Tag';
handles.axes_OCT_Win_02.Tag = 'axes_OCT_Win_02.Tag';
handles.axes_OCT_Win_03.Tag = 'axes_OCT_Win_03.Tag';
end





function [c0, c1, mat] = cylVis(p0, p1, r)
global im0


x0 = min(p0(1),p1(1)) - r/im0.Hvox(1);
x1 = max(p0(1),p1(1)) + r/im0.Hvox(1);
y0 = min(p0(2),p1(2)) - r/im0.Hvox(2);
y1 = max(p0(2),p1(2)) + r/im0.Hvox(2);
z0 = min(p0(3),p1(3)) - r/im0.Hvox(3);
z1 = max(p0(3),p1(3)) + r/im0.Hvox(3);
c0 = floor([x0,y0,z0]);
c1 = ceil([x1,y1,z1]);
Vsize = c1 - c0 + 1;
vSMax = max(Vsize(1), Vsize(2));
Vsize = [vSMax vSMax Vsize(3)];
c1 = c0 + Vsize - 1;

sz = size(im0.I);

if c0(1) < 1; c0(1) = 1; c1(1) = Vsize(1); end
if c0(2) < 1; c0(2) = 1; c1(2) = Vsize(2); end
if c0(3) < 1; c0(3) = 1; c1(3) = Vsize(3); end

if c1(1) > sz(1); c1(1) = sz(1); c0(1) = c1(1) - Vsize(1) + 1; end
if c1(2) > sz(2); c1(2) = sz(2); c0(2) = c1(2) - Vsize(2) + 1; end
if c1(3) > sz(3); c1(3) = sz(3); c0(3) = c1(3) - Vsize(3) + 1; end

mat = zeros(Vsize);

v = p1 - p0;       %edge vector
vScld = (v.*im0.Hvox);  %scaled edge vector
vScldMag = norm(vScld);

[X, Y, Z] = meshgrid(c0(1):c1(1), c0(2):c1(2), c0(3):c1(3));
X = (X - p0(1))* im0.Hvox(1);
Y = (Y - p0(2))* im0.Hvox(2);
Z = (Z - p0(3))* im0.Hvox(3);

d = sqrt((Y*vScld(3) - Z*vScld(2)).^2 + (Z*vScld(1) - X*vScld(3)).^2 + (X*vScld(2) - Y*vScld(1)).^2)...
    ./ vScldMag;
v0Prj = (X*vScld(1) + Y*vScld(2) + Z*vScld(3))/vScldMag;

mat(d<=r & v0Prj<=vScldMag & v0Prj>=0) = 999;
end

% --- Executes on mouse press over axes background.
function axes_TP_Ang_ButtonDownFcn(~, ~, handles)
axesButtonDownFcn(handles)
end

% --- Executes on mouse press over axes background.
function axes_OCT_Flow_ButtonDownFcn(~, ~, handles) 
axesButtonDownFcn(handles)
end

function axesButtonDownFcn(handles)
global V im0
slctType = get(gcf, 'SelectionType');

switch slctType
    case 'normal'
        pos = get(gca, 'CurrentPoint');
        p = round([pos(1,1), pos(1,2)]);
        cylMode = get(handles.cylMode, 'value');
        
        if cylMode
            e = selectedEdgeIndx(V.projVctrs, p, 4);
            if ~isempty(e)
                im0.gui.edgeID = e;
                im0.gui.segmentID = im0.edgeSegN(e);
                updateGui_edgeInfo(handles)
                
                
                computeEdgeFlow(handles, e);
                
                if ~isempty(im0.selectedEdgesList) &&...
                        ~isempty(intersect(im0.selectedEdgesList(2,:), e))
                    l = find(im0.selectedEdgesList(2,:)== e);
                    set(handles.selectedEdgesList, 'Value', l);
                end
                
                
                pt = im0.cylNodes(e, :);
                im0.gui.mntrPt = round((pt(1:3) + pt(4:6)) / 2);
                
                im0.gui.positionFlag(2) = true;
                highlightEdges(handles, e)
                
            end
            
        else
            
            switch im0.gui.viewPlane
                case 'XY'
                    im0.gui.mntrPt([1,2]) = p;
                case 'YZ'
                    im0.gui.mntrPt(3) = p(2);
                    im0.gui.mntrPt(2) = (im0.szOCT(1) + 1) - p(1);
                case 'XZ'
                    im0.gui.mntrPt(3) = p(2);
                    im0.gui.mntrPt(1) = (im0.szOCT(1) + 1) - p(1);
            end
            
            im0.gui.positionFlag(2) = true;
            if all(im0.gui.positionFlag)
                update_Monitors(handles, im0.gui.edgeID, true)
                update_TPMA_image(handles)
                update_OCTV_image(handles)
                positionFlag_OCTFlow(handles, im0.gui.mntrPt)
                positionFlag_TPAng(handles, im0.gui.mntrPt)
                overlayGraph(handles)
            else
                update_Monitors(handles, im0.gui.edgeID, true)
            end
        end
        
    case 'alt'
        axP = get(gca,'Position');
        p = get(gcf,'CurrentPoint');
        p = ((p - axP(1:2)) .* ([500, 500]./axP(3:4)));
        if axP(1) > 20
            p = p + [622, 389];
        else
            p = p + [98, 389];
        end
        
        handles.axes_TP_Ang.UIContextMenu.Position = p;
        handles.axes_TP_Ang.UIContextMenu.Visible = 'on';
end
end


%%-------------------------------------------------------------------------
function [XX, YY, ZZ] = edgeGraph(vNodes, nodePos)

x0 = nodePos(vNodes(:,1),1);
x1 = nodePos(vNodes(:,2),1);
y0 = nodePos(vNodes(:,1),2);
y1 = nodePos(vNodes(:,2),2);
z0 = nodePos(vNodes(:,1),3);
z1 = nodePos(vNodes(:,2),3);

XX = [x0'; x1']; XX(3,:) = nan;
YY = [y0'; y1']; YY(3,:) = nan;
ZZ = [z0'; z1']; ZZ(3,:) = nan;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function outVol = intensProj(vol, numFrames, plane, type)
% intensProj() - computes the intensity projection for each plane in a volume
%
% Syntax: outVol = intensProj(vol, numFrames, axis, type)
%
% Inputs:
%   vol - 3D matrix
%   numFrames - number of planes from which the projection is computed
%   axis - ('X','Y' or 'Z') The axis along which the projection is computed
%   type - ('max' or 'min') projection
%
% Outputs:
%   outVol - 3D matrix of all projectons

s = size(vol);  %% volume size

switch plane
    case 'XY'
        outVol = zeros( s(1), s(2), (s(3)-numFrames+1) );
        for i = 1:(s(3)-numFrames+1)  %% index of first frame in projection vol
            f = i+numFrames-1; %% index of last frame in projection vol
            switch type
                case 'max'
                    maxP = squeeze(max(vol(:,:,i:f),[],3));
                    minP = squeeze(min(vol(:,:,i:f),[],3));
                    minP_abs = abs(minP);
                    maxP(minP_abs(:) > maxP(:)) = minP(minP_abs(:) > maxP(:));
                    outVol(:,:,i) = maxP;
                case 'min'
                    maxP = squeeze(max(vol(:,:,i:f),[],3));
                    minP = squeeze(min(vol(:,:,i:f),[],3));
                    minP_abs = abs(minP);
                    minP(minP_abs(:) > maxP(:)) = maxP(minP_abs(:) > maxP(:));
                    outVol(:,:,i) = minP;
                case 'mean'
                    meanP = squeeze(mean(vol(:,:,i:f),3));
                    outVol(:,:,i) = meanP;
            end
        end
        if s(3) == numFrames
            outVol = squeeze(outVol);
        end
        %%%%%%%%%%%
        
    case 'YZ'
        outVol = zeros(s(3), s(1), (s(2)-numFrames+1));
        for i = 1:(s(2)-numFrames+1)  %% index of first frame in projection vol
            f = i+numFrames-1; %% index of last frame in projection vol
            switch type
                case 'max'
                    maxP = squeeze(max(vol(:,i:f,:),[],2));
                    minP = squeeze(min(vol(:,i:f,:),[],2));
                    minP_abs = abs(minP);
                    maxP(minP_abs(:) > maxP(:)) = minP(minP_abs(:) > maxP(:));
                    %outVol(:,i,:) = maxP;
                    outVol(:,:,i) = rot90(maxP, -1);
                case 'min'
                    maxP = squeeze(max(vol(:,i:f,:),[],2));
                    minP = squeeze(min(vol(:,i:f,:),[],2));
                    minP_abs = abs(minP);
                    minP(minP_abs(:) > maxP(:)) = maxP(minP_abs(:) > maxP(:));
                    %outVol(:,i,:) = minP;
                    outVol(:,:,i) = rot90(minP,-1);
                case 'mean'
                    meanP = squeeze(mean(vol(:,i:f,:),2));
                    %outVol(:,i,:) = meanP;
                    outVol(:,:,i) = rot90(meanP,-1);
            end
        end
        if s(2) == numFrames
            outVol = squeeze(outVol);
        end
        %%%%%%%%%%%
        
    case 'XZ'
        outVol = zeros(s(3), s(2), (s(1)-numFrames+1));
        for i = 1:(s(1)-numFrames+1)  %% index of first frame in projection vol
            f = i+numFrames-1; %% index of last frame in projection vol
            switch type
                case 'max' %#ok<*BDSCA>
                    maxP = squeeze(max(vol(i:f,:,:),[],1));
                    minP = squeeze(min(vol(i:f,:,:),[],1));
                    minP_abs = abs(minP);
                    maxP(minP_abs(:) > maxP(:)) = minP(minP_abs(:) > maxP(:));
                    %outVol(i,:,:) = maxP;
                    outVol(:,:,i) = rot90(maxP,-1);
                case 'min'
                    maxP = squeeze(max(vol(i:f,:,:),[],1));
                    minP = squeeze(min(vol(i:f,:,:),[],1));
                    minP_abs = abs(minP);
                    minP(minP_abs(:) > maxP(:)) = maxP(minP_abs(:) > maxP(:));
                    %outVol(i,:,:) = minP;
                    outVol(:,:,i) = rot90(minP,-1);
                case 'mean'
                    meanP = squeeze(mean(vol(i:f,:,:),1));
                    %outVol(i,:,:) = meanP;
                    outVol(:,:,i) = rot90(meanP,-1);
            end
        end
        if s(1) == numFrames
            outVol = squeeze(outVol);
        end
        %%%%%%%%%%%
end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function vctrs = vInfo(vNodes, nodePos, plane)
global im0
%%-- [edgVec edgVec edg0Node edg0Node edgeMag edgeID]
l = length(vNodes);

if plane == 'XY'
    p0 = nodePos(vNodes(:,1),[1,2]);
    p1 = nodePos(vNodes(:,2),[1,2]);
    v = p1 - p0;
    mg = vecnorm(v, 2, 2);
    vctrs = [v, p0, mg, (1:l)'];
    
elseif plane == 'YZ'
    p0 = nodePos(vNodes(:,1),[2,3]); p0(:,1) = (im0.szOCT(1) + 1) - p0(:,1);
    p1 = nodePos(vNodes(:,2),[2,3]); p1(:,1) = (im0.szOCT(1) + 1) - p1(:,1);
    v = p1 - p0;
    mg = vecnorm(v, 2, 2);
    vctrs = [v,p0,mg,(1:l)'];
    
elseif plane == 'XZ'
    p0 = nodePos(vNodes(:,1),[1,3]); p0(:,1) = (im0.szOCT(1) + 1) - p0(:,1);
    p1 = nodePos(vNodes(:,2),[1,3]); p1(:,1) = (im0.szOCT(1) + 1) - p1(:,1);
    v = p1 - p0;
    mg = vecnorm(v, 2, 2);
    vctrs = [v,p0,mg,(1:l)'];
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function indx = selectedEdgeIndx(V, currentPoint, range)

rho = currentPoint - V(:,3:4);
d = vecnorm((V(:,1).*rho(:,2) - V(:,2).*rho(:,1)), 2, 2)./ V(:,5);

rhoProj = dot(V(:,1:2),rho(:,1:2),2) ./ V(:,5);
d(:,2) = V(:,6);

i = d(d(:,1) <= range & rhoProj <= V(:,5) & rhoProj >= 0, :);

if length(i) == 1
    indx = i(1,2);
elseif length(i) > 1
    i = i(i(:,1) == min(i(:,1)), 2);
    indx = i(1);
else
    indx = [];
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function computeEdgeFlow(handles, edgeID)

if isempty(edgeID)
    set(handles.numOfDatapoints, 'String', 'Ø');
    set(handles.angleToZ,'String', 'Ø');
    set(handles.edgeFlowV,'String', 'Ø');
    set(handles.maxVelZ,'String', 'Ø');
    set(handles.meanVelZ,'String', 'Ø');
    set(handles.minVelZ,'String', 'Ø');
    set(handles.cylLength,'String', 'Ø');
    set(handles.sgmntMeanFlowV,'String', 'Ø');
    return
end
%-----------------------------------------------------------

global im0

vScld = (im0.nodePos_um(im0.nodeEdges(edgeID, 2), :) -...
    im0.nodePos_um(im0.nodeEdges(edgeID, 1), :));
vScldMag = norm(vScld);
cosPhZ = dot([0,0,1],vScld)/vScldMag;
flowVel = getDataPt(edgeID, 'allButZeros');

%%% Qi Update 2023 March Start
maxFlowVel_Z = [];
minFlowVel_Z = [];
meanFlowVel_Z = [];
stdFlowVel_Z = [];
meanFlowEdge = [];
stdFlowEdge = [];
phZ = [];
%%% Qi Update 2023 March End

if ~isempty(flowVel)
    phZ = acos(abs(cosPhZ)) * 180/pi();
    meanFlowVel_Z = mean(flowVel);
    stdFlowVel_Z = std(flowVel);
    
    if meanFlowVel_Z > 0
        maxFlowVel_Z = max(flowVel(:));
        minFlowVel_Z = min(flowVel(:));
    else
        maxFlowVel_Z = min(flowVel(:));
        minFlowVel_Z = max(flowVel(:));
    end
    ctOff = im0.cutOff(edgeID);
    if cosPhZ ~= 0
        if meanFlowVel_Z > 0
            flo = flowVel(flowVel > 0);
            flo = flo(flo >= ctOff);
        else
            flo = flowVel(flowVel < 0);
            flo = flo(flo <= ctOff);
        end
        if isempty(flo); flo = 0; end
        meanFlowEdge = mean(flo)/abs(cosPhZ);
        stdFlowEdge = std(flo/abs(cosPhZ));
    end
end

%%-------------------------------------------------------------------------
histogram(handles.axes_Histogram, flowVel, 30, 'Orientation', 'horizontal')
enableDefaultInteractivity(handles.axes_Histogram)

if length(flowVel) > 1
    xLm = xlim(handles.axes_Histogram);
    mad = median(flowVel);
    hold(handles.axes_Histogram, 'on')
    switch handles.menu_showMedian.Checked
        case 'on'
    plot(handles.axes_Histogram, [0, xLm(2)],[mad, mad], 'color', 'k',...
        'LineStyle', '--', 'LineWidth', 1)
    end
    switch handles.menu_showCutoff.Checked
        case 'on'
    plot(handles.axes_Histogram, [0, xLm(2)],[ctOff, ctOff], 'color', 'r',...
        'LineStyle', '-', 'LineWidth', 1)
    end
    hold(handles.axes_Histogram, 'off')
end
%%-------------------------------------------------------------------------
im0.gui.edgeFlowV = meanFlowEdge;
if ~isempty(meanFlowEdge)
    im0.edgeMeanFlow(edgeID) = meanFlowEdge;
    %%% Qi Update 2023 March Start
    im0.edgeStdFlow(edgeID) = stdFlowEdge;
    im0.edgeMeanFlowVel_Z(edgeID) = meanFlowVel_Z;
    im0.edgeStdFlowVel_Z(edgeID) = stdFlowVel_Z;
    im0.cosPhZ(edgeID) = cosPhZ;
    im0.phZ(edgeID) = phZ;
    %%% Qi Update 2023 March End
else
    meanFlowEdge = 0;
    %%% Qi Update 2023 March Start
    stdFlowEdge = 0;
    meanFlowVel_Z = 0;
    stdFlowVel_Z = 0;
    %%% Qi Update 2023 March End
    im0.edgeMeanFlow(edgeID) = meanFlowEdge;
    %%% Qi Update 2023 March Start
    im0.edgeStdFlow(edgeID) = stdFlowEdge;
    im0.edgeMeanFlowVel_Z(edgeID) = meanFlowVel_Z;
    im0.edgeStdFlowVel_Z(edgeID) = stdFlowVel_Z;
    im0.cosPhZ(edgeID) = cosPhZ;
    im0.phZ(edgeID) = phZ;
    %%% Qi Update 2023 March End
end

if ~isempty(im0.selectedEdgesList)
    im0.selectedEdgesList(3, im0.selectedEdgesList(2,:) == edgeID) = meanFlowEdge;
    updateGui_Listbox(handles)
    if ~isempty(intersect(im0.selectedEdgesList(1,:), im0.gui.segmentID))
        
        f = im0.selectedEdgesList(3, (im0.selectedEdgesList(1,:) == im0.gui.segmentID...
            & im0.selectedEdgesList(2,:) ~= im0.gui.edgeID));
        f = [f, meanFlowEdge];
        
        SgmntMean = mean(f);
        set(handles.sgmntMeanFlowV,'String',[num2str(SgmntMean), ' mm/s']);
    else
        set(handles.sgmntMeanFlowV,'String', 'Ø');
    end
else
    set(handles.sgmntMeanFlowV,'String', 'Ø');
end

if ~isempty(flowVel)
    set(handles.numOfDatapoints, 'String', num2str(length(flowVel)));
    set(handles.angleToZ,'String', [num2str(phZ),'°']);
    if ~isempty(meanFlowEdge)
        set(handles.edgeFlowV,'String', [num2str(meanFlowEdge), ' mm/s']);
    else
        set(handles.edgeFlowV,'String', 'Ø');
    end
    set(handles.maxVelZ,'String', [num2str(maxFlowVel_Z), ' mm/s']);
    set(handles.meanVelZ,'String', [num2str(meanFlowVel_Z), ' mm/s']);
    set(handles.minVelZ,'String', [num2str(minFlowVel_Z), ' mm/s']);
    set(handles.cylLength,'String', [num2str(vScldMag), ' μm']);
else
    set(handles.numOfDatapoints, 'String', 'Ø');
    set(handles.angleToZ,'String', 'Ø');
    set(handles.edgeFlowV,'String', 'Ø');
    set(handles.maxVelZ,'String', 'Ø');
    set(handles.meanVelZ,'String', 'Ø');
    set(handles.minVelZ,'String', 'Ø');
    set(handles.cylLength,'String', 'Ø');
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function updateSelectedEdgesList(handles)
global im0
if ~isempty(im0.selectedEdgesList)
    e = im0.selectedEdgesList(2,:);
    sg = im0.edgeSegN(e)';
    types = im0.segVesType(sg);
    flowV = im0.edgeMeanFlow(e);
    im0.selectedEdgesList = [sg; e; flowV; types];
end
updateGui_Listbox(handles)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Executes on slider movement.
function slider_Callback(hObject, ~, handles) 
% hObject    handle to slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global im0
sliderval = round(get(hObject,'Value'));
im0.gui.frameNum = sliderval + ceil(im0.gui.numOfFrames/2) - 1;

update_TPMA_image(handles);
update_OCTV_image(handles);
overlayGraph(handles)
if all(im0.gui.positionFlag)
    positionFlag_OCTFlow(handles, im0.gui.mntrPt)
    positionFlag_TPAng(handles, im0.gui.mntrPt)
end

set(handles.frameNum, 'string', num2str(im0.gui.frameNum))

switch im0.gui.viewPlane
    case 'XY'
        im0.gui.mntrPt(3) = im0.gui.frameNum;
    case 'YZ'
        im0.gui.mntrPt(1) = im0.gui.frameNum;
    case 'XZ'
        im0.gui.mntrPt(2) = im0.gui.frameNum;
end

cylMode = get(handles.cylMode, 'value');
if cylMode
    update_Monitors(handles, im0.gui.edgeID, false)
else
    update_Monitors(handles, im0.gui.edgeID, true)
end
end

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider_CreateFcn(hObject, ~, ~)
% hObject    handle to slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
end



function alpha_Callback(hObject, ~, handles)
% hObject    handle to alpha (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global im0
im0.gui.alpha = str2double(get(hObject,'String'));
edgeID = str2double(get(handles.edgeID, 'string'));
update_Monitors(handles, edgeID, true)
end
% Hints: get(hObject,'String') returns contents of alpha as text
%        str2double(get(hObject,'String')) returns contents of alpha as a double


% --- Executes during object creation, after setting all properties.
function alpha_CreateFcn(hObject, ~, ~) 
% hObject    handle to alpha (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


% --- Executes on selection change in selectedEdgesList.
function selectedEdgesList_Callback(hObject,~, handles)  
global im0
if isempty(im0.selectedEdgesList)
    im0.gui.edgeID = [];
    im0.gui.segmentID = [];
    updateGui_edgeInfo(handles)
    computeEdgeFlow(handles, []);
    highlightEdges(handles, [])
    return
end

indx = get(hObject,'Value');
im0.gui.edgeID = im0.selectedEdgesList(2, indx);
im0.gui.segmentID = im0.selectedEdgesList(1, indx);
updateGui_edgeInfo(handles)
computeEdgeFlow(handles, im0.gui.edgeID);

pt = im0.nodePos(im0.nodeEdges(im0.gui.edgeID, :), :);
im0.gui.mntrPt = round(sum(pt, 1) / 2);
switch im0.gui.viewPlane
    case 'XY'
        im0.gui.frameNum = im0.gui.mntrPt(3);
    case 'YZ'
        im0.gui.frameNum = im0.gui.mntrPt(1);
    case 'XZ'
        im0.gui.frameNum = im0.gui.mntrPt(2);
end

set(handles.frameNum, 'string', num2str(im0.gui.frameNum))
sliderVal = im0.gui.frameNum - ceil(im0.gui.numOfFrames/2) + 1;
if sliderVal < 1; sliderVal = 1; end
set(handles.slider,'value', sliderVal)

im0.gui.positionFlag(2) = true;
highlightEdges(handles, im0.gui.edgeID)
end

% hObject    handle to selectedEdgesList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns selectedEdgesList contents as cell array
%        contents{get(hObject,'Value')} returns selected item from selectedEdgesList

function positionFlag_TPAng(handles, midP)
global im0
sz = size(im0.I);
c0 = midP - (im0.gui.winSize/2);    c1 = midP + (im0.gui.winSize/2);

switch im0.gui.viewPlane
    case 'XY'
        X = [c0(1), c1(1), 1, sz(1); c0(1), c1(1), sz(1), 1];
        Y = [1, sz(2), c0(2), c1(2); sz(2), 1, c0(2), c1(2)];
        line(handles.axes_TP_Ang, X, Y, 'color', 'g')
    case 'YZ'
        X = [(im0.szOCT(1) + 1)-c0(2), (im0.szOCT(1) + 1)-c1(2), 1, sz(2);...
            (im0.szOCT(1) + 1)-c0(2), (im0.szOCT(1) + 1)-c1(2), sz(2), 1];
        Y = [1, sz(3), c0(3), c1(3); sz(3), 1, c0(3), c1(3)];
        line(handles.axes_TP_Ang, X, Y, 'color', 'g')
    case 'XZ'
        X = [(im0.szOCT(1) + 1)-c0(1), (im0.szOCT(1) + 1)-c1(1), 1, sz(1);...
            (im0.szOCT(1) + 1)-c0(1), (im0.szOCT(1) + 1)-c1(1), sz(1), 1];
        Y = [1, sz(3), c0(3), c1(3); sz(3), 1, c0(3), c1(3)];
        line(handles.axes_TP_Ang, X, Y, 'color', 'g')
end
end


function positionFlag_OCTFlow(handles, midP)
global im0
sz = size(im0.I);
c0 = midP - (im0.gui.winSize/2);    c1 = midP + (im0.gui.winSize/2);


switch im0.gui.viewPlane
    case 'XY'
        X = [c0(1), c1(1), 1, sz(1); c0(1), c1(1), sz(1), 1];
        Y = [1, sz(2), c0(2), c1(2); sz(2), 1, c0(2), c1(2)];
        line(handles.axes_OCT_Flow, X, Y, 'color', 'g')
        
    case 'YZ'
        X = [(im0.szOCT(1) + 1)-c0(2), (im0.szOCT(1) + 1)-c1(2), 1, sz(2);...
            (im0.szOCT(1) + 1)-c0(2), (im0.szOCT(1) + 1)-c1(2), sz(2), 1];
        Y = [1, sz(3), c0(3), c1(3); sz(3), 1, c0(3), c1(3)];
        line(handles.axes_OCT_Flow, X, Y, 'color', 'g')
        
    case 'XZ'
        X = [(im0.szOCT(1) + 1)-c0(1), (im0.szOCT(1) + 1)-c1(1), 1, sz(1);...
            (im0.szOCT(1) + 1)-c0(1), (im0.szOCT(1) + 1)-c1(1), sz(1), 1];
        Y = [1, sz(3), c0(3), c1(3); sz(3), 1, c0(3), c1(3)];
        line(handles.axes_OCT_Flow, X, Y, 'color', 'g')
        
end
end



% --- Executes during object creation, after setting all properties.
function selectedEdgesList_CreateFcn(hObject, ~, ~)  
% hObject    handle to selectedEdgesList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end



function TPthreshold_Callback(hObject, ~, handles)  
% hObject    handle to TPthreshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global im0
im0.gui.TPthreshold = str2double(get(hObject,'String'));
update_TPMA_image(handles)
overlayGraph(handles)
end
% Hints: get(hObject,'String') returns contents of TPthreshold as text
%        str2double(get(hObject,'String')) returns contents of TPthreshold as a double


% --- Executes during object creation, after setting all properties.
function TPthreshold_CreateFcn(hObject, ~, ~)  
% hObject    handle to TPthreshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end



function colorbarUpperLim_Callback(hObject, ~, handles)  
% hObject    handle to colorbarUpperLim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global im0
l = get(handles.colorbarCheck,'Value');
upper = str2double(get(hObject,'String'));
if l
    set(handles.colorbalLowerLim,'string',num2str(-upper));
end
lower = str2double(get(handles.colorbalLowerLim,'String'));

if (upper == 0 && lower == 0) || (lower >= upper)
    set(handles.colorbarUpperLim,'string',num2str(im0.gui.colorbarUpperLim));
    set(handles.colorbalLowerLim,'string',num2str(im0.gui.colorbalLowerLim));
    errordlg('Limits are not valid','Input Error');
    return
end

im0.gui.colorbarUpperLim = upper;
im0.gui.colorbalLowerLim = lower;

caxis(handles.axes_Colorbar, [lower upper])
caxis(handles.axes_OCT_Flow, [lower upper])
caxis(handles.axes_OCT_Win_01, [lower upper])
caxis(handles.axes_OCT_Win_02, [lower upper])
caxis(handles.axes_OCT_Win_03, [lower upper])
end

% Hints: get(hObject,'String') returns contents of colorbarUpperLim as text
%        str2double(get(hObject,'String')) returns contents of colorbarUpperLim as a double


% --- Executes during object creation, after setting all properties.
function colorbarUpperLim_CreateFcn(hObject, ~, ~)  
% hObject    handle to colorbarUpperLim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end



function colorbalLowerLim_Callback(hObject, ~, handles)  
% hObject    handle to colorbalLowerLim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global im0
l = get(handles.colorbarCheck,'Value');
lower = str2double(get(hObject,'String'));
if l
    upper = -lower;
    set(handles.colorbarUpperLim,'string',num2str(upper));
end
upper = str2double(get(handles.colorbarUpperLim, 'string'));

if (upper == 0 && lower == 0) || (lower >= upper)
    set(handles.colorbarUpperLim,'string',num2str(im0.gui.colorbarUpperLim));
    set(handles.colorbalLowerLim,'string',num2str(im0.gui.colorbalLowerLim));
    errordlg('Limits are not valid','Input Error');
    return
end


im0.gui.colorbarUpperLim = upper;
im0.gui.colorbalLowerLim = lower;


caxis(handles.axes_Colorbar, [lower upper])
caxis(handles.axes_OCT_Flow, [lower upper])
caxis(handles.axes_OCT_Win_01, [lower upper])
caxis(handles.axes_OCT_Win_02, [lower upper])
caxis(handles.axes_OCT_Win_03, [lower upper])
end
% Hints: get(hObject,'String') returns contents of colorbalLowerLim as text
%        str2double(get(hObject,'String')) returns contents of colorbalLowerLim as a double


% --- Executes during object creation, after setting all properties.
function colorbalLowerLim_CreateFcn(hObject, ~, ~)  
% hObject    handle to colorbalLowerLim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


% --- Executes on button press in colorbarCheck.
function colorbarCheck_Callback(hObject, ~, ~)  
global im0
im0.gui.colorbarCheck = get(hObject,'Value');
end
% hObject    handle to colorbarCheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of colorbarCheck



function cylDiam_Callback(hObject, ~, handles)  
% hObject    handle to cylDiam (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global im0
input = strsplit(get(hObject,'String'));
if length(input) > 1
    val = str2double(cell2mat(input(1)));
else
    val = str2double(get(hObject,'String'));
end

im0.cylDia(im0.gui.edgeID) = val;

set(hObject, 'String', [num2str(val), ' μm'])
computeEdgeFlow(handles, im0.gui.edgeID)
update_Monitors(handles, im0.gui.edgeID, true)
end
% Hints: get(hObject,'String') returns contents of cylDiam as text
%        str2double(get(hObject,'String')) returns contents of cylDiam as a double


% --- Executes during object creation, after setting all properties.
function cylDiam_CreateFcn(hObject, ~, ~)  
% hObject    handle to cylDiam (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end



function cylLength_Callback(hObject, ~, ~)
l = get(hObject,'String');
set(hObject, 'String', [l, ' μm'])
end
% hObject    handle to cylLength (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of cylLength as text
%        str2double(get(hObject,'String')) returns contents of cylLength as a double


% --- Executes during object creation, after setting all properties.
function cylLength_CreateFcn(hObject, ~, ~)  
% hObject    handle to cylLength (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end



% --- Executes during object creation, after setting all properties.
function edgeDiam_Est_CreateFcn(hObject, ~, ~)  
% hObject    handle to edgeDiam_Est (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


% --- Executes on button press in win03Up.
function win03Up_Callback(~, ~, handles)  
global im0
im0.cylNodes(im0.gui.edgeID,[3,6]) = im0.cylNodes(im0.gui.edgeID,[3,6]) - im0.gui.step;
update_Monitors(handles, im0.gui.edgeID, true)
computeEdgeFlow(handles, im0.gui.edgeID)

end
% hObject    handle to win03Up (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in win03Down.
function win03Down_Callback(~, ~, handles)  
global im0
im0.cylNodes(im0.gui.edgeID,[3,6]) = im0.cylNodes(im0.gui.edgeID,[3,6]) + im0.gui.step;
update_Monitors(handles, im0.gui.edgeID, true)
computeEdgeFlow(handles, im0.gui.edgeID)
end
% hObject    handle to win03Down (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in win03Left.
function win03Left_Callback(~, ~, handles)  
global im0
im0.cylNodes(im0.gui.edgeID,[1,4]) = im0.cylNodes(im0.gui.edgeID,[1,4]) + im0.gui.step;
update_Monitors(handles, im0.gui.edgeID, true)
computeEdgeFlow(handles, im0.gui.edgeID)
end
% hObject    handle to win03Left (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in win03Right.
function win03Right_Callback(~, ~, handles)  
global im0
im0.cylNodes(im0.gui.edgeID,[1,4]) = im0.cylNodes(im0.gui.edgeID,[1,4]) - im0.gui.step;
update_Monitors(handles, im0.gui.edgeID, true)
computeEdgeFlow(handles, im0.gui.edgeID)
end
% hObject    handle to win03Right (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in win03RotAntiClock.
function win03RotAntiClock_Callback(~, ~, handles)  
global im0
p0 = im0.cylNodes(im0.gui.edgeID, 1:3);
p1 = im0.cylNodes(im0.gui.edgeID, 4:6);
pm = (p0 + p1)/2;

theta = -im0.gui.rotAngl * pi/180;
RM = [cos(theta), 0, sin(theta); 0, 1, 0; -sin(theta), 0, cos(theta)];

p0 = (RM * (p0 - pm)')' + pm;
p1 = (RM * (p1 - pm)')' + pm;

im0.cylNodes(im0.gui.edgeID, :) = [p0, p1];

update_Monitors(handles, im0.gui.edgeID, true)
computeEdgeFlow(handles, im0.gui.edgeID)
end
% hObject    handle to win03RotAntiClock (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in win03RotClock.
function win03RotClock_Callback(~, ~, handles)  
global im0
p0 = im0.cylNodes(im0.gui.edgeID, 1:3);
p1 = im0.cylNodes(im0.gui.edgeID, 4:6);
pm = (p0 + p1)/2;

theta = im0.gui.rotAngl * pi/180;
RM = [cos(theta), 0, sin(theta); 0, 1, 0; -sin(theta), 0, cos(theta)];

p0 = (RM * (p0 - pm)')' + pm;
p1 = (RM * (p1 - pm)')' + pm;

im0.cylNodes(im0.gui.edgeID, :) = [p0, p1];

update_Monitors(handles, im0.gui.edgeID, true)
computeEdgeFlow(handles, im0.gui.edgeID)
end
% hObject    handle to win03RotClock (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in win02Up.
function win02Up_Callback(~, ~, handles)  
global im0
im0.cylNodes(im0.gui.edgeID,[3,6]) = im0.cylNodes(im0.gui.edgeID,[3,6]) - im0.gui.step;
update_Monitors(handles, im0.gui.edgeID, true)
computeEdgeFlow(handles, im0.gui.edgeID)
end
% hObject    handle to win02Up (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in win02Down.
function win02Down_Callback(~, ~, handles)  
global im0
im0.cylNodes(im0.gui.edgeID,[3,6]) = im0.cylNodes(im0.gui.edgeID,[3,6]) + im0.gui.step;
update_Monitors(handles, im0.gui.edgeID, true)
computeEdgeFlow(handles, im0.gui.edgeID)
end
% hObject    handle to win02Down (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in win02Left.
function win02Left_Callback(~, ~, handles)  
global im0
im0.cylNodes(im0.gui.edgeID,[2,5]) = im0.cylNodes(im0.gui.edgeID,[2,5]) + im0.gui.step;
update_Monitors(handles, im0.gui.edgeID, true)
computeEdgeFlow(handles, im0.gui.edgeID)
end
% hObject    handle to win02Left (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in win02Right.
function win02Right_Callback(~, ~, handles)  
global im0
im0.cylNodes(im0.gui.edgeID,[2,5]) = im0.cylNodes(im0.gui.edgeID,[2,5]) - im0.gui.step;
update_Monitors(handles, im0.gui.edgeID, true)
computeEdgeFlow(handles, im0.gui.edgeID)
end
% hObject    handle to win02Right (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in win02RotAntiClock.
function win02RotAntiClock_Callback(~, ~, handles)  
global im0
p0 = im0.cylNodes(im0.gui.edgeID, 1:3);
p1 = im0.cylNodes(im0.gui.edgeID, 4:6);
pm = (p0 + p1)/2;

theta = im0.gui.rotAngl * pi/180;
RM = [1, 0, 0; 0, cos(theta), -sin(theta); 0, sin(theta), cos(theta)];

p0 = (RM * (p0 - pm)')' + pm;
p1 = (RM * (p1 - pm)')' + pm;

im0.cylNodes(im0.gui.edgeID, :) = [p0, p1];

update_Monitors(handles, im0.gui.edgeID, true)
computeEdgeFlow(handles, im0.gui.edgeID)
end
% hObject    handle to win02RotAntiClock (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in win02RotClock.
function win02RotClock_Callback(~, ~, handles)  
global im0
p0 = im0.cylNodes(im0.gui.edgeID, 1:3);
p1 = im0.cylNodes(im0.gui.edgeID, 4:6);
pm = (p0 + p1)/2;

theta = -im0.gui.rotAngl * pi/180;
RM = [1, 0, 0; 0, cos(theta), -sin(theta); 0, sin(theta), cos(theta)];

p0 = (RM * (p0 - pm)')' + pm;
p1 = (RM * (p1 - pm)')' + pm;

im0.cylNodes(im0.gui.edgeID, :) = [p0, p1];

update_Monitors(handles, im0.gui.edgeID, true)
computeEdgeFlow(handles, im0.gui.edgeID)
end
% hObject    handle to win02RotClock (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in win01Up.
function win01Up_Callback(~, ~, handles)  
global im0
im0.cylNodes(im0.gui.edgeID,[2,5]) = im0.cylNodes(im0.gui.edgeID,[2,5]) - im0.gui.step;
update_Monitors(handles, im0.gui.edgeID, true)
computeEdgeFlow(handles, im0.gui.edgeID)
end
% hObject    handle to win01Up (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in win01Down.
function win01Down_Callback(~, ~, handles)  
global im0
im0.cylNodes(im0.gui.edgeID,[2,5]) = im0.cylNodes(im0.gui.edgeID,[2,5]) + im0.gui.step;
update_Monitors(handles, im0.gui.edgeID, true)
computeEdgeFlow(handles, im0.gui.edgeID)
end
% hObject    handle to win01Down (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in win01Left.
function win01Left_Callback(~, ~, handles)  
global im0
im0.cylNodes(im0.gui.edgeID,[1,4]) = im0.cylNodes(im0.gui.edgeID,[1,4]) - im0.gui.step;
update_Monitors(handles, im0.gui.edgeID, true)
computeEdgeFlow(handles, im0.gui.edgeID)
end
% hObject    handle to win01Left (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in win01Right.
function win01Right_Callback(~, ~, handles)  
global im0
im0.cylNodes(im0.gui.edgeID,[1,4]) = im0.cylNodes(im0.gui.edgeID,[1,4]) + im0.gui.step;
update_Monitors(handles, im0.gui.edgeID, true)
computeEdgeFlow(handles, im0.gui.edgeID)
end
% hObject    handle to win01Right (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in win01RotAntiClock.
function win01RotAntiClock_Callback(~, ~, handles)  
global im0
p0 = im0.cylNodes(im0.gui.edgeID, 1:3);
p1 = im0.cylNodes(im0.gui.edgeID, 4:6);
pm = (p0 + p1)/2;

theta = -im0.gui.rotAngl * pi/180;
RM = [cos(theta), -sin(theta), 0; sin(theta), cos(theta), 0; 0, 0, 1];

p0 = (RM * (p0 - pm)')' + pm;
p1 = (RM * (p1 - pm)')' + pm;

im0.cylNodes(im0.gui.edgeID, :) = [p0, p1];

update_Monitors(handles, im0.gui.edgeID, true)
computeEdgeFlow(handles, im0.gui.edgeID)
end
% hObject    handle to win01RotAntiClock (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in win01RotClock.
function win01RotClock_Callback(~, ~, handles)  
global im0
p0 = im0.cylNodes(im0.gui.edgeID, 1:3);
p1 = im0.cylNodes(im0.gui.edgeID, 4:6);
pm = (p0 + p1)/2;

theta = im0.gui.rotAngl * pi/180;
RM = [cos(theta), -sin(theta), 0; sin(theta), cos(theta), 0; 0, 0, 1];

p0 = (RM * (p0 - pm)')' + pm;
p1 = (RM * (p1 - pm)')' + pm;

im0.cylNodes(im0.gui.edgeID, :) = [p0, p1];

update_Monitors(handles, im0.gui.edgeID, true)
computeEdgeFlow(handles, im0.gui.edgeID)
end
% hObject    handle to win01RotClock (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in addEdge.
function addEdge_Callback(~, ~, handles)  
global im0
if isempty(im0.selectedEdgesList) ||isempty(intersect(im0.gui.edgeID, im0.selectedEdgesList(2,:)))
    type = im0.segVesType(im0.gui.segmentID);
    im0.selectedEdgesList = [im0.selectedEdgesList, [im0.gui.segmentID; im0.gui.edgeID; im0.gui.edgeFlowV; type]];
    
    fv = round(im0.selectedEdgesList(3,:), 4);
    list = [im0.selectedEdgesList(1:2,:); fv];
    set(handles.selectedEdgesList,'string',num2str(list'));
    %set(handles.selectedEdgesList,'string',num2str(im0.selectedEdgesList'));
    im0.edgeFlag(im0.gui.edgeID) = 1;
end
[~, l] = size(im0.selectedEdgesList);
set(handles.selectedEdgesList, 'Value', l);
end
% hObject    handle to addEdge (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function moveListSelection(handles, ask, dltCrnt)

if dltCrnt && ask
    answer = questdlg('Are you sure you want to remove edge from list?', 'Remove Edge');
    
    switch answer
        case 'No'
            return
        case 'Cancel'
            return
    end
end

global im0

indx = get(handles.selectedEdgesList, 'Value');

if dltCrnt; im0.selectedEdgesList(:, indx) = []; end

switch handles.cntxtMenu_lst_nex.Checked
    case 'on'
        if ~dltCrnt; indx = indx + 1; end
        if indx > size(im0.selectedEdgesList, 2)
            indx = size(im0.selectedEdgesList, 2);
        end
    case 'off'
        indx = indx - 1;
         if indx == 0; indx = 1; end
end

set(handles.selectedEdgesList, 'Value', indx);

if isempty(im0.selectedEdgesList)
    im0.gui.edgeID = [];
    im0.gui.segmentID = [];
    updateGui_edgeInfo(handles)
    updateGui_Listbox(handles)
    highlightEdges(handles, [])
    return
end
im0.gui.edgeID = im0.selectedEdgesList(2, indx);
im0.gui.segmentID = im0.selectedEdgesList(1, indx);

computeEdgeFlow(handles, im0.gui.edgeID)
updateGui_edgeInfo(handles)

pt = im0.nodePos(im0.nodeEdges(im0.gui.edgeID, :), :);
im0.gui.mntrPt = round(sum(pt, 1) / 2);
switch im0.gui.viewPlane
    case 'XY'
        im0.gui.frameNum = im0.gui.mntrPt(3);
    case 'YZ'
        im0.gui.frameNum = im0.gui.mntrPt(1);
    case 'XZ'
        im0.gui.frameNum = im0.gui.mntrPt(2);
end

set(handles.frameNum, 'string', num2str(im0.gui.frameNum))
sliderVal = im0.gui.frameNum - ceil(im0.gui.numOfFrames/2) + 1;
set(handles.slider,'value', sliderVal)

highlightEdges(handles, im0.gui.edgeID)
end

% --- Executes on button press in rejectEdge.
function rejectEdge_Callback(~, ~, handles)  
%%
global im0

if ~isempty(intersect(im0.selectedEdgesList(2,:), im0.gui.edgeID))
    moveListSelection(handles, true, true)
end

im0.edgeFlag(im0.gui.edgeID) = -1;

end
% hObject    handle to rejectEdge (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in rstCylPos.
function rstCylPos_Callback(~, ~, handles)  
global im0
p0 = im0.nodeEdges(im0.gui.edgeID, 1);  p0 = im0.nodePos(p0,:);
p1 = im0.nodeEdges(im0.gui.edgeID, 2);  p1 = im0.nodePos(p1,:);
im0.cylNodes(im0.gui.edgeID,:) = [p0 p1];
update_Monitors(handles, im0.gui.edgeID, true)
computeEdgeFlow(handles, im0.gui.edgeID)
end
% hObject    handle to rstCylPos (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function step_Callback(hObject, ~, ~)  
global im0
s = get(hObject,'String');
im0.gui.step = str2double(s);
set(hObject, 'String', [s, ' pixel'])
end
% hObject    handle to step (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of step as text
%        str2double(get(hObject,'String')) returns contents of step as a double


% --- Executes during object creation, after setting all properties.
function step_CreateFcn(hObject, ~, ~)  
% hObject    handle to step (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end



function winSize_Callback(hObject, ~, handles)  
% hObject    handle to winSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global im0
sz = get(hObject,'String');
if str2double(sz) < 2; sz = '2'; end
im0.gui.winSize = str2double(sz);
set(hObject, 'String', [sz, ' pixel'])
edge = str2double(get(handles.edgeID, 'string'));
update_Monitors(handles, edge, true)
end
% Hints: get(hObject,'String') returns contents of winSize as text
%        str2double(get(hObject,'String')) returns contents of winSize as a double


% --- Executes during object creation, after setting all properties.
function winSize_CreateFcn(hObject, ~, ~)  
% hObject    handle to winSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end



function edgeLinewidth_monitors_Callback(hObject, ~, handles)  
% hObject    handle to edgeLinewidth_monitors (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global im0
w = get(hObject,'String');
im0.gui.edgeLinewidth_monitors = str2double(w);
set(hObject, 'String', [w, ' pt'])
update_Monitors(handles, im0.gui.edgeID, true)
end
% Hints: get(hObject,'String') returns contents of edgeLinewidth_monitors as text
%        str2double(get(hObject,'String')) returns contents of edgeLinewidth_monitors as a double


% --- Executes during object creation, after setting all properties.
function edgeLinewidth_monitors_CreateFcn(hObject, ~, ~)  
% hObject    handle to edgeLinewidth_monitors (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end



function nodePointWidth_monitors_Callback(hObject, ~, handles)  
% hObject    handle to nodePointWidth_monitors (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global im0
w = get(hObject,'String');
im0.gui.nodePointWidth_monitors = str2double(w);
set(hObject, 'String', [w, ' pt'])
update_Monitors(handles, im0.gui.edgeID, true)
end
% Hints: get(hObject,'String') returns contents of nodePointWidth_monitors as text
%        str2double(get(hObject,'String')) returns contents of nodePointWidth_monitors as a double


% --- Executes during object creation, after setting all properties.
function nodePointWidth_monitors_CreateFcn(hObject, ~, ~)  
% hObject    handle to nodePointWidth_monitors (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end



function numOfFrames_Callback(hObject, ~, handles)  
% hObject    handle to numOfFrames (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global im0
global V
im0.gui.numOfFrames = str2double(get(hObject,'String'));

    %%%%
hwait = waitbar(0.1,'Processing XY Projection (TP Angiogram)');
V.TP_XY = intensProj(im0.I, im0.gui.numOfFrames, 'XY', im0.gui.projType);
waitbar(1/6,hwait,'Processing YZ Projection (TP Angiogram)');
V.TP_YZ = intensProj(im0.I, im0.gui.numOfFrames, 'YZ', im0.gui.projType);
waitbar(2/6,hwait,'Processing XZ Projection (TP Angiogram)');
V.TP_XZ = intensProj(im0.I, im0.gui.numOfFrames, 'XZ', im0.gui.projType);
%%%%%%%%%%
waitbar(3/6,hwait,'Processing XY Projection (Flow Velocity)');
V.OCT_Flow_XY = intensProj(im0.OCT_Flow, im0.gui.numOfFrames, 'XY', im0.gui.projType);
waitbar(4/6,hwait,'Processing YZ Projection (Flow Velocity)');
V.OCT_Flow_YZ = intensProj(im0.OCT_Flow, im0.gui.numOfFrames, 'YZ', im0.gui.projType);
waitbar(5/6,hwait,'Processing XZ Projection (Flow Velocity)');
V.OCT_Flow_XZ = intensProj(im0.OCT_Flow, im0.gui.numOfFrames, 'XZ', im0.gui.projType);
%%%%%%%%%%
waitbar(0.98,hwait,'Fininshing up');
pause(0.25)
waitbar(1,hwait,'Done!');
pause(0.25)
close(hwait)
    %%%%

switch im0.gui.viewPlane
    case 'XY'
        sliderLim = max(size(V.TP_XY, 3), size(V.OCT_Flow_XY, 3));
    case 'YZ'
        sliderLim = max(size(V.TP_YZ, 3), size(V.OCT_Flow_YZ, 3));
    case 'XZ'
        sliderLim = max(size(V.TP_XZ, 3), size(V.OCT_Flow_XZ, 3));
end
sliderVal = round(get(handles.slider,'Value'));
if sliderVal > sliderLim; sliderVal = sliderLim; end
set(handles.slider,'Value', sliderVal)
set(handles.slider,'max', sliderLim)

update_TPMA_image(handles)
update_OCTV_image(handles)
overlayGraph(handles)
end

% Hints: get(hObject,'String') returns contents of numOfFrames as text
%        str2double(get(hObject,'String')) returns contents of numOfFrames as a double


% --- Executes during object creation, after setting all properties.
function numOfFrames_CreateFcn(hObject, ~, ~)  
% hObject    handle to numOfFrames (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end



function frameNum_Callback(hObject, ~, handles)  
global im0 V
val = round(str2double(get(hObject,'String')));
val(val < 1) = 1;
switch im0.gui.viewPlane
    case 'XY'
        val(val > im0.gui.maxVolumeSz(3)) = im0.gui.maxVolumeSz(3);
        frm = val - floor((im0.gui.numOfFrames - 1) / 2);
        frm(frm < ceil(im0.gui.numOfFrames/2)) = ceil(im0.gui.numOfFrames/2);
        frm(frm > V.prjVolumeSz(3) && frm <= im0.gui.maxVolumeSz(3)) = V.prjVolumeSz(3);
    case 'YZ'
        val(val > im0.gui.maxVolumeSz(1)) = im0.gui.maxVolumeSz(1);
        frm = val - floor((im0.gui.numOfFrames - 1) / 2);
        frm(frm < ceil(im0.gui.numOfFrames/2)) = ceil(im0.gui.numOfFrames/2);
        frm(frm > V.prjVolumeSz(1) && frm <= im0.gui.maxVolumeSz(1)) = V.prjVolumeSz(1);
    case 'XZ'
        val(val > im0.gui.maxVolumeSz(2)) = im0.gui.maxVolumeSz(2);
        frm = val - floor((im0.gui.numOfFrames - 1) / 2);
        frm(frm < ceil(im0.gui.numOfFrames/2)) = ceil(im0.gui.numOfFrames/2);
        frm(frm > V.prjVolumeSz(2) && frm <= im0.gui.maxVolumeSz(2)) = V.prjVolumeSz(2);
end
set(hObject, 'String', num2str(val))

im0.gui.frameNum = val;
set(handles.slider,'value', frm)

update_TPMA_image(handles)
update_OCTV_image(handles);
overlayGraph(handles)
end
% Hints: get(hObject,'String') returns contents of frameNum as text
%        str2double(get(hObject,'String')) returns contents of frameNum as a double


% --- Executes during object creation, after setting all properties.
function frameNum_CreateFcn(hObject, ~, ~)  
% hObject    handle to frameNum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


% --- Executes on selection change in projType.
function projType_Callback(hObject, eventdata, handles)  
global im0
contents = cellstr(get(hObject,'String'));
im0.gui.projType = contents{get(hObject,'Value')};
numOfFrames_Callback( (handles.numOfFrames), eventdata, handles )
end
% hObject    handle to projType (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns projType contents as cell array
%        contents{get(hObject,'Value')} returns selected item from projType


% --- Executes during object creation, after setting all properties.
function projType_CreateFcn(hObject, ~, ~)  
% hObject    handle to projType (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


% --- Executes on selection change in viewPlane.
function viewPlane_Callback(hObject, ~, handles)  
global im0 V
% hObject    handle to viewPlane (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
contents = cellstr(get(hObject,'String'));
newViewPlane = contents{get(hObject,'Value')};
sz = max (size(im0.I), size(im0.OCT_Flow));

xLimit = xlim(handles.axes_OCT_Flow);
yLimit = ylim(handles.axes_OCT_Flow);

p1 = ceil((sum(xLimit)-1)/2);
p2 = ceil((sum(yLimit)-1)/2);

switch im0.gui.viewPlane
    case 'XY'
        midP = [p1, p2, im0.gui.frameNum];
        
    case 'YZ'
        midP = [im0.gui.frameNum, p1, p2];
        
    case 'XZ'
        midP = [p1, im0.gui.frameNum, p2];
end


switch newViewPlane
    case 'XY'
        sliderLim = max(size(V.TP_XY, 3), size(V.OCT_Flow_XY, 3));
        set(handles.slider,'max',sliderLim)
        set(handles.slider,'value', (midP(3) - ceil(im0.gui.numOfFrames/2) + 1))
        im0.gui.frameNum = midP(3);
        im0.gui.mntrPt(3) = im0.gui.frameNum;
        
        xlim(handles.axes_TP_Ang, [1, sz(1)])
        ylim(handles.axes_TP_Ang, [1, sz(2)])
        xlim(handles.axes_OCT_Flow, [1, sz(1)])
        ylim(handles.axes_OCT_Flow, [1, sz(2)])
        
    case 'YZ'
        sliderLim = max(size(V.TP_YZ, 3), size(V.OCT_Flow_YZ, 3));
        set(handles.slider,'max',sliderLim)
        set(handles.slider,'value', (midP(1) - ceil(im0.gui.numOfFrames/2) + 1))
        im0.gui.frameNum = midP(1);
        im0.gui.mntrPt(1) = im0.gui.frameNum;
        
        xlim(handles.axes_TP_Ang, [1, sz(2)])
        ylim(handles.axes_TP_Ang, [1, sz(3)])
        xlim(handles.axes_OCT_Flow, [1, sz(2)])
        ylim(handles.axes_OCT_Flow, [1, sz(3)])
        
    case 'XZ'
        sliderLim = max(size(V.TP_XZ, 3), size(V.OCT_Flow_XZ, 3));
        set(handles.slider,'max',sliderLim)
        set(handles.slider,'value', (midP(2) - ceil(im0.gui.numOfFrames/2) + 1))
        im0.gui.frameNum = midP(2);
        im0.gui.mntrPt(2) = im0.gui.frameNum;
        
        xlim(handles.axes_TP_Ang, [1, sz(1)])
        ylim(handles.axes_TP_Ang, [1, sz(3)])
        xlim(handles.axes_OCT_Flow, [1, sz(1)])
        ylim(handles.axes_OCT_Flow, [1, sz(3)])   
end

im0.gui.viewPlane = newViewPlane;
V.vctrs = vInfo(im0.nodeEdges, im0.nodePos, im0.gui.viewPlane);
update_OCTV_image(handles);
update_TPMA_image(handles);
overlayGraph(handles)

if all(im0.gui.positionFlag)
    positionFlag_OCTFlow(handles, im0.gui.mntrPt)
    positionFlag_TPAng(handles, im0.gui.mntrPt)
end

end



% Hints: contents = cellstr(get(hObject,'String')) returns viewPlane contents as cell array
%        contents{get(hObject,'Value')} returns selected item from viewPlane


% --- Executes during object creation, after setting all properties.
function viewPlane_CreateFcn(hObject, ~, ~)  
% hObject    handle to viewPlane (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


function edgeLinewidth_axes_Callback(hObject, ~, handles)  
% hObject    handle to edgeLinewidth_axes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global im0
w = get(hObject,'String');
if str2double(w) < 0.001
    w = '0.001';
end
im0.gui.edgeLinewidth_axes = str2double(w);
set(hObject, 'String', [w, ' pt'])
update_TPMA_image(handles)
update_OCTV_image(handles)
overlayGraph(handles)
end
% Hints: get(hObject,'String') returns contents of edgeLinewidth_axes as text
%        str2double(get(hObject,'String')) returns contents of edgeLinewidth_axes as a double


% --- Executes during object creation, after setting all properties.
function edgeLinewidth_axes_CreateFcn(hObject, ~, ~)  
% hObject    handle to edgeLinewidth_axes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


%%-------------------------------------------------------------------------
function nodePointWidth_axes_Callback(hObject, ~, handles)  
% hObject    handle to nodePointWidth_axes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global im0
w = get(hObject,'String');
if str2double(w) < 0.001
    w = '0.001';
end
im0.gui.nodePointWidth_axes = str2double(w);
set(hObject, 'String', [w, ' pt'])
update_TPMA_image(handles)
update_OCTV_image(handles)
overlayGraph(handles)
end
% Hints: get(hObject,'String') returns contents of nodePointWidth_axes as text
%        str2double(get(hObject,'String')) returns contents of nodePointWidth_axes as a double


% --- Executes during object creation, after setting all properties.
function nodePointWidth_axes_CreateFcn(hObject, ~, ~)  
% hObject    handle to nodePointWidth_axes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end



%%-------------------------------------------------------------------------
function edgeDiam_Est_Callback(hObject, ~, ~)  
% hObject    handle to edgeDiam_Est (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global im0
d = get(hObject,'String');
im0.edgeDiaEst(im0.gui.edgeID) = str2double(d);
set(hObject, 'String', [d, ' μm'])
end
% Hints: get(hObject,'String') returns contents of edgeDiam_Est as text
%        str2double(get(hObject,'String')) returns contents of edgeDiam_Est as a double


% --- Executes on button press in radiobutton_XY.
function radiobutton_XY_Callback(~, ~, handles)  
% hObject    handle to radiobutton_XY (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global im0
im0.gui.angWind = 'XY';
update_Monitors(handles, im0.gui.edgeID, true)
end
% Hint: get(hObject,'Value') returns toggle state of radiobutton_XY


% --- Executes on button press in radiobutton_YZ.
function radiobutton_YZ_Callback(~, ~, handles)  
% hObject    handle to radiobutton_YZ (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global im0
im0.gui.angWind = 'YZ';
update_Monitors(handles, im0.gui.edgeID, true)
end
% Hint: get(hObject,'Value') returns toggle state of radiobutton_YZ


% --- Executes on button press in radiobutton_XZ.
function radiobutton_XZ_Callback(~, ~, handles)  
% hObject    handle to radiobutton_XZ (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global im0
im0.gui.angWind = 'XZ';
update_Monitors(handles, im0.gui.edgeID, true)
end
% Hint: get(hObject,'Value') returns toggle state of radiobutton_XZ


% --- Executes on button press in inspectMode.
function inspectMode_Callback(~, ~, handles)
global im0
handles.cntxtMenu_inspect.Checked = 'on';
handles.cntxtMenu_cylinder.Checked = 'off';
handles.cntxtMenu_mntrsOff.Checked = 'off';
update_Monitors(handles, im0.gui.edgeID, true)
end
% hObject    handle to inspectMode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of inspectMode


% --------------------------------------------------------------------
function Untitled_1_Callback(hObject, eventdata, handles)
end
% hObject    handle to Untitled_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function exprtEdgeFlw_Callback(hObject, eventdata, handles)  
global im0

if ~isfield(im0, 'T')
    prcssEdgFlw_Callback(hObject, eventdata, handles)
end
T = im0.T;
[file,path] = uiputfile('FlowData00.mat');
save([path file], 'T');
end
% hObject    handle to exprtEdgeFlw (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Untitled_3_Callback(hObject, eventdata, handles)
end
% hObject    handle to Untitled_3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function prcssSgmntFlwV_Callback(~, ~, ~)

global im0


% T = table(sgmnt,edgeNo,edgeFlwVel,cylDiam);
% rawData = table2array(T);

% SE = im0.selectedEdgesList(:, im0.selectedEdgesList(3,:) ~= 0);
lst = im0.selectedEdgesList;
im0.selectedEdgesList(:, im0.selectedEdgesList(3,:) == 0) = [];
sgmntIndx = unique(im0.selectedEdgesList(1,:));


l = length(sgmntIndx);

sgmnt = zeros(l,1);
maxFlowVel = zeros(l,1);
meanFlowVel = zeros(l,1);
STD = zeros(l,1);
STErr = zeros(l,1);
numOfEdges = zeros(l,1);
sgmntDiam = zeros(l,1);
r = 0;
hwait = waitbar(0.15,'Processing Segnmt Flow v');

for sg = sgmntIndx
    r = r+1;
    
    sgmntFlow = im0.selectedEdgesList(3, im0.selectedEdgesList(1,:) == sg);
    
    if isempty(sgmntFlow)
        continue
    end
    
    meanFlow = mean(sgmntFlow);
    if meanFlow > 0
        maxFlow = max(sgmntFlow);
    else
        maxFlow = min(sgmntFlow);
    end
    
    s = std(sgmntFlow);
    e = s/sqrt(length(sgmntFlow));
    
    
    sgmnt (r,1)= sg;
    maxFlowVel (r,1)= maxFlow;
    meanFlowVel (r,1)= meanFlow;
    STD (r,1)= s;
    STErr (r,1)= e;
    numOfEdges (r,1) = length(sgmntFlow);
    sgmntDiam(r,1) = im0.segDiam(sg);
    
    waitbar(r/l, hwait, 'Processing Segnmt Flow v');
end

waitbar(1, hwait, 'Segmnt Flow v Process Complete');
data = table(sgmnt, maxFlowVel, meanFlowVel, STD, STErr, numOfEdges, sgmntDiam);
im0.data = data;
im0.selectedEdgesList = lst;

pause(0.2)
close(hwait)
end



% hObject    handle to prcssSgmntFlwV (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --------------------------------------------------------------------
function prcssEdgFlw_Callback(~, ~, ~)
global im0
edgeNo = im0.selectedEdgesList(2,:)';
edgeFlwVel = im0.selectedEdgesList(3,:)';
cylDiam = im0.cylDia(edgeNo);
sgmnt = im0.selectedEdgesList(1,:)';
T = table(sgmnt,edgeNo,edgeFlwVel,cylDiam);
im0.T = T;
end
% hObject    handle to prcssEdgFlw (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function exprtSgmntFlw_Callback(hObject, eventdata, handles)  
global im0

if ~isfield(im0, 'data')
    prcssSgmntFlwV_Callback(hObject, eventdata, handles)
end
data = im0.data;

[file, path] = uiputfile('SgmntFlowData_Session00.mat');
save([path file], 'data');
% hObject    handle to exprtSgmntFlw (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
end

% --------------------------------------------------------------------
function flwV3DMtrx_Callback(hObject, eventdata, handles)  
global im0

if ~isfield(im0, 'data')
    prcssSgmntFlwV_Callback(hObject, eventdata, handles)
end
data = im0.data;

edgeDia = im0.edgeDiaEst;

edgeNodes = im0.nodeEdges;
nodePos = im0.nodePos;
sz = size(im0.I);
flowData = table2array(data);
edgeSgmnts = im0.edgeSegN;
Hvox = im0.Hvox;

sgmnts = flowData(:,1);
mat = zeros(sz);
[XX, YY, ZZ] = meshgrid(1:sz(1),1:sz(2),1:sz(3));
itr = 0;


hwait = waitbar(0,'Processing matrix');
l = length(sgmnts);
for i = sgmnts'
    itr = itr + 1;
    edges = find(edgeSgmnts == i);
    
    r = edgeDia(edges);
    r = median(r);
    flowV = flowData(flowData(:,1) == i, 3);
    
    for j = edges'
        nodes = edgeNodes(j,:);
        p0 = nodePos(nodes(1),:);
        p1 = nodePos(nodes(2),:);
        v = p1-p0;       %edge vector
        vScld = (v .* Hvox);  %scaled edge vector
        vScldMag = norm(vScld);
        
        X = (XX - p0(1))* Hvox(1);
        Y = (YY - p0(2))* Hvox(2);
        Z = (ZZ - p0(3))* Hvox(3);
        d = sqrt((Y*vScld(3) - Z*vScld(2)).^2 + (Z*vScld(1) -...
            X*vScld(3)).^2 + (X*vScld(2) - Y*vScld(1)).^2)./ vScldMag;
        v0Prj = (X*vScld(1) + Y*vScld(2) + Z*vScld(3))/vScldMag;
        mat(d<=r & v0Prj<=vScldMag & v0Prj>=0) = flowV;
    end
    
    waitbar(itr/l, hwait, 'Processing matrix');
end
waitbar(1, hwait, 'Done!');

pause(0.3)
close(hwait)

[file, path] = uiputfile('Flow3DMat00.mat');
save([path file], 'mat');

% hObject    handle to flwV3DMtrx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
end

% --------------------------------------------------------------------
function srtEdgeLst_Callback(~, ~, handles)  
global im0
updateGui_Listbox(handles)
if ~isempty(im0.selectedEdgesList)
    val = get(handles.selectedEdgesList, 'Value');
    e = im0.selectedEdgesList(2, val);
    
    lst = im0.selectedEdgesList';
    lst = sortrows(lst);
    im0.selectedEdgesList = lst';
    updateGui_Listbox(handles)
    
    nuVal = find(im0.selectedEdgesList(2,:) == e);
    set(handles.selectedEdgesList, 'Value', nuVal);
end
% hObject    handle to srtEdgeLst (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

end
% --------------------------------------------------------------------
function reCompFloV_Callback(hObject, eventdata, handles)  
global im0
method = choosedialog;
if method == 1
    cutoffMethodDialog %- this is the input gui
    uiwait(cutoffMethodDialog)
    if (im0.cuoffMethd.canceled); return; end
    recomputeFlow_edgeCutOff(hObject, eventdata, handles)
    updateGui_Listbox(handles)
elseif method == 2
    cutoffMethodDialog %- this is the input gui
    uiwait(cutoffMethodDialog)
    if (im0.cuoffMethd.canceled); return; end
    recomputeFlow_segmentCutOff
    updateGui_Listbox(handles)
end
end
% hObject    handle to reCompFloV (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


function choice = choosedialog

    options = {'cutoff -edge based', 'cutoff -segment based'};

    d = dialog('Position',[300 300 250 150],'Name','Select One');
    txt = uicontrol('Parent',d,'Style','text',...
           'Position',[20 80 210 45],...
           'String','Select recomputing method'); %#ok<NASGU>
       
    popup = uicontrol('Parent',d, 'Style','popup',...
           'Position',[40 70 170 25], 'String',options,...
           'Callback',@popup_callback); %#ok<NASGU>
       
    button = uicontrol('Parent',d, 'Position',[89 20 70 25],...
           'String','OK', 'Callback','delete(gcf)'); %#ok<NASGU>
       
    choice = 1;
    uiwait(d);
   
       function popup_callback(popup, event) %#ok<INUSD>
           choice = popup.Value;
%           idx = popup.Value;
%           popup_items = popup.String;
%           choice = char(popup_items(idx,:));
       end
   
end

% ------------------------------------------------------------------------
function recomputeFlow_edgeCutOff(~, ~, ~)

global im0
edges = im0.selectedEdgesList(2, :);
l = size(im0.cuoffMethd.out, 1);

msg = {'processing...', ['number of cycles is: ', num2str(l)]};
hw = waitbar(0, msg);
pause(0.3)

for itr = 1:l
    msg = {'processing...', ['cycle ', num2str(itr), 'out of', num2str(l)]};
    waitbar(0, hw, msg);
    
    values = num2cell(im0.cuoffMethd.out(itr, :));
    [idMthd, diamLrLim, diamUprLim, ctOffMthd, ctOff] = deal(values{:});
    
    ll = length(edges);
    for j = 1:ll
        e = edges(j);
        
        if idMthd == 1   %- diameter range
            d = im0.segDiam(im0.edgeSegN(e));
            if d < diamLrLim || d > diamUprLim; continue; end
            compute(ctOffMthd)
            
        elseif idMthd == 2  %- for segments labeled as arteries only
            vsslType = im0.segVesType(im0.edgeSegN(e));
            if vsslType ~= 1; continue; end
            compute(ctOffMthd)
            
        elseif idMthd == 3  %- for segments labeled as capillaries only
            vsslType = im0.segVesType(im0.edgeSegN(e));
            if vsslType ~= 2; continue; end
            compute(ctOffMthd)
            
        elseif idMthd == 4  %- for segments labeled as veins only
            vsslType = im0.segVesType(im0.edgeSegN(e));
            if vsslType ~= 3; continue; end
            compute(ctOffMthd)
        end
        
        waitbar(j/ll, hw, msg);
    end 
end
waitbar(1, hw, 'done!');
pause(0.3)
close(hw)

    function compute(cutOffMethod)
        
            if cutOffMethod == 1  %- top one half
                dataPt = getDataPt(e, 'cutAt', 0);
                ctOffValue = quantile(dataPt, 0.5); %- data median (Q2 value)
                dataMean = mean(dataPt);
                if dataMean > 0
                    flowV = dataPt(dataPt >= ctOffValue);
                else
                    flowV = dataPt(dataPt <= ctOffValue);
                end
                vScld = (im0.nodePos_um(im0.nodeEdges(e, 2), :) - im0.nodePos_um(im0.nodeEdges(e, 1), :));
                vScldMag = norm(vScld);
                cosPhZ = dot([0,0,1],vScld)/vScldMag;
                
                meanFlowV = mean(flowV)/abs(cosPhZ);
                im0.edgeMeanFlow(e) = meanFlowV;
                im0.selectedEdgesList(3, j) = meanFlowV;
                im0.cutOff(e) = ctOffValue;
                
            elseif cutOffMethod == 2  %- top one quarter (Q3 value)
                dataPt = getDataPt(e, 'cutAt', 0);
                ctOffValue = quantile(dataPt, 0.25);
                dataMean = mean(dataPt);
                if dataMean > 0
                    flowV = dataPt(dataPt >= ctOffValue);
                else
                    flowV = dataPt(dataPt <= ctOffValue);
                end
                vScld = (im0.nodePos_um(im0.nodeEdges(e, 2), :) - im0.nodePos_um(im0.nodeEdges(e, 1), :));
                vScldMag = norm(vScld);
                cosPhZ = dot([0,0,1],vScld)/vScldMag;
                meanFlowV = mean(flowV)/abs(cosPhZ);
                im0.edgeMeanFlow(e) = meanFlowV;
                im0.selectedEdgesList(3, j) = meanFlowV;
                im0.cutOff(e) = ctOffValue;
                
            elseif cutOffMethod == 3  %- fixed cutoff value
                dataPt = getDataPt(e, 'cutAt', 0);
                ctOffValue = ctOff;
                dataMean = mean(dataPt);
                if dataMean > 0
                    flowV = dataPt(dataPt >= ctOffValue);
                else
                    flowV = dataPt(dataPt <= ctOffValue);
                end
                vScld = (im0.nodePos_um(im0.nodeEdges(e, 2), :) - im0.nodePos_um(im0.nodeEdges(e, 1), :));
                vScldMag = norm(vScld);
                cosPhZ = dot([0,0,1],vScld)/vScldMag;
                meanFlowV = mean(flowV)/abs(cosPhZ);
                im0.edgeMeanFlow(e) = meanFlowV;
                im0.selectedEdgesList(3, j) = meanFlowV;
                im0.cutOff(e) = ctOffValue;
                %%% Qi Update 2023 March
                phZ = acos(abs(cosPhZ)) * 180/pi();
                meanFlowVel_Z = mean(flowV);
                stdFlowVel_Z = std(flowV);
                stdFlowV = std(flowV/abs(cosPhZ));
                im0.edgeStdFlow(e) = stdFlowV;
                im0.edgeMeanFlowVel_Z(e) = meanFlowVel_Z;
                im0.edgeStdFlowVel_Z(e) = stdFlowVel_Z;
                im0.cosPhZ(e) = cosPhZ;
                im0.phZ(e) = phZ;
                %%% Qi Update 2023 March
            end
    end
end



% --------------------------------------------------------------------
function recomputeFlow_segmentCutOff (~, ~, ~)
%%
global im0
sgmnts = unique(im0.selectedEdgesList(1, :));
l = size(im0.cuoffMethd.out, 1);

msg = {'processing...', ['number of cycles is: ', num2str(l)]};
hw = waitbar(0, msg);
pause(0.3)

for itr = 1:l
    %%
    msg = {'processing...', ['cycle ', num2str(itr), 'out of', num2str(l)]};
    waitbar(0, hw, msg);
    
    values = num2cell(im0.cuoffMethd.out(itr, :));
    [idMthd, diamLrLim, diamUprLim, ctOffMthd, ctOff] = deal(values{:});
    
    ll = length(sgmnts);
    for j = 1:ll
        s = sgmnts(j);
        e = im0.selectedEdgesList(2, im0.selectedEdgesList(1,:) == s);
        
        if idMthd == 1   %- diameter range
            d = im0.segDiam(s);
            if d < diamLrLim || d > diamUprLim; continue; end
            compute(ctOffMthd)
            
        elseif idMthd == 2  %- for segments labeled as arteries only
            vsslType = im0.segVesType(s);
            if vsslType ~= 1; continue; end
            compute(ctOffMthd)
            
        elseif idMthd == 3  %- for segments labeled as capillaries only
            vsslType = im0.segVesType(s);
            if vsslType ~= 2; continue; end
            compute(ctOffMthd)
            
        elseif idMthd == 4  %- for segments labeled as veins only
            vsslType = im0.segVesType(s);
            if vsslType ~= 3; continue; end
            compute(ctOffMthd)
        end
        
        waitbar(j/ll, hw, msg);
    end
end
waitbar(1, hw, 'done!');
pause(0.3)
close(hw)

    function compute(cutOffMethod)
        
            if cutOffMethod == 1  %- top one half
                dataPt = [];
                for ee = e
                    eDataPt = getDataPt(ee, 'cutAt', 0);
                    vScld = (im0.nodePos_um(im0.nodeEdges(ee, 2), :) - im0.nodePos_um(im0.nodeEdges(ee, 1), :));
                    vScldMag = norm(vScld);
                    cosPhZ = abs(dot([0,0,1],vScld)/vScldMag);
                    eDataPt = eDataPt / cosPhZ;
                    
                    dataPt = [dataPt; eDataPt]; %#ok<AGROW>
                end
                ctOffValue = quantile(dataPt, 0.5); %- data median (Q2 value)
                dataMean = mean(dataPt);
                if dataMean > 0
                    direction = 1;
                else
                    direction = -1;
                end
                dataPt = abs(dataPt);
                flowV = dataPt(dataPt >= ctOffValue);
                meanFlowV = mean(flowV) * direction;
                im0.segMeanFlowV(s) = meanFlowV;
                im0.edgeMeanFlow(e) = meanFlowV;
                im0.selectedEdgesList(3, im0.selectedEdgesList(1,:) == s) = meanFlowV;
                im0.cutOff(e) = ctOffValue;
                
            elseif cutOffMethod == 2  %- top one quarter (Q3 value)
                dataPt = [];
                for ee = e
                    eDataPt = getDataPt(ee, 'cutAt', 0);
                    vScld = (im0.nodePos_um(im0.nodeEdges(ee, 2), :) - im0.nodePos_um(im0.nodeEdges(ee, 1), :));
                    vScldMag = norm(vScld);
                    cosPhZ = abs(dot([0,0,1],vScld)/vScldMag);
                    eDataPt = eDataPt / cosPhZ;
                    
                    dataPt = [dataPt; eDataPt]; %#ok<AGROW>
                end
                ctOffValue = quantile(dataPt, 0.25); 
                dataMean = mean(dataPt);
                if dataMean > 0
                    direction = 1;
                else
                    direction = -1;
                end
                dataPt = abs(dataPt);
                flowV = dataPt(dataPt >= ctOffValue);
                meanFlowV = mean(flowV) * direction;
                im0.segMeanFlowV(s) = meanFlowV;
                im0.edgeMeanFlow(e) = meanFlowV;
                im0.selectedEdgesList(3, im0.selectedEdgesList(1,:) == s) = meanFlowV;
                im0.cutOff(e) = ctOffValue;
                
            elseif cutOffMethod == 3  %- fixed cutoff value
                dataPt = [];
                for ee = e
                    eDataPt = getDataPt(ee, 'cutAt', 0);
                    vScld = (im0.nodePos_um(im0.nodeEdges(ee, 2), :) - im0.nodePos_um(im0.nodeEdges(ee, 1), :));
                    vScldMag = norm(vScld);
                    cosPhZ = abs(dot([0,0,1],vScld)/vScldMag);
                    eDataPt = eDataPt / cosPhZ;
                    
                    dataPt = [dataPt; eDataPt]; %#ok<AGROW>
                end
                ctOffValue = ctOff;
                dataMean = mean(dataPt);
                if dataMean > 0
                    direction = 1;
                else
                    direction = -1;
                end
                dataPt = abs(dataPt);
                flowV = dataPt(dataPt >= ctOffValue);
                meanFlowV = mean(flowV) * direction;
                im0.segMeanFlowV(s) = meanFlowV;
                im0.edgeMeanFlow(e) = meanFlowV;
                im0.selectedEdgesList(3, im0.selectedEdgesList(1,:) == s) = meanFlowV;
                im0.cutOff(e) = ctOffValue;
                
            end
    end
end


% --------------------------------------------------------------------
function statsMnu_Callback(hObject, eventdata, handles)  
global im0
if ~isfield(im0, 'data')
    prcssSgmntFlwV_Callback(hObject, eventdata, handles)
end
stats
% hObject    handle to statsMnu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
end


function cutOff_Callback(hObject, ~, handles)  
% hObject    handle to cutOff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global im0
ctOff = str2double(get(hObject,'String'));
im0.cutOff(im0.gui.edgeID) = ctOff;
computeEdgeFlow(handles, im0.gui.edgeID)
% Hints: get(hObject,'String') returns contents of cutOff as text
%        str2double(get(hObject,'String')) returns contents of cutOff as a double
end

% --- Executes during object creation, after setting all properties.
function cutOff_CreateFcn(hObject, ~, ~)  
% hObject    handle to cutOff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

end
% --------------------------------------------------------------------
function save_TP_img_Callback(~, ~, handles)  
img = getframe(handles.axes_TP_Ang);
[file,path] = uiputfile('TPAngImg00.jpg');
imwrite(img.cdata, [path,file]);

% hObject    handle to save_TP_img (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

end
% --------------------------------------------------------------------
function saveOCTImg_Callback(~, ~, handles)  
img = getframe(handles.axes_OCT_Flow);
[file,path] = uiputfile('OCTFlowImg00.jpg');
imwrite(img.cdata, [path,file]);
% hObject    handle to saveOCTImg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


end
function SgmntDiam_Est_Callback(hObject, ~, ~)  
% hObject    handle to SgmntDiam_Est (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global im0
d = get(hObject,'String');
im0.segDiam(im0.gui.segmentID) = str2double(d);
set(hObject, 'String', [d, ' μm'])
end
% --- Executes during object creation, after setting all properties.
function SgmntDiam_Est_CreateFcn(hObject, ~, ~)  
% hObject    handle to SgmntDiam_Est (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

end

% --------------------------------------------------------------------
function sgmntType_Callback(hObject, ~, handles)  
% hObject    handle to sgmntType (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global im0

if isempty(im0.gui.edgeID)
    errordlg('Please select an edge or a segment first', 'Input Error')
    set(hObject, 'String', 'Ø')
    set(handles.txt_vslType, 'String', '    Ø')
    return
end

type = str2double(get(hObject,'String'));
if type > 3 || type < 0
    errordlg({'Vessel type could only be 1,2,3 or 0    ',...
        '1 for Artery', '2 for Capillary', '3 for Vein',...
        '0 for unspecified vessel type'},'Input Error')
    set(hObject, 'String', num2str(im0.segVesType(im0.gui.segmentID)))
    return
end

im0.segVesType(im0.gui.segmentID) = type;
if type == 0
    set(handles.txt_vslType, 'String', 'unspecified')
elseif type == 1
    set(handles.txt_vslType, 'String', 'Artery')
elseif type == 2
    set(handles.txt_vslType, 'String', 'Capillary')
elseif type == 3
    set(handles.txt_vslType, 'String', 'Vein')
end

updateSelectedEdgesList(handles)
% Hints: get(hObject,'String') returns contents of sgmntType as text
%        str2double(get(hObject,'String')) returns contents of sgmntType as a double
end

% --- Executes during object creation, after setting all properties.
function sgmntType_CreateFcn(hObject, ~, ~)  
% hObject    handle to sgmntType (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

% --------------------------------------------------------------------
function menu_voxSize_Callback(~, ~, ~)  
global im0

prompt = {'X dimension (μm):','Y dimension (μm):', 'Z dimension (μm)'};
dlgtitle = 'Voxel size';
dims = [1 20];
definput = {num2str(im0.Hvox(1)),num2str(im0.Hvox(2)), num2str(im0.Hvox(3))};
answer = inputdlg(prompt,dlgtitle,dims,definput);

if isempty(answer)
    return
end

X = str2double(answer{1});
Y = str2double(answer{2});
Z = str2double(answer{3});

im0.Hvox = [X, Y, Z];
end
% hObject    handle to menu_voxSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function rotAngle_Callback(hObject, ~, ~)
global im0
a = get(hObject,'String');
im0.gui.rotAngl = str2double(a);
set(hObject, 'String', [a, '°'])
% hObject    handle to rotAngle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
end
% Hints: get(hObject,'String') returns contents of rotAngle as text
%        str2double(get(hObject,'String')) returns contents of rotAngle as a double


% --- Executes during object creation, after setting all properties.
function rotAngle_CreateFcn(hObject, ~, ~)  
% hObject    handle to rotAngle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

%--------------------------------------
function dataPt = getDataPt(edgeID, method, cutOff)
%%- method can be 'all' which would return all datapoints in the cylinder
%%- or can be 'cutAt' which applies a threshold at zero.
%%- or can be 'cutRange'
%%- or can be 'allButZeros'

switch method
    case 'cutAt'
        if size(cutOff) ~= 1
            error('err:one', ['(ಠ_ಠ) heyyy budd!!!\n',...
                'I think you didn''t put proper cutoff value.\n',...
                'The value should be just one number for method ''cutAt'''])
            return %#ok<UNRCH>
        end
        
    case 'cutRange'
        if size(cutOff) ~= 2
            error('err:one', ['(ಠ_ಠ) heyyy budd!!!\n',...
                'I think you didn''t put proper cutoff value.\n',...
                'The value should be two numbers for method ''cutRange'''])
            return %#ok<UNRCH>
        end
end


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

sz = im0.szOCT;

if c0(1) < 1; c0(1) = 1; c1(1) = Vsize(1); end
if c0(2) < 1; c0(2) = 1; c1(2) = Vsize(2); end
if c0(3) < 1; c0(3) = 1; c1(3) = Vsize(3); end

if c1(1) > sz(1); c1(1) = sz(1); c0(1) = c1(1) - Vsize(1) + 1; end
if c1(2) > sz(2); c1(2) = sz(2); c0(2) = c1(2) - Vsize(2) + 1; end
if c1(3) > sz(3); c1(3) = sz(3); c0(3) = c1(3) - Vsize(3) + 1; end

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
    
    switch method
        case 'allButZeros'
            dataPt = dataPt(dataPt ~= 0);
            
        case 'cutAt'
            if cutOff == 0
                dataPt = dataPt(dataPt ~= 0);
            end
            dataMean = mean(dataPt);
            if dataMean > 0
                dataPt = dataPt(dataPt >= cutOff);
            else
                dataPt = dataPt(dataPt <= cutOff);
            end
            
        case 'cutRange'
            dataPt = dataPt(dataPt >= cutOff(1) & dataPt <= cutOff(2));
    end
end
if isempty(dataPt); dataPt=0; end
end


% --- Executes on button press in segHist_gui.
function segHist_gui_Callback(~, ~, ~)  
sgDataptHist

% hObject    handle to segHist_gui (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
end


% --------------------------------------------------------------------
function menu_labelSegments_Callback(~, ~, handles)
global im0

if isempty(im0.VsslLbeling.pial)
    menu_loadPial_Callback()
    if isempty(im0.VsslLbeling.pial)
        errordlg('process was aborted')
        return
    end
end
if isempty(im0.VsslLbeling.artsAndVens)
    menu_loadArtAndVen_Callback()
    if isempty(im0.VsslLbeling.artsAndVens)
        errordlg('process was aborted')
        return
    end
end


if ~isfield(im0, 'sgmntBrOrderArt') || ~isfield(im0, 'sgmntBrOrderVein')
    if ~isfield(im0, 'edgeBRorderVeins') || ~isfield(im0, 'edgeBRorderArt')
        menu_loadBrOr_Callback()
    else
        sgmntIDs = 1: length(im0.segDiam);
        im0.sgmntBrOrderArt = zeros(1, max(sgmntIDs(:)));
        im0.sgmntBrOrderVein = zeros(1, max(sgmntIDs(:)));
        
        for sg = sgmntIDs
            brOr = im0.edgeBRorderArt(im0.edgeSegN == sg);
            brOr = unique(brOr);
            im0.sgmntBrOrderArt(sg) = brOr;

            brOr = im0.edgeBRorderVeins(im0.edgeSegN == sg);
            brOr = unique(brOr);
            im0.sgmntBrOrderVein(sg) = brOr;
        end
    end
end


prompt = {'max capillary diameter (μm):'};
dlgtitle = 'Input diameter limit';
dims = [1 30];
definput = {'2'};
answer = cell2mat(inputdlg(prompt, dlgtitle, dims, definput));
lmt = str2double(answer);

if isempty(lmt) || lmt < 0
    errordlg({'process was aborted', 'due to improper diameter limit'}, 'input error')
    return
end

labelSegments(handles, lmt)
msgbox('labeling Completed');
% hObject    handle to menu_labelSegments (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
end




%%--
function producePlots(plt)
%%
%%-- plt is an integer or an array of integers which refers to the plot
%%-- requested.
%%-- 1 for Diamerter vs absolute flow V (all)
%%-- 2 for Diamerter vs absolute flow V (Arteries)
%%-- 3 for Diamerter vs absolute flow V (Veins)
%%-- 4 for Diamerter vs absolute flow V (Cap)
%%-- 5 for Branching Order From Artery Side
%%-- 6 for Branching Order From Vein Side
%%
global im0

if ~isfield(im0, 'data')
    prcssSgmntFlwV_Callback()
end

sgmntT = table2array(im0.data);
sgmntIDs = sgmntT(:,1)';
sgmntFloV = sgmntT(:,3)';
sgmntDiam = im0.segDiam(sgmntIDs);
sgmntType = im0.segVesType(sgmntIDs);
artSgmntBrOrdr = im0.sgmntBrOrderArt(sgmntIDs);
veinSgmntBrOrdr = im0.sgmntBrOrderVein(sgmntIDs);


%%
%%-- (1) Diamerter vs absolute flow V (all)
if ~isempty(intersect(1, plt))
    figure('Name','Diamerter vs absolute flow V (all)','NumberTitle','off');
    
    indx = sgmntType == 1;
    diams = sgmntDiam(indx);    flowVs = abs(sgmntFloV(indx));
    scatter(diams, flowVs, 13, 'filled', 'r')
    hold on
    
    indx = sgmntType == 3;
    diams = sgmntDiam(indx);    flowVs = abs(sgmntFloV(indx));
    scatter(diams, flowVs, 13, 'filled', 'b')
    
    indx = sgmntType == 2;
    diams = sgmntDiam(indx);    flowVs = abs(sgmntFloV(indx));
    scatter(diams, flowVs, 13, 'filled', 'g')
    
    legend('Arteries', 'Veins', 'Capillaries')
    hold off
    xlabel('Sgmnt Diam')
    ylabel('abs Velocity (mm/s)')
    title('Diam vs Flow V (All)')
end
%%

%%-- (2) Diamerter vs absolute flow V (Arteries)
if ~isempty(intersect(2, plt))
    indx = sgmntType == 1;
    diams = sgmntDiam(indx);    flowVs = abs(sgmntFloV(indx));
    
    %%%
    intrvl = [8, 10, 14, 18, 26]; i = 0;
    avgDiams = [];  avgFlowVs = []; err = [];
    for f = intrvl
        indx = find(diams > i & diams <= f);
        avgDiams = [avgDiams, mean(unique(diams(indx)))]; %#ok<AGROW>
        avgFlowVs = [avgFlowVs, mean(flowVs(indx))]; %#ok<AGROW>
        err = [err, std(flowVs(indx))/sqrt(length(flowVs(indx)))]; %#ok<AGROW>
        i = f;
    end
    %%%
    
    figure('Name','Diamerter vs absolute flow V (Arteries)','NumberTitle','off');
    scatter(diams, flowVs, 13, 'filled')
    hold on
    errorbar(avgDiams, avgFlowVs, err, 'o', 'MarkerSize', 10, 'LineWidth', 1.5)
    hold off
    xlabel('Artery Diameters')
    ylabel('abs Velocity (mm/s)')
    title('Diam vs Flow V (Arteries)')
end
%%

%%-- (3) Diamerter vs absolute flow V (Veins)

if ~isempty(intersect(3, plt))
    indx = sgmntType == 3;
    diams = sgmntDiam(indx);    flowVs = abs(sgmntFloV(indx));
    
    %%%
    intrvl = [10, 12, 14, 20, 28, 45]; i = 0;
    avgDiams = [];  avgFlowVs = []; err = [];
    for f = intrvl
        indx = find(diams > i & diams <= f);
        avgDiams = [avgDiams, mean(unique(diams(indx)))]; %#ok<AGROW>
        avgFlowVs = [avgFlowVs, mean(flowVs(indx))]; %#ok<AGROW>
        err = [err, std(flowVs(indx))/sqrt(length(flowVs(indx)))]; %#ok<AGROW>
        i = f;
    end
    %%%
    
    figure('Name','Diamerter vs absolute flow V (Veins)','NumberTitle','off');
    scatter(diams, flowVs, 13, 'filled')
    hold on
    errorbar(avgDiams, avgFlowVs, err, 'o', 'MarkerSize', 10, 'LineWidth', 1.5)
    hold off
    xlabel('Vein Diameters')
    ylabel('abs Velocity (mm/s)')
    title('Diam vs Flow V (Veins)')
end
%%

%%-- (4) Diamerter vs absolute flow V (Cap)

if ~isempty(intersect(4, plt))
    
    indx = sgmntType == 2;
    diams = sgmntDiam(indx);    flowVs = abs(sgmntFloV(indx));
    
    %%%
    intrvl = (1:14); i = 0;
    avgDiams = [];  avgFlowVs = []; err = [];
    for f = intrvl
        indx = find(diams > i & diams <= f);
        avgDiams = [avgDiams, mean(unique(diams(indx)))]; %#ok<AGROW>
        avgFlowVs = [avgFlowVs, mean(flowVs(indx))]; %#ok<AGROW>
        err = [err, std(flowVs(indx))/sqrt(length(flowVs(indx)))]; %#ok<AGROW>
        i = f;
    end
    %%%
    
    figure('Name','Diamerter vs absolute flow V (Cap)','NumberTitle','off');
    scatter(diams, flowVs, 13, 'filled')
    hold on
    errorbar(avgDiams, avgFlowVs, err, 'o', 'MarkerSize', 10, 'LineWidth', 1.5)
    hold off
    xlim([1.2, max(diams)+0.2])
    xlabel('Capillary Diameters')
    ylabel('abs Velocity (mm/s)')
    title('Diam vs Flow V (Cap)')
end
%%

%%-- (5) Branching Order From Artery Side
if ~isempty(intersect(5, plt))
    
    maxBrOrdr = 6;
    
    indx = sgmntType == 2;
    brOrdr = artSgmntBrOrdr(indx);  flowVs = abs(sgmntFloV(indx));
    
    flowVs(brOrdr > maxBrOrdr) = [];
    brOrdr(brOrdr > maxBrOrdr) = [];
    
    %%%
    intrvl = (1:maxBrOrdr); i = 0;
    avgOrdr = [];  avgFlowVs = []; err = [];
    for f = intrvl
        indx = find(brOrdr > i & brOrdr <= f);
        avgOrdr = [avgOrdr, mean(unique(brOrdr(indx)))]; %#ok<AGROW>
        avgFlowVs = [avgFlowVs, mean(flowVs(indx))]; %#ok<AGROW>
        err = [err, std(flowVs(indx))/sqrt(length(flowVs(indx)))]; %#ok<AGROW>
        i = f;
    end
    %%%
    
    
    figure('Name','Branching Order From Artery Side','NumberTitle','off');
    scatter(brOrdr, flowVs, 13, 'filled')
    hold on
    errorbar(avgOrdr, avgFlowVs, err, 'o', 'MarkerSize', 10, 'LineWidth', 1.5)
    hold off
    xlabel('Branch Order')
    ylabel('abs velocity (mm/s)')
    title('Brnch Ordr From Artery Side (Cap)')
    xlim([0.8, maxBrOrdr+0.2])
    xticks(1:maxBrOrdr)
end
%%

%%-- (6) Branching Order From Vein Side

if ~isempty(intersect(6, plt))
    
    maxBrOrdr = 6;
    
    indx = sgmntType == 2;
    brOrdr = veinSgmntBrOrdr(indx);    flowVs = abs(sgmntFloV(indx));
    
    flowVs(brOrdr > maxBrOrdr) = [];
    brOrdr(brOrdr > maxBrOrdr) = [];
    
    %%%
    intrvl = (1:maxBrOrdr); i = 0;
    avgOrdr = [];  avgFlowVs = []; err = [];
    for f = intrvl
        indx = find(brOrdr > i & brOrdr <= f);
        avgOrdr = [avgOrdr, mean(unique(brOrdr(indx)))]; %#ok<AGROW>
        avgFlowVs = [avgFlowVs, mean(flowVs(indx))]; %#ok<AGROW>
        err = [err, std(flowVs(indx))/sqrt(length(flowVs(indx)))]; %#ok<AGROW>
        i = f;
    end
    %%%
    
    figure('Name','Branching Order From Vein Side','NumberTitle','off');
    scatter(brOrdr, flowVs, 13, 'filled')
    hold on
    errorbar(avgOrdr, avgFlowVs, err, 'o', 'MarkerSize', 10, 'LineWidth', 1.5)
    hold off
    xlabel('Branch Order')
    ylabel('abs velocity (mm/s)')
    title('Brnch Ordr From Vein Side (Cap)')
    xlim([0.8, maxBrOrdr+0.2])
    xticks(1:maxBrOrdr)
end
end


% --------------------------------------------------------------------
function menuTitle_Plot_Callback(hObject, eventdata, handles)
% hObject    handle to menuTitle_Plot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
end

% --------------------------------------------------------------------
function diamFlowV_all_Callback(~, ~, ~)  
producePlots(1)
end
% hObject    handle to diamFlowV_all (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function diamFlowV_art_Callback(~, ~, ~)  
producePlots(2)
end
% hObject    handle to diamFlowV_art (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function diamFlowV_vein_Callback(~, ~, ~)  
producePlots(3)
end
% hObject    handle to diamFlowV_vein (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function diamFlowV_cap_Callback(~, ~, ~)  
producePlots(4)
end
% hObject    handle to diamFlowV_cap (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menu_brOrd_artSide_Callback(~, ~, ~)  
producePlots(5)
end
% hObject    handle to menu_brOrd_artSide (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menu_brOrd_veinSide_Callback(~, ~, ~)  
producePlots(6)
end
% hObject    handle to menu_brOrd_veinSide (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menu_plotAll_Callback(~, ~, ~)  
producePlots(1:6)
end
% hObject    handle to menu_plotAll (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menu_loadBrOr_Callback(~, ~, ~)  
global im0
[fN, fP] = uigetfile('*.mat','Select Branching Order File');
if fN == 0
    return
end
load([fP fN]);
im0.sgmntBrOrderArt = sgmntBrOrderArt;
im0.sgmntBrOrderVein = sgmntBrOrderVein;

msgbox('file loaded');
% hObject    handle to menu_loadBrOr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
end


% --------------------------------------------------------------------
function relabelArtVein_Callback(~, ~, handles)  
%%
global im0
sgmnts = im0.selectedEdgesList(1, :);
types = im0.selectedEdgesList(4, :);
flows = im0.selectedEdgesList(3, :);

artSg = sgmnts(types ~= 2 & flows < 0);
veinSg = sgmnts(types ~= 2 & flows > 0);

im0.segVesType(artSg) = 1;
im0.segVesType(veinSg) = 3;

updateSelectedEdgesList(handles)

if ~isempty(im0.VsslLbeling.pial)
    for sg = artSg
        im0.VsslLbeling.pial(im0.VsslLbeling.pial(:,1) == sg, 3) = 1;
    end
    for sg = veinSg
        im0.VsslLbeling.pial(im0.VsslLbeling.pial(:,1) == sg, 3) = 3;
    end
end
if ~isempty(im0.VsslLbeling.artsAndVens)
    for sg = artSg
        im0.VsslLbeling.artsAndVens(im0.VsslLbeling.artsAndVens(:,1) == sg, 3) = 1;
    end
    for sg = veinSg
        im0.VsslLbeling.artsAndVens(im0.VsslLbeling.artsAndVens(:,1) == sg, 3) = 3;
    end
end

msgbox('Operation Completed')
% hObject    handle to relabelArtVein (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
end



%%--
function updateGui_Listbox(handles)
global im0
if ~isempty(im0.selectedEdgesList)
    lst = im0.selectedEdgesList;
    lst(3,:) = round(lst(3,:), 4);
else
    lst = [];
end

lstVal = get(handles.selectedEdgesList,'Value');
set(handles.selectedEdgesList,'string',num2str(lst'));
set(handles.selectedEdgesList, 'Value', lstVal);
end


% --------------------------------------------------------------------
function menu_edgeAngleFilter_Callback(~, ~, handles)  
prompt = {'maximum angle to the Z axis (°):'};
dlgtitle = 'Input angle limit';
dims = [1 37];
definput = {'35'};
answer = cell2mat(inputdlg(prompt,dlgtitle,dims,definput));
phi = str2double(answer);

edgeFilter(phi)
updateGui_Listbox(handles)
msgbox('Operation Completed');
% hObject    handle to menu_edgeAngleFilter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
end

function edgeFilter(phi)
%%
global im0
f0 = 0;
f1 = size(im0.OCT_Flow, 3);

p0 = im0.nodeEdges(:, 1);
p1 = im0.nodeEdges(:, 2);


indx_ZFltr = (im0.nodePos(p0, 3) > f0 & im0.nodePos(p0, 3) <= f1 &...
    im0.nodePos(p1, 3) > f0 & im0.nodePos(p1, 3) <= f1);

vctrs = im0.nodePos_um(p1, :) - im0.nodePos_um(p0, :);
vctrsMag = vecnorm(vctrs, 2, 2);
Z = zeros(length(vctrs), 3); Z(:, 3) = 1;
cosPhZ = abs(dot(Z, vctrs,2) ./ vctrsMag);
phZ = acos(cosPhZ) * (180/pi());
indx_PhiFltr = abs(phZ) <= phi;

indx = indx_ZFltr & indx_PhiFltr;

edges = find(indx == true)';
if ~isempty(im0.selectedEdgesList)
    prcssEdges = intersect(edges, im0.selectedEdgesList(2,:));
else
    prcssEdges = [];
end

if ~isempty(prcssEdges)
    l = length(prcssEdges);
    indx = zeros(1, l);
    for i = 1:l
        indx(i) = find(edges == prcssEdges(i));
    end
    edges(indx) = [];
end


sgmnts = im0.edgeSegN(edges)';
types = im0.segVesType(sgmnts);

ext = zeros(4, length(edges));
ext([1,2,4], :) = [sgmnts; edges; types];

im0.selectedEdgesList = [im0.selectedEdgesList, ext];
end


% --------------------------------------------------------------------
function loadSession_ClickedCallback(~, ~, handles)  
openingDlg(handles)
end
% hObject    handle to loadSession (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



% --------------------------------------------------------------------
function refreshGui(handles)
%this function update displayed volues on the gui panels
global im0 V

l = 34;
hwait = waitbar(0, 'updating GUI panels');

%%-- update view plane for the TPMA window & OCT velocity windo
if im0.gui.viewPlane == 'XY'
    set(handles.viewPlane,'value', 1)
elseif im0.gui.viewPlane == 'YZ'
    set(handles.viewPlane,'value', 2)
else
    set(handles.viewPlane,'value', 3)
end
waitbar(1/l, hwait, 'updating GUI panels')

%%-- update rpojection type for all two window and four monitors
if im0.gui.projType == 'max'
    set(handles.viewPlane,'value', 1)
elseif im0.gui.projType == 'min'
    set(handles.projType,'value', 2)
else
    set(handles.projType,'value', 3)
end
waitbar(2/l, hwait, 'updating GUI panels')

%%-- update slider parameters
switch im0.gui.viewPlane
    case 'XY'
        sliderLim = max(size(V.TP_XY, 3), size(V.OCT_Flow_XY, 3));
    case 'YZ'
        sliderLim = max(size(V.TP_YZ, 3), size(V.OCT_Flow_YZ, 3));
    case 'XZ'
        sliderLim = max(size(V.TP_XZ, 3), size(V.OCT_Flow_XZ, 3));
end
sliderVal = im0.gui.frameNum - ceil(im0.gui.numOfFrames/2) + 1;
set(handles.slider,'max', sliderLim)
set(handles.slider,'min', 1)
set(handles.slider,'value', sliderVal)
set(handles.slider,'sliderstep', [1/sliderLim, 0.1])
waitbar(3/l, hwait, 'updating GUI panels')

%%-- update abs angle to Z for selected edge if there is one
waitbar(4/l, hwait, 'updating GUI panels')

%%-- update calculated edge flow velocity for selected edge if theris one
set(handles.edgeFlowV,'string',num2str(im0.gui.edgeFlowV))
waitbar(5/l, hwait, 'updating GUI panels')

%%--
updateGui_Listbox(handles)
waitbar(6/l, hwait, 'updating GUI panels')

%%--
if ~isempty(im0.selectedEdgesList) && ~isempty(im0.gui.edgeID)
    val = find(im0.selectedEdgesList(2,:)== im0.gui.edgeID);
    if ~isempty(val)
        set(handles.selectedEdgesList, 'Value', val)
    end
end
waitbar(7/l, hwait, 'updating GUI panels')

%%--
computeEdgeFlow(handles, im0.gui.edgeID)
waitbar(8/l, hwait, 'updating GUI panels')

%%--
waitbar(9/l, hwait, 'updating GUI panels')

%%--
waitbar(10/l, hwait, 'updating GUI panels')

%%--
waitbar(12/l, hwait, 'updating GUI panels')

%%--
waitbar(13/l, hwait, 'updating GUI panels')
switch im0.gui.angWind
    case 'XY'
        set(handles.radiobutton_XY, 'Value', 1)
    case 'YZ'
        set(handles.radiobutton_YZ, 'Value', 1)
    case 'XZ'
        set(handles.radiobutton_XZ, 'Value', 1)
end

%%--
waitbar(14/l, hwait, 'updating GUI panels')

%%--
waitbar(15/l, hwait, 'updating GUI panels')

%%--
updateGui_edgeInfo(handles)
waitbar(16/l, hwait, 'updating GUI panels')

%%-- update the number of the center frame of the displayed projection
%%-- images
set(handles.frameNum,'string',num2str(im0.gui.frameNum))
waitbar(17/l, hwait, 'updating GUI panels')

%%-- update upper limit value for the colorbar
set(handles.colorbarUpperLim,'string',num2str(im0.gui.colorbarUpperLim))
waitbar(18/l, hwait, 'updating GUI panels')

%%-- update lower limit value for the colorbar
set(handles.colorbalLowerLim,'string',num2str(im0.gui.colorbalLowerLim))
waitbar(19/l, hwait, 'updating GUI panels')

%%-- update the TPMA image threashold value
set(handles.TPthreshold,'string',num2str(im0.gui.TPthreshold))
waitbar(20/l, hwait, 'updating GUI panels')

%%-- update the number of frames projected in the current images on all
%%-- windows and monitors
set(handles.numOfFrames,'string',num2str(im0.gui.numOfFrames))
waitbar(21/l, hwait, 'updating GUI panels')

%%-- update the alpha value (transparency value) of the cylinder image
%%-- projected on the monitors
set(handles.alpha,'string',num2str(im0.gui.alpha))
waitbar(22/l, hwait, 'updating GUI panels')

%%-- update the width value of graph edges ploted on the TPMA and OCT
%%-- windows
set(handles.edgeLinewidth_axes,'string',...
    [num2str(im0.gui.edgeLinewidth_axes), ' pt'])
waitbar(23/l, hwait, 'updating GUI panels')

%%-- update the width value of graph nodes ploted on the TPMA and OCT
%%-- windows
set(handles.nodePointWidth_axes,'string',...
    [num2str(im0.gui.nodePointWidth_axes), ' pt'])
waitbar(24/l, hwait, 'updating GUI panels')

%%-- update the number of steps taken when moving the cylinder in
%%-- anydirection using the arrow pushbuttons
set(handles.step,'string', [num2str(im0.gui.step), ' pixel'])
waitbar(25/l, hwait, 'updating GUI panels')

%%-- update the number of pixels in image edge projected onto the monitors
%%-- (the monitors are square)
set(handles.winSize,'string', [num2str(im0.gui.winSize), ' pixel'])
waitbar(26/l, hwait, 'updating GUI panels')

%%-- update edge width of edge ploted on the monitors
set(handles.edgeLinewidth_monitors,'string',...
    [num2str(im0.gui.edgeLinewidth_monitors), ' pt'])
waitbar(27/l, hwait, 'updating GUI panels')

%%-- update the width of the nodes ploted on the monitors
set(handles.nodePointWidth_monitors,'string',...
    [num2str(im0.gui.nodePointWidth_monitors), ' pt'])
waitbar(28/l, hwait, 'updating GUI panels')

%%-- update the check box state which keeps the upper and lower colorbar
%%-- limits symmetric about 0 when it's checked and free when its not
set(handles.colorbarCheck, 'Value', im0.gui.colorbarCheck)
waitbar(29/l, hwait, 'updating GUI panels')

%%--
if all(im0.gui.graphDispState(1:2))
    handles.menu_dispGraph_TPMAPandOCTVP.Checked = 'on';
    handles.menu_dispGraph_TPMAP.Checked = 'on';
    handles.cntxtMenu_dispGraph_TPMA.Checked = 'on';
    handles.menu_dispGraph_OCTVP.Checked = 'on';
    handles.cntxtMenu_dispGraph_OCTVP.Checked = 'on';
else
    handles.menu_dispGraph_TPMAPandOCTVP.Checked = 'off';
    if im0.gui.graphDispState(1)
        handles.menu_dispGraph_TPMAP.Checked = 'on';
        handles.cntxtMenu_dispGraph_TPMA.Checked = 'on';
    else
        handles.menu_dispGraph_TPMAP.Checked = 'off';
        handles.cntxtMenu_dispGraph_TPMA.Checked = 'off';
    end
    if im0.gui.graphDispState(2)
        handles.menu_dispGraph_OCTVP.Checked = 'on';
        handles.cntxtMenu_dispGraph_OCTVP.Checked = 'on';
    else
        handles.menu_dispGraph_OCTVP.Checked = 'off';
        handles.cntxtMenu_dispGraph_OCTVP.Checked = 'off';
    end
end
waitbar(30/l, hwait, 'updating GUI panels')

%%--
xlim(handles.axes_TP_Ang, im0.gui.axesXlim)
ylim(handles.axes_TP_Ang, im0.gui.axesYlim)

xlim(handles.axes_OCT_Flow, im0.gui.axesXlim)
ylim(handles.axes_OCT_Flow, im0.gui.axesYlim)
waitbar(31/l, hwait, 'updating GUI panels')

%%--
if im0.gui.positionFlag(1)
    handles.menu_posFlag.Checked = 'on';
    handles.cntxtMenu_posFlag.Checked = 'on';
else
    handles.menu_posFlag.Checked = 'off';
    handles.cntxtMenu_posFlag.Checked = 'off';
end
waitbar(32/l, hwait, 'updating GUI panels')

%%--
switch im0.gui.selectionType
    case 'selectOneEdge'
        handles.cntxtMenu_selectOneEdge.Checked = 'on';
        handles.cntxtMenu_selectOneSgmnt.Checked = 'off';
        handles.cntxtMenu_selectMultEdges.Checked = 'off';
        handles.cntxtMenu_selectMultSgmnts.Checked = 'off';
        
    case 'selectOneSgmnt'
        handles.cntxtMenu_selectOneEdge.Checked = 'off';
        handles.cntxtMenu_selectOneSgmnt.Checked = 'on';
        handles.cntxtMenu_selectMultEdges.Checked = 'off';
        handles.cntxtMenu_selectMultSgmnts.Checked = 'off';
        
    case 'selectMultEdges'
        handles.cntxtMenu_selectOneEdge.Checked = 'off';
        handles.cntxtMenu_selectOneSgmnt.Checked = 'off';
        handles.cntxtMenu_selectMultEdges.Checked = 'on';
        handles.cntxtMenu_selectMultSgmnts.Checked = 'off';
        
    case 'selectMultSgmnts'
        handles.cntxtMenu_selectOneEdge.Checked = 'off';
        handles.cntxtMenu_selectOneSgmnt.Checked = 'off';
        handles.cntxtMenu_selectMultEdges.Checked = 'off';
        handles.cntxtMenu_selectMultSgmnts.Checked = 'on';
end
waitbar(33/l, hwait, 'updating GUI panels')



%%--
set(handles.rotAngle, 'String', [num2str(im0.gui.rotAngl), '°'])
waitbar(34/l, hwait, 'updating GUI panels')



waitbar(1,hwait,'Done!')
pause(0.2)
close(hwait)
end



%%-------------------------------------------------------------------------
function initialiseVar()
%initialiseVar() initialises all variables needed in the processes of the
%software.

l = 46;
global im0 V


hwait = waitbar(0, 'initializing variables');


if ~isfield(im0,'gui')
    im0.gui = struct();
end


if ~isfield(im0,'edgeMeanFlow')
    im0.edgeMeanFlow = zeros(1, length(im0.nodeEdges));
end
waitbar(1/l, hwait, 'initializing variables')

%%-------------------------------------------------------------------------
if ~isfield(im0,'edgeDiaEst')
    if isfield(im0,'nodeDiam')
        im0.edgeDiaEst = median(im0.nodeDiam(im0.nodeEdges), 2);
    else
        im0.edgeDiaEst = zeros(length(im0.nodeEdges), 1);
    end
end
waitbar(2/l, hwait, 'initializing variables')

%%-------------------------------------------------------------------------
if ~isfield(im0,'cylDia') || isempty(im0.cylDia)
    im0.cylDia = im0.edgeDiaEst;
end
waitbar(3/l, hwait, 'initializing variables')

%%-------------------------------------------------------------------------
if ~isfield(im0,'segDiam') || (max(im0.segDiam(:)) <= 1)
    if ~isfield(im0,'segNedges')
        im0.segDiam = [];
        im0.segNedges = [];
    else
        ll = length(im0.segNedges);
        for ii = 1:ll
            diams = im0.nodeDiam(im0.nodeSegN == ii);
            im0.segDiam(ii) = median(diams);
        end
    end
end
waitbar(4/l, hwait, 'initializing variables')

%%-------------------------------------------------------------------------
if ~isfield(im0,'segMeanFlow')
    if ~isfield(im0,'nodeSegN')
        im0.nodeSegN = [];
    end
    
    if max(im0.edgeMeanFlow) == 0
        im0.segMeanFlowV = zeros(1,max(im0.nodeSegN(:)));
    else
        ll = length(im0.segNedges);
        for ii = 1:ll
            flowVs = im0.edgeMeanFlow(im0.edgeSegN == ii);
            im0.segMeanFlowV(ii) = mean(flowVs);
        end
    end
end
if ~isfield(im0,'edgeFlag')
    im0.edgeFlag = zeros(1,length(im0.nodeEdges));
end
waitbar(5/l, hwait, 'initializing variables')

%%-------------------------------------------------------------------------
if ~isfield(im0, 'cutOff')
    im0.cutOff = zeros(1,length(im0.nodeEdges));
end
waitbar(6/l, hwait, 'initializing variables')

%%-------------------------------------------------------------------------
if ~isfield(im0.gui,'viewPlane')
    im0.gui.viewPlane = 'XY';
end

waitbar(7/l, hwait, 'initializing variables')


%  12-Apr-2023 Qi & Allen mods start 
if ~isfield(im0,'edgeStdFlow')
    im0.edgeStdFlow = zeros(1,length(im0.nodeEdges));
end

if ~isfield(im0,'edgeMeanFlowVel_Z')
    im0.edgeMeanFlowVel_Z = zeros(1,length(im0.nodeEdges));
end

if ~isfield(im0,'edgeStdFlowVel_Z')
    im0.edgeStdFlowVel_Z = zeros(1,length(im0.nodeEdges));
end

if ~isfield(im0,'cosPhZ')
    im0.cosPhZ = zeros(1,length(im0.nodeEdges));
end

if ~isfield(im0,'phZ')
    im0.phZ = zeros(1,length(im0.nodeEdges));
end
 %  12-Apr-2023 Qi & Allen mods end 
%%-------------------------------------------------------------------------
if ~isfield(im0,'selectedEdgesList')
    im0.selectedEdgesList = [];
end
waitbar(8/l, hwait, 'initializing variables')

%%-------------------------------------------------------------------------
waitbar(9/l, hwait, 'initializing variables')

%%-------------------------------------------------------------------------
if ~isfield(im0.gui,'edgeFlowV')
    im0.gui.edgeFlowV = [];
end
waitbar(10/l, hwait, 'initializing variables')

%%-------------------------------------------------------------------------
waitbar(11/l, hwait, 'initializing variables')

%%-------------------------------------------------------------------------
if ~isfield(im0.gui,'edgeID')
    im0.gui.edgeID =[];
end
waitbar(12/l, hwait, 'initializing variables')

%%-------------------------------------------------------------------------
if ~isfield(im0.gui,'segmentID')
    im0.gui.segmentID = [];
end
waitbar(13/l, hwait, 'initializing variables')

%%-------------------------------------------------------------------------
if ~isfield(im0,'VsslLbeling')
    im0.VsslLbeling.pial = [];
    im0.VsslLbeling.artsAndVens = [];
end
waitbar(14/l, hwait, 'initializing variables')

%%-------------------------------------------------------------------------
if ~isfield(im0,'sgmntBrOrderArt')
    im0.sgmntBrOrderArt = zeros(1, length(im0.segDiam));
end
waitbar(15/l, hwait, 'initializing variables')

%%-------------------------------------------------------------------------
if ~isfield(im0,'sgmntBrOrderVein')
    im0.sgmntBrOrderVein = zeros(1, length(im0.segDiam));
end
waitbar(16/l, hwait, 'initializing variables')

%%-------------------------------------------------------------------------
if ~isfield(im0.gui,'edgeFlagTemp')
    im0.gui.edgeFlagTemp = [];
end
waitbar(17/l, hwait, 'initializing variables')

%%-------------------------------------------------------------------------
waitbar(18/l, hwait, 'initializing variables')

%%-------------------------------------------------------------------------
waitbar(19/l, hwait, 'initializing variables')

%%-------------------------------------------------------------------------
waitbar(20/l, hwait, 'initializing variables')

%%-------------------------------------------------------------------------
if ~isfield(im0.gui,'frameNum')
    im0.gui.frameNum = 6;
end
waitbar(21/l, hwait, 'initializing variables')

%%-------------------------------------------------------------------------
if ~isfield(im0.gui,'colorbarUpperLim')
    im0.gui.colorbarUpperLim = 5;
end
waitbar(22/l, hwait, 'initializing variables')

%%-------------------------------------------------------------------------
if ~isfield(im0.gui,'colorbalLowerLim')
    im0.gui.colorbalLowerLim = -5;
end
waitbar(23/l, hwait, 'initializing variables')

%%-------------------------------------------------------------------------
if ~isfield(im0.gui,'TPthreshold')
    im0.gui.TPthreshold = 0.1;
end
waitbar(24/l, hwait, 'initializing variables')

%%-------------------------------------------------------------------------
if ~isfield(im0.gui,'numOfFrames')
    im0.gui.numOfFrames = 11;
end
waitbar(25/l, hwait, 'initializing variables')

%%-------------------------------------------------------------------------
if ~isfield(im0.gui,'projType')
    im0.gui.projType = 'max';
end
waitbar(26/l, hwait, 'initializing variables')

%%-------------------------------------------------------------------------
if ~isfield(im0.gui,'alpha')
    im0.gui.alpha = 0.5;
end
waitbar(27/l, hwait, 'initializing variables')

%%-------------------------------------------------------------------------
if ~isfield(im0.gui,'edgeLinewidth_axes')
    im0.gui.edgeLinewidth_axes = 0.5;
end
waitbar(28/l, hwait, 'initializing variables')

%%-------------------------------------------------------------------------
if ~isfield(im0.gui,'nodePointWidth_axes')
    im0.gui.nodePointWidth_axes = 4;
end
waitbar(29/l, hwait, 'initializing variables')

%%-------------------------------------------------------------------------
if ~isfield(im0.gui,'step')
    im0.gui.step = 1;
end
waitbar(30/l, hwait, 'initializing variables')

%%-------------------------------------------------------------------------
if ~isfield(im0.gui,'angWind')
    im0.gui.angWind = 'XY';
end
waitbar(31/l, hwait, 'initializing variables')

%%-------------------------------------------------------------------------
if ~isfield(im0.gui, 'mntrPt')
   im0.gui.mntrPt = [256, 256, 1]; 
end
waitbar(32/l, hwait, 'initializing variables')

%%-------------------------------------------------------------------------
if ~isfield(im0.gui,'winSize')
    im0.gui.winSize = 40;
end
waitbar(33/l, hwait, 'initializing variables')

%%-------------------------------------------------------------------------
if ~isfield(im0.gui,'edgeLinewidth_monitors')
    im0.gui.edgeLinewidth_monitors = 1;
end
waitbar(34/l, hwait, 'initializing variables')

%%-------------------------------------------------------------------------
if ~isfield(im0.gui,'nodePointWidth_monitors')
    im0.gui.nodePointWidth_monitors = 5;
end
waitbar(35/l, hwait, 'initializing variables')

%%-------------------------------------------------------------------------
if ~isfield(im0.gui,'axesXlim') || im0.gui.axesXlim(2) > size(im0.I, 1) ||...
        im0.gui.axesYlim(2) > size(im0.I, 2)
    im0.gui.axesXlim = [1, size(im0.I, 1)];
    im0.gui.axesYlim = [1, size(im0.I, 2)];
end
waitbar(36/l, hwait, 'initializing variables')

%%-------------------------------------------------------------------------
if ~isfield(im0,'szOCT')
    im0.szOCT = size(im0.OCT_Flow);
end
waitbar(37/l, hwait, 'initializing variables')

%%-------------------------------------------------------------------------
if ~isfield(im0.gui,'colorbalCheck')
    im0.gui.colorbarCheck = 1;
end
waitbar(38/l, hwait, 'initializing variables')

%%-------------------------------------------------------------------------
if ~isfield(im0, 'cylNodes') || isempty(im0.cylNodes)
    l = length(im0.nodeEdges);
    edgeP = im0.nodeEdges(1:l,:);
    p0 = im0.nodePos(edgeP(:,1)',:);
    p1 = im0.nodePos(edgeP(:,2)',:);
    im0.cylNodes = [p0,p1];
end
waitbar(39/l, hwait, 'initializing variables')

%%-------------------------------------------------------------------------
% if ~isfield(im0.gui,'slider')
%     im0.gui.slider = 6;
% end
waitbar(40/l, hwait, 'initializing variables')

%%-------------------------------------------------------------------------
if ~isfield(im0.gui,'flowClrMp')
    cm1a = (32:-1:1)' * [0 1 0] / 32;
    cm1a(:,3) = 1;
    cm1b = (32:-1:1)' * [0 0 1] / 32;
    cm2 = hot(64);
    im0.gui.flowClrMp = [cm1a; cm1b; cm2];
end
waitbar(41/l, hwait, 'initialization complete')
pause(0.3)

%%-------------------------------------------------------------------------
if ~isfield(im0.gui,'graphDispState')
    im0.gui.graphDispState = [true, true];
end
waitbar(42/l, hwait, 'initializing variables')

%%-------------------------------------------------------------------------
if ~isfield(im0.gui, 'positionFlag')
    im0.gui.positionFlag = [false, false];
end
waitbar(43/l, hwait, 'initializing variables')

%%-------------------------------------------------------------------------
if ~isfield(im0.gui, 'selectionType')
    im0.gui.selectionType = 'selectOneEdge';
end
waitbar(44/l, hwait, 'initializing variables')

%%-------------------------------------------------------------------------
if ~isfield(im0.gui, 'maxVolumeSz')
    im0.gui.maxVolumeSz = max(size(im0.I), size(im0.OCT_Flow));
end
waitbar(45/l, hwait, 'initializing variables')

%%-------------------------------------------------------------------------
if ~isfield(im0.gui, 'rotAngl')
    im0.gui.rotAngl = 10;
end
waitbar(46/l, hwait, 'initializing variables')



%%-------------------------------------------------------------------------
l = 6;
im0.I = max(uint8(im0.I),1);
waitbar(0.1 , hwait, {'Processing TPM angiogram projections', 'XY Projection'})
V.TP_XY = intensProj(im0.I, im0.gui.numOfFrames, 'XY', im0.gui.projType);

waitbar(1/l, hwait, {'Processing TPM angiogram projections', 'YZ Projection'})
V.TP_YZ = intensProj(im0.I, im0.gui.numOfFrames, 'YZ', im0.gui.projType);

waitbar(2/l, hwait, {'Processing TPM angiogram projections', 'XZ Projection'})
V.TP_XZ = intensProj(im0.I, im0.gui.numOfFrames, 'XZ', im0.gui.projType);
%%-------------------------------------------------------------------------
waitbar(3/l, hwait, {'Processing OCT velocity projections', 'XY Projection'})
V.OCT_Flow_XY = intensProj(im0.OCT_Flow, im0.gui.numOfFrames, 'XY', im0.gui.projType);

waitbar(4/l, hwait, {'Processing OCT velocity projections', 'YZ Projection'})
V.OCT_Flow_YZ = intensProj(im0.OCT_Flow, im0.gui.numOfFrames, 'YZ', im0.gui.projType);

waitbar(5/l, hwait, {'Processing OCT velocity projections', 'XZ Projection'})
V.OCT_Flow_XZ = intensProj(im0.OCT_Flow, im0.gui.numOfFrames, 'XZ', im0.gui.projType);
%%-------------------------------------------------------------------------
waitbar(0.95,hwait,'recording graph projections')
V.vctrs = vInfo(im0.nodeEdges, im0.nodePos, im0.gui.viewPlane);
[V.X, V.Y, V.Z] = edgeGraph(im0.nodeEdges, im0.nodePos);
V.projVctrs = [];
%%-------------------------------------------------------------------------
if ~isfield(V, 'prjVolumeSz')
    
    V.prjVolSz_OCT = [size(V.OCT_Flow_YZ, 3), size(V.OCT_Flow_XZ, 3), size(V.OCT_Flow_XY, 3)];
    V.prjVolSz_TPM = [size(V.TP_YZ, 3), size(V.TP_XZ, 3), size(V.TP_XY, 3)];
    V.prjVolumeSz = max(V.prjVolSz_OCT, V.prjVolSz_TPM);
end


waitbar(1,hwait,'Done!')
pause(0.2)
close(hwait)

end


% --------------------------------------------------------------------
function menuParent_view_Callback(~, ~, ~)
% hObject    handle to menuParent_view (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
end


% --------------------------------------------------------------------
function menu_graph_Callback(~, ~, ~)
% hObject    handle to menu_graph (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
end



% --------------------------------------------------------------------
function menu_dispGraph_TPMAPandOCTVP_Callback(hObject, ~, handles)  
global im0
val = hObject.Checked;
switch val
    case 'on'
        im0.gui.graphDispState(1:2) = false;
        handles.menu_dispGraph_TPMAPandOCTVP.Checked = 'off';
        handles.menu_dispGraph_TPMAP.Checked = 'off';
        handles.cntxtMenu_dispGraph_TPMA.Checked = 'off';
        handles.menu_dispGraph_OCTVP.Checked = 'off';
        handles.cntxtMenu_dispGraph_OCTVP.Checked = 'off';
    case 'off'
        im0.gui.graphDispState(1:2) = true;
        handles.menu_dispGraph_TPMAPandOCTVP.Checked = 'on';
        handles.menu_dispGraph_TPMAP.Checked = 'on';
        handles.cntxtMenu_dispGraph_TPMA.Checked = 'on';
        handles.menu_dispGraph_OCTVP.Checked = 'on';
        handles.cntxtMenu_dispGraph_OCTVP.Checked = 'on';
end
update_TPMA_image(handles)
update_OCTV_image(handles)
overlayGraph(handles)
end
% hObject    handle to menu_dispGraph_TPMAPandOCTVP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menu_dispGraph_TPMAP_Callback(hObject, ~, handles)
global im0
val = hObject.Checked;
switch val
    case 'on'
        im0.gui.graphDispState(1) = false;
        handles.menu_dispGraph_TPMAP.Checked = 'off';
        handles.cntxtMenu_dispGraph_TPMA.Checked = 'off';
        handles.menu_dispGraph_TPMAPandOCTVP.Checked = 'off';
    case 'off'
        im0.gui.graphDispState(1) = true;
        handles.menu_dispGraph_TPMAP.Checked = 'on';
        handles.cntxtMenu_dispGraph_TPMA.Checked = 'on';
        if all(im0.gui.graphDispState)
            handles.menu_dispGraph_TPMAPandOCTVP.Checked = 'on';
        end
end

update_TPMA_image(handles)
update_OCTV_image(handles)
overlayGraph(handles)
end
% hObject    handle to menu_dispGraph_TPMAP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menu_dispGraph_OCTVP_Callback(hObject, ~, handles)
global im0
val = hObject.Checked;
switch val
    case 'on'
        im0.gui.graphDispState(2) = false;
        handles.menu_dispGraph_OCTVP.Checked = 'off';
        handles.cntxtMenu_dispGraph_OCTVP.Checked = 'off';
        handles.menu_dispGraph_TPMAPandOCTVP.Checked = 'off';
    case 'off'
        im0.gui.graphDispState(2) = true;
        handles.menu_dispGraph_OCTVP.Checked = 'on';
        handles.cntxtMenu_dispGraph_OCTVP.Checked = 'on';
        if all(im0.gui.graphDispState)
            handles.menu_dispGraph_TPMAPandOCTVP.Checked = 'on';
        end
end

update_TPMA_image(handles)
update_OCTV_image(handles)
overlayGraph(handles)
end
% hObject    handle to menu_dispGraph_OCTVP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function cntxtMenu_dispGraph_TPMA_Callback(hObject, eventdata, handles)
menu_dispGraph_TPMAP_Callback(hObject, eventdata, handles)
end
% hObject    handle to cntxtMenu_dispGraph_TPMA (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function cntxtMenu_dispGraph_OCTVP_Callback(hObject, eventdata, handles)  
menu_dispGraph_OCTVP_Callback(hObject, eventdata, handles)
end
% hObject    handle to cntxtMenu_dispGraph_OCTVP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function cntxtMenu_posFlag_Callback(hObject, eventdata, handles)  
menu_posFlag_Callback(hObject, eventdata, handles)
end
% hObject    handle to cntxtMenu_posFlag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menu_posFlag_Callback(hObject, ~, handles)
global im0
val = hObject.Checked;
switch val
    case 'on'
        im0.gui.positionFlag(1) = false;
        handles.menu_posFlag.Checked = 'off';
        handles.cntxtMenu_posFlag.Checked = 'off';
        
        update_TPMA_image(handles)
        update_OCTV_image(handles)
        overlayGraph(handles)
        
    case 'off'
        im0.gui.positionFlag(1) = true;
        handles.menu_posFlag.Checked = 'on';
        handles.cntxtMenu_posFlag.Checked = 'on';
        
        if all(im0.gui.positionFlag)
            positionFlag_OCTFlow(handles, im0.gui.mntrPt)
            positionFlag_TPAng(handles, im0.gui.mntrPt)
        end
end
end
% hObject    handle to menu_posFlag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function cntxtMenu_clrSlction_Callback(~, ~, handles)
global im0
im0.gui.positionFlag(2) = false;
im0.edgeFlag(im0.edgeFlag == 2) = 0;
update_TPMA_image(handles)
update_OCTV_image(handles)
overlayGraph(handles)
update_Monitors(handles, [], true)
im0.gui.edgeFlowV = [];
im0.gui.edgeID = [];
im0.gui.segmentID = [];

set(handles.numOfDatapoints, 'string', 'Ø');
set(handles.angleToZ,'string', 'Ø');
set(handles.edgeFlowV,'string', 'Ø');
set(handles.maxVelZ,'string', 'Ø');
set(handles.meanVelZ,'string', 'Ø');
set(handles.minVelZ,'string', 'Ø');
set(handles.cylLength,'string', 'Ø');

updateGui_edgeInfo(handles)

histogram(handles.axes_Histogram, [], 30, 'Orientation', 'horizontal')
enableDefaultInteractivity(handles.axes_Histogram)

end
% hObject    handle to cntxtMenu_clrSlction (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function cntxtMenuParent_slctTyp_Callback(~, ~, ~)
% hObject    handle to cntxtMenuParent_slctTyp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
end

% --------------------------------------------------------------------
function cntxtMenu_selectOneEdge_Callback(~, ~, handles)
global im0
im0.gui.selectionType = 'selectOneEdge';
handles.cntxtMenu_selectOneEdge.Checked = 'on';
handles.cntxtMenu_selectOneSgmnt.Checked = 'off';
handles.cntxtMenu_selectMultEdges.Checked = 'off';
handles.cntxtMenu_selectMultSgmnts.Checked = 'off';
% hObject    handle to cntxtMenu_selectOneEdge (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
end

% --------------------------------------------------------------------
function cntxtMenu_selectOneSgmnt_Callback(~, ~, handles)
global im0
im0.gui.selectionType = 'selectOneSgmnt';
handles.cntxtMenu_selectOneEdge.Checked = 'off';
handles.cntxtMenu_selectOneSgmnt.Checked = 'on';
handles.cntxtMenu_selectMultEdges.Checked = 'off';
handles.cntxtMenu_selectMultSgmnts.Checked = 'off';
% hObject    handle to cntxtMenu_selectOneSgmnt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
end


% --------------------------------------------------------------------
function cntxtMenu_selectMultEdges_Callback(~, ~, handles)
global im0
im0.gui.selectionType = 'selectMultEdges';
handles.cntxtMenu_selectOneEdge.Checked = 'off';
handles.cntxtMenu_selectOneSgmnt.Checked = 'off';
handles.cntxtMenu_selectMultEdges.Checked = 'on';
handles.cntxtMenu_selectMultSgmnts.Checked = 'off';
% hObject    handle to cntxtMenu_selectMultEdges (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
end

% --------------------------------------------------------------------
function cntxtMenu_selectMultSgmnts_Callback(~, ~, handles)
global im0
im0.gui.selectionType = 'selectMultSgmnts';
handles.cntxtMenu_selectOneEdge.Checked = 'off';
handles.cntxtMenu_selectOneSgmnt.Checked = 'off';
handles.cntxtMenu_selectMultEdges.Checked = 'off';
handles.cntxtMenu_selectMultSgmnts.Checked = 'on';
% hObject    handle to cntxtMenu_selectMultSgmnts (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
end


% --------------------------------------------------------------------
function cntxtMenu_clrPosFlg_Callback(~, ~, handles)
global im0
im0.gui.positionFlag(2) = false;

update_Monitors(handles, [], true)
update_TPMA_image(handles)
update_OCTV_image(handles)
overlayGraph(handles)

% hObject    handle to cntxtMenu_clrPosFlg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
end

%%---------------------------------------------------------------------
function updateGui_edgeInfo(handles)
global im0
if ~isempty(im0.gui.edgeID)
    set(handles.edgeID, 'string', num2str(im0.gui.edgeID))
    set(handles.segmentID, 'string', num2str(im0.gui.segmentID))
    set(handles.cylDiam, 'String',...
        [num2str(im0.cylDia(im0.gui.edgeID)), ' μm']);
    set(handles.SgmntDiam_Est, 'String',...
        [num2str(im0.segDiam(im0.edgeSegN(im0.gui.edgeID))), ' μm'])
    set(handles.edgeDiam_Est, 'string',...
        [num2str(im0.edgeDiaEst(im0.gui.edgeID)), ' μm'])
    set(handles.cutOff, 'String', num2str(im0.cutOff(im0.gui.edgeID)))
    set(handles.brOrdr_art, 'string', num2str(im0.sgmntBrOrderArt(im0.gui.segmentID)))
    set(handles.brOrdr_vein, 'string', num2str(im0.sgmntBrOrderVein(im0.gui.segmentID)))
    
    type = im0.segVesType(im0.gui.segmentID);
    set(handles.sgmntType, 'String', num2str(type))
    if type == 0
        set(handles.txt_vslType, 'String', 'unspecified')
    elseif type == 1
        set(handles.txt_vslType, 'String', 'Artery')
    elseif type == 2
        set(handles.txt_vslType, 'String', 'Capillary')
    elseif type == 3
        set(handles.txt_vslType, 'String', 'Vein')
    end
    
else
    set(handles.edgeID, 'string', 'Ø')
    set(handles.segmentID, 'string', 'Ø')
    set(handles.cylDiam, 'String', 'Ø');
    set(handles.SgmntDiam_Est, 'String', 'Ø')
    set(handles.edgeDiam_Est, 'string', 'Ø')
    set(handles.sgmntType, 'String', 'Ø')
    set(handles.txt_vslType, 'String', '    Ø')
    set(handles.cutOff, 'String', 'Ø')
    set(handles.brOrdr_art, 'string', 'Ø')
    set(handles.brOrdr_vein, 'string', 'Ø')
end
end


% --------------------------------------------------------------------
function cntxtMenuParent_mntrMode_Callback(hObject, eventdata, handles)
% hObject    handle to cntxtMenuParent_mntrMode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
end

% --------------------------------------------------------------------
function cntxtMenu_inspect_Callback(hObject, ~, handles)
global im0
hObject.Checked = 'on';
handles.cntxtMenu_cylinder.Checked = 'off';
handles.cntxtMenu_mntrsOff.Checked = 'off';
set(handles.inspectMode, 'Value', 1)
update_Monitors(handles, im0.gui.edgeID, true)
% hObject    handle to cntxtMenu_inspect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
end

% --------------------------------------------------------------------
function cntxtMenu_cylinder_Callback(hObject, ~, handles)
global im0
hObject.Checked = 'on';
handles.cntxtMenu_inspect.Checked = 'off';
handles.cntxtMenu_mntrsOff.Checked = 'off';
set(handles.cylMode, 'Value', 1)
update_Monitors(handles, im0.gui.edgeID, true)
% hObject    handle to cntxtMenu_cylinder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
end

% --------------------------------------------------------------------
function cntxtMenu_mntrsOff_Callback(hObject, ~, handles)
global im0
hObject.Checked = 'on';
handles.cntxtMenu_inspect.Checked = 'off';
handles.cntxtMenu_cylinder.Checked = 'off';
set(handles.inspectMode, 'Value', 0)
set(handles.cylMode, 'Value', 0)

%%-- window_01
imagesc(handles.axes_OCT_Win_01, zeros(im0.gui.winSize),...
    [im0.gui.colorbalLowerLim im0.gui.colorbarUpperLim])
colormap(handles.axes_OCT_Win_01, im0.gui.flowClrMp);
xticks(handles.axes_OCT_Win_01, [])
yticks(handles.axes_OCT_Win_01, [])
handles.axes_OCT_Win_01.Color = 'k';

%%-- window_02
imagesc(handles.axes_OCT_Win_02, zeros(im0.gui.winSize),...
    [im0.gui.colorbalLowerLim im0.gui.colorbarUpperLim])
colormap(handles.axes_OCT_Win_02, im0.gui.flowClrMp);
xticks(handles.axes_OCT_Win_02, [])
yticks(handles.axes_OCT_Win_02, [])
handles.axes_OCT_Win_02.Color = 'k';

%%-- window_03
imagesc(handles.axes_OCT_Win_03, zeros(im0.gui.winSize),...
    [im0.gui.colorbalLowerLim im0.gui.colorbarUpperLim])
colormap(handles.axes_OCT_Win_03, im0.gui.flowClrMp);
xticks(handles.axes_OCT_Win_03, [])
yticks(handles.axes_OCT_Win_03, [])
handles.axes_OCT_Win_03.Color = 'k';

%%-- wind_TPAng
imagesc(handles.axes_TPAng_win, zeros(im0.gui.winSize),...
    [im0.gui.TPthreshold*32 32])
colormap(handles.axes_TPAng_win, gray);
xticks(handles.axes_TPAng_win, [])
yticks(handles.axes_TPAng_win, [])
handles.axes_TPAng_win.Color = 'k';
% hObject    handle to cntxtMenu_mntrsOff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
end


% --- Executes on button press in cylMode.
function cylMode_Callback(~, ~, handles)
global im0
handles.cntxtMenu_inspect.Checked = 'off';
handles.cntxtMenu_cylinder.Checked = 'on';
handles.cntxtMenu_mntrsOff.Checked = 'off';
update_Monitors(handles, im0.gui.edgeID, true)
% hObject    handle to cylMode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of cylMode
end

% --------------------------------------------------------------------
function toolBar_save_ClickedCallback(~, ~, ~)  
global im0

if isfield(im0, 'fName')
    answer = questdlg('Are you sure you want to overwrite the existing file?', ...
        'Save', ...
        'Yes','No','Cancel','No');
    
    switch answer
        case 'Yes'
            hwait = waitbar(0.2, 'Saving Session');
            save(im0.fName, 'im0')
            waitbar(0.95,hwait,'Finishing');
            pause(0.25)
            waitbar(1,hwait,'Done!');
            pause(0.25)
            close(hwait)
            
        case 'No'
           menu_saveAs_Callback()
            
        case 'Cancel'
            errordlg('File was not saved!', 'alert')
            return
    end
else
    menu_saveAs_Callback()
end
end
% hObject    handle to toolBar_save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menu_save_Callback(~, ~, ~)
toolBar_save_ClickedCallback()
% hObject    handle to menu_save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
end


% --------------------------------------------------------------------
function menu_saveAs_Callback(~, ~, ~)
global im0

if isfield(im0, 'fName')
    [fP, fN, ex] = fileparts(im0.fName);
    savePath = [fP, '/', fN, '_01', ex];
    [f, p] = uiputfile(savePath);
    if f == 0
        errordlg('File was not saved!', 'alert')
        return
    end
    
    im0.fName = fullfile(p, f);
    hwait = waitbar(0.2,'Saving Session');
    save(im0.fName, 'im0')
    waitbar(0.95,hwait,'Finishing');
    pause(0.25)
    waitbar(1,hwait,'Done!');
    pause(0.25)
    close(hwait)
    
else
    [f, p] = uiputfile('Session001.mat');
    if f == 0
        errordlg('File was not saved!', 'alert')
        return
    end
    im0.fName = fullfile(p, f);
    
    hwait = waitbar(0.2,'Saving Session');
    save(im0.fName, 'im0')
    waitbar(0.95,hwait, 'Finishing');
    pause(0.25)
    waitbar(1, hwait, 'Done!');
    pause(0.25)
    close(hwait)
    
end

% hObject    handle to menu_saveAs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
end


% --------------------------------------------------------------------
function menu_loadSession_Callback(~, ~, handles)
openingDlg(handles)
% hObject    handle to menu_loadSession (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
end


% --------------------------------------------------------------------
function menu_refreshLst_Callback(~, ~, handles)
updateGui_Listbox(handles)
% hObject    handle to menu_refreshLst (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
end


% --------------------------------------------------------------------
function cntxMenu_refreshLst_Callback(~, ~, handles)
updateGui_Listbox(handles)
% hObject    handle to cntxMenu_refreshLst (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
end

% --------------------------------------------------------------------
function cntxtMenu_lst_Callback(hObject, eventdata, handles)
% hObject    handle to cntxtMenu_lst (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
end

% --------------------------------------------------------------------
function cntxtMenu_lst_nex_Callback(hObject, eventdata, handles)
handles.cntxtMenu_lst_nex.Checked = 'on';
handles.cntxtMenu_lst_prev.Checked = 'off';
% hObject    handle to cntxtMenu_lst_nex (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
end

% --------------------------------------------------------------------
function cntxtMenu_lst_prev_Callback(hObject, eventdata, handles)
handles.cntxtMenu_lst_nex.Checked = 'off';
handles.cntxtMenu_lst_prev.Checked = 'on';
% hObject    handle to cntxtMenu_lst_prev (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
end

% --------------------------------------------------------------------
function highlightEdges(handles, edgeID)
global im0
if ~isempty(edgeID)
    switch im0.gui.selectionType
        case 'selectOneEdge'
            if ~isempty(im0.gui.edgeFlagTemp) && isempty(intersect(edgeID, im0.gui.edgeFlagTemp(1,:)))
                im0.edgeFlag(im0.gui.edgeFlagTemp(1,:)) = im0.gui.edgeFlagTemp(2,:);

                im0.gui.edgeFlagTemp = [edgeID; im0.edgeFlag(edgeID)];
                im0.edgeFlag(edgeID) = 2;

            elseif ~isempty(im0.gui.edgeFlagTemp) && ~isempty(intersect(edgeID, im0.gui.edgeFlagTemp(1,:)))
                lgc = ismember(im0.gui.edgeFlagTemp(1,:), edgeID);
                im0.edgeFlag(im0.gui.edgeFlagTemp(1,lgc)) = im0.gui.edgeFlagTemp(2, lgc);
            else
                im0.gui.edgeFlagTemp = [edgeID; im0.edgeFlag(edgeID)];
                im0.edgeFlag(edgeID) = 2;
            end
            
            
        case 'selectOneSgmnt'
            if ~isempty(im0.gui.edgeFlagTemp) 
                im0.edgeFlag(im0.gui.edgeFlagTemp(1,:)) = im0.gui.edgeFlagTemp(2,:);
            end
            sg = im0.edgeSegN(edgeID);
            edges = find(im0.edgeSegN == sg)';

            if ~isempty(im0.gui.edgeFlagTemp) && ~isempty(intersect(edges, im0.gui.edgeFlagTemp(1,:)))
                im0.edgeFlag(im0.gui.edgeFlagTemp(1,:)) = im0.gui.edgeFlagTemp(2,:);
                im0.gui.edgeFlagTemp = [];
            else
                flg = im0.edgeFlag(edges);  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Qi 5-May-2023
                if size(flg, 2) == 1
                    flg = flg';
                end

                im0.gui.edgeFlagTemp = [edges; flg];
                im0.edgeFlag(edges) = 2;
            end
            
            
        case 'selectMultEdges'
            if ~isempty(im0.gui.edgeFlagTemp) && any( intersect(edgeID, im0.gui.edgeFlagTemp(1,:)) )
                im0.edgeFlag(edgeID) = im0.gui.edgeFlagTemp(2, im0.gui.edgeFlagTemp(1,:) == edgeID);
            else
                im0.gui.edgeFlagTemp(:, end+1) = [edgeID;  im0.edgeFlag(edgeID)'];
                im0.edgeFlag(edgeID) = 2;
            end

        case 'selectMultSgmnts'
            sg = im0.edgeSegN(edgeID);
            edges = find(im0.edgeSegN == sg);
            if ~isempty(im0.gui.edgeFlagTemp) && intersect(edges, im0.gui.edgeFlagTemp(1,:))
                lgc = ismember(im0.gui.edgeFlagTemp(1,:), edges);
                im0.edgeFlag(im0.gui.edgeFlagTemp(1,lgc)) = im0.gui.edgeFlagTemp(2, lgc);
            else
                im0.gui.edgeFlagTemp = [im0.gui.edgeFlagTemp, [edges;  im0.edgeFlag(edges)']];
                im0.edgeFlag(edges) = 2;
            end

            im0.edgeFlag(im0.edgeSegN == sg) = 2;


    end
end

if all(im0.gui.positionFlag) && ~isempty(edgeID)
    update_Monitors(handles, edgeID, true)
    update_TPMA_image(handles)
    update_OCTV_image(handles)
    positionFlag_OCTFlow(handles, im0.gui.mntrPt)
    positionFlag_TPAng(handles, im0.gui.mntrPt)
    overlayGraph(handles)
else
    update_Monitors(handles, edgeID, true)
    update_TPMA_image(handles)
    update_OCTV_image(handles)
    overlayGraph(handles)
end

end

% --------------------------------------------------------------------
function menu_loadPial_Callback(~, ~, ~)
%%
global im0
[fN, fP] = uigetfile('*.txt','Select Pial Vessels TXT File');
if fN == 0; return; end
array = table2array(readtable([fP fN]));
im0.VsslLbeling.pial = array;

im0.segVesType(array(:, 1)') = array(:, 3)';

msgbox('file loaded');
% hObject    handle to menu_loadPial (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
end

% --------------------------------------------------------------------
function menu_loadArtAndVen_Callback(~, ~, ~)
global im0
[fN, fP] = uigetfile('*.txt','Art. and Ven. Vessels TXT File');
if fN == 0; return; end
array = table2array(readtable([fP fN]));
im0.VsslLbeling.artsAndVens = array;
im0.segVesType(array(:, 1)') = array(:, 3)';

msgbox('file loaded');
% hObject    handle to menu_loadArtAndVen (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
end

%---------------------------------------------------------------------
function calcBranchOrdr()

global im0
%[fname pname] = uigetfile('*.seed','Select seed file with Im structure');
%load([pname fname],'-mat'); % load 'im2' structure

%%

choice = questdlg('Do you want to recalculate segment distance matrices?',...
    'Long or quick calculation?','Yes, recalculate all',...
    'No, just load pial vessels','No, just load pial vessels');
if strcmp(choice, 'Yes, recalculate all') % recalculate matrices...
    
    nGroups = length(unique(im0.nodeGrp));
    disp(['Number of groups is ' num2str(nGroups) ]);
    disp('Grp. number      num. elements');
    for i=1:nGroups
        grpElements(i) = length(find(im0.nodeGrp==i)); %#ok<AGROW> 
        disp([num2str(i) '         ' num2str( grpElements(i) )  ]);
    end
    [mg, mi] = max(grpElements); %#ok<ASGLU> 
    answer = inputdlg(['Select group 1 to ' num2str(nGroups)],'Select which group to process',1,{num2str(mi(1))});
    grp = str2double(answer{1});
    im0.grpStat.selectedGroupNumber = grp;
    % find list of segments that belong to selected group
    totNumSegments = length(im0.segDiam); %#ok<NASGU> % number of segments
    foo1 = im0.nodeGrp(im0.segEndNodes(:,1));
    foo2 = im0.nodeGrp(im0.segEndNodes(:,2));
    segLst = find(foo1==grp | foo2==grp); % indexes of segments which are in group grp
    numSegments = length(segLst);
    im0.grpStat.segLst = segLst;
    im0.grpStat.numSegments = numSegments;

    %%
    % start from largest artery (vessel types: A=1, C=2, V=3)
    %idxArt = find(im0.segVesType==1);
    %[foo, maxArtIdx] = sort(im0.segDiam(idxArt),'descend');
    %maxArtIdx = idxArt(maxArtIdx); % descending indexes of arterial segments
    %%
    segEndNodes = im0.segEndNodes(segLst,:); % use only segments that belong to group grp
    endNodes = unique(segEndNodes(:));
    endNodes = sort(endNodes);
    numEndNodes = length(endNodes);
    lookuptable(1:numEndNodes,1) = endNodes;
    im0.grpStat.lookuptable = lookuptable;
    im0.grpStat.numEndNodes = numEndNodes;

    %%
    %         % create connectivity matrix
    %         segConnectivityMartix = zeros(numEndNodes);
    %         %segConnectivityMartix(numSegments,numSegments)=1;
    %         hp = waitbar(0,'Populating first order segment connectivity matrix'); 
    %         for ii=1:numEndNodes-1,
    %             waitbar(ii/(numEndNodes-1),hp);
    %             node1 = lookuptable(ii);
    %             %%segConnectivityMartix(ii,ii)=1;
    %             %%disp(num2str(ii));
    %             for jj=ii+1:numEndNodes,
    %                 node2 = lookuptable(jj);
    %                 lst = find(im0.segEndNodes(:,1) == node1);
    %                 lst2 = find(im0.segEndNodes(:,2) == node2);
    %                 c = intersect(lst, lst2);
    %                 if ~isempty(c),
    %                     segConnectivityMartix(ii,jj)=im0.segLen_um(c(1));
    %                     segConnectivityMartix(jj,ii)=segConnectivityMartix(ii,jj);
    %                 else
    %                     lst = find(im0.segEndNodes(:,1) == node2);
    %                     lst2 = find(im0.segEndNodes(:,2) == node1);
    %                     c = intersect(lst, lst2);
    %                     if ~isempty(c),
    %                         segConnectivityMartix(ii,jj)=im0.segLen_um(c(1));
    %                         segConnectivityMartix(jj,ii)=segConnectivityMartix(ii,jj);
    %                     end;
    %                 end;
    %             end;
    %         end;
    %         close(hp);
    %         segConnectivityMartix_1stOrder = segConnectivityMartix;

    %%
    % create matrix E(1:numSegments,1:3) with elements (node1, node2, segmentLEngth)
    segLen_um = im0.segLen_um(segLst);
    E = zeros(2*numSegments,3); % 2 times for two directions of arrow
    for ii=1:numSegments
        nodes = segEndNodes(ii,:);
        E(2*ii-1,1)=find(lookuptable == nodes(1));
        E(2*ii-1,2)=find(lookuptable == nodes(2));
        E(2*ii-1,3)=segLen_um(ii);

        E(2*ii,1)=E(2*ii-1,2);
        E(2*ii,2)=E(2*ii-1,1);
        E(2*ii,3)=E(2*ii-1,3);
    end

    %%
    % calculate shortest distance between nodes
    m=2*numSegments; % number of oriented edges
    n = numEndNodes;
    %[m,n,E] = grValidation(E); % E data validation

    % ================ Initial values ===============
    dSP=ones(n)*inf; % initial distances
    dSP((E(:,2)-1)*n+E(:,1))=E(:,3);
    %dSP0=dSP;
    % ========= The main cycle of Floyd-Warshall algorithm =========
    hp = waitbar(0,'Calculating minimum distance for all graph nodes...'); 
    for j=1:n
      waitbar(j/n,hp);
      i=setdiff((1:n),j);
      dSP(i,i)=min(dSP(i,i),repmat(dSP(i,j),1,n-1)+repmat(dSP(j,i),n-1,1));
    end
    close(hp);

    %%
    % find shortest path from node1 ('s') to node2 ('t')

    % populate matrix nSP (numEndNodes x numEndNodes), where nSP(i,j) is index of node
    % before 'j' in shortest path from i to j

    nSP = zeros(numEndNodes);
    dSP1=dSP;
    dSP1(1:n+1:n^2)=0; % modified dSP

    hp = waitbar(0,'Calculating shortest paths for all pairs of graph nodes...'); 
    for ii=1:numEndNodes-1
        waitbar(ii/(numEndNodes-1),hp);
        s=ii;
        %s=1;
        for jj=ii+1:numEndNodes
            t=jj;
            %t=915;

            sp=[];
            %s=s(1);
            %t=t(1);
            if isinf(dSP(s,t)) % t is not accessible from s
                disp('Grrrrr! This pair of nodes is not connected');
            else

                l=ones(m,1); % label for each arrow
                sp=t; % final vertex
                while ~(sp(1)==s)
                  nv=find((E(:,2)==sp(1))&l); % nv contains all labeled arrows (segments) pointing at node sp(1)
                  vnv=abs((dSP1(s,sp(1))-dSP1(s,E(nv,1)))'-E(nv,3))<eps*1e4; % valided arrows (I changed 1e3 to 1e4)
                  l(nv(~vnv))=0; % labels of not valided arrows
                  if all(~vnv) % invalided arrows
                    l(find((E(:,1)==sp(1))&(E(:,2)==sp(2))))=0; 
                    sp=sp(2:end); % one step back
                  else
                    nv=nv(vnv); % rested valided arrows
                    sp=[E(nv(1),1) sp]; % add one vertex to shortest path
                  end
                end
            end
            nSP(ii,jj) = sp(end-1);
            nSP(jj,ii) = sp(2);

        end
    end

    for ii=1:numEndNodes % make sure that distances on diagonal are zero
        dSP(ii,ii)=0;
    end
    im0.grpStat.dSP = dSP;
    im0.grpStat.nSP = nSP;

    close(hp);

    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Above was calculation based on distances. Now for branching order
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%
    % create matrix E(1:numSegments,1:3) with elements (node1, node2, segmentLEngth)

    %Ebr = zeros(2*numSegments,3); % 2 times for two directions of arrow
    Ebr = E;
    Ebr(:,3) = 1; % this will assure that we calculate shortest paths based on branching order

    %%
    % calculate shortest distance between nodes
    m=2*numSegments; % number of oriented edges
    n = numEndNodes;
    %[m,n,E] = grValidation(E); % E data validation

    % ================ Initial values ===============
    dSPbr=ones(n)*inf; % initial distances
    dSPbr((Ebr(:,2)-1)*n+Ebr(:,1))=Ebr(:,3);
    %dSP0=dSP;
    % ========= The main cycle of Floyd-Warshall algorithm =========
    hp = waitbar(0,'Calculating minimum distance for all graph nodes...'); 
    for j=1:n
      waitbar(j/n,hp);
      i=setdiff((1:n),j);
      dSPbr(i,i)=min(dSPbr(i,i),repmat(dSPbr(i,j),1,n-1)+repmat(dSPbr(j,i),n-1,1));
    end
    close(hp);

    %%
    % find shortest path from node1 ('s') to node2 ('t')

    % populate matrix nSP (numEndNodes x numEndNodes), where nSP(i,j) is index of node
    % before 'j' in shortest path from i to j

    nSPbr = zeros(numEndNodes);
    dSP1br=dSPbr;
    dSP1br(1:n+1:n^2)=0; % modified dSP

    hp = waitbar(0,'Calculating shortest paths for all pairs of graph nodes...'); 
    for ii=1:numEndNodes-1
        waitbar(ii/(numEndNodes-1),hp);
        s=ii;
        %s=1;
        for jj=ii+1:numEndNodes
            t=jj;
            %t=915;

            sp=[];
            %s=s(1);
            %t=t(1);
            if isinf(dSPbr(s,t)) % t is not accessible from s
                disp('Grrrrr! This pair of nodes is not connected');
            else

                l=ones(m,1); % label for each arrow
                sp=t; % final vertex
                while ~(sp(1)==s)
                  nv=find((Ebr(:,2)==sp(1))&l); % nv contains all labeled arrows (segments) pointing at node sp(1)
                  vnv=abs((dSP1br(s,sp(1))-dSP1br(s,Ebr(nv,1)))'-Ebr(nv,3))<eps*1e4; % valided arrows (I changed 1e3 to 1e4)
                  l(nv(~vnv))=0; % labels of not valided arrows
                  if all(~vnv) % invalided arrows
                    l(find((Ebr(:,1)==sp(1))&(Ebr(:,2)==sp(2))))=0; 
                    sp=sp(2:end); % one step back
                  else
                    nv=nv(vnv); % rested valided arrows
                    sp=[Ebr(nv(1),1) sp]; % add one vertex to shortest path
                  end
                end
            end
            nSPbr(ii,jj) = sp(end-1);
            nSPbr(jj,ii) = sp(2);

        end
    end

    for ii=1:numEndNodes % make sure that distances on diagonal are zero
        dSPbr(ii,ii)=0;
    end
    im0.grpStat.dSPbr = dSPbr;
    im0.grpStat.nSPbr = nSPbr;

    close(hp);

end % recalculate matrices...



% NOW, IMPORT NEW PIAL VESSEL BRANCH ORDER AND CALCULATE SEGMENT DISTANCES
% TO THESE PIAL SEGMENTS...

%[fname pname] = uigetfile('*.seed','Select seed file with Im structure');
%load([pname fname],'-mat'); % load 'im2' structure
%%
assignBranchOrdr()
end


function assignBranchOrdr()
global im0
%%
if ~isfield(im0,'grpStat')
    disp('No info on group statistics!');
    disp('Run Calc. Segment Distances first!');
    return;
end
%%
grp = im0.grpStat.selectedGroupNumber; 
segLst = im0.grpStat.segLst;
numSegments = im0.grpStat.numSegments;
lookuptable = im0.grpStat.lookuptable;
numEndNodes = im0.grpStat.numEndNodes;
dSP = im0.grpStat.dSP;
nSP = im0.grpStat.nSP;
dSPbr = im0.grpStat.dSPbr;
nSPbr = im0.grpStat.nSPbr;



%%

% First, supply segment numbers in pial arteries and veins together
% with branch order. Algorithm will calculate branching order for each 
% vessel segment in selected connected group of segments


%Bill Sep 10 2015. avoid making the text file with pial segment seeds by
%having these segments labeled in the GUI and saved in im0.grpStat.pialseg
% only 3 col.

if isfield(im0.grpStat, 'pialseg')
pialseg = im0.grpStat.pialseg;
else
% load file with info about pial segments
[ff pp] = uigetfile('*.txt','Select File name with Pial Segments info');
pialseg = load([pp ff],'-ascii'); % 3 columns, segment number, branch order, vessel type (A=1, 2=C, 3=V)
end
%%
indxDlt = [];
for ii=1:size(pialseg,1) % find end nodes of pial segments and convert them based on lookuptable
    lgc = ~isempty(find(lookuptable == im0.segEndNodes(pialseg(ii,1),1), 1)) &&...
        ~isempty(find(lookuptable == im0.segEndNodes(pialseg(ii,1),2), 1));
    if lgc
        pialseg(ii,4) = find(lookuptable == im0.segEndNodes(pialseg(ii,1),1));
        pialseg(ii,5) = find(lookuptable == im0.segEndNodes(pialseg(ii,1),2));
    else
        indxDlt = [indxDlt, ii];
    end
end
%%
pialseg(indxDlt, :) = [];
im0.grpStat.pialseg = pialseg;



% first calculate branch order for each segment in selected graph group -
% this is not neccessary for PO2, but it is important for testing the
% algorithm
segBranchOrder = zeros(numSegments,6); 
% columns are: true segment number, branch order from A, length from A, BR
% order from V, length from V, diameter in 'um'

% BRANCHING ORDER FROM PIAL ARTERIES
maxBR = max(dSPbr(:));
for ii=1:numSegments
    segBranchOrder(ii,1) = segLst(ii);
    segBranchOrder(ii,6) = im0.segDiam(segLst(ii)); % segment diameters
    segnode1 = find(lookuptable == im0.segEndNodes(segLst(ii),1));
    segnode2 = find(lookuptable == im0.segEndNodes(segLst(ii),2));
    % find minimum branch distance to each node of pial arterial segments
    lstArt = find(pialseg(:,3)==1);
    minBR = maxBR+1;
    sameAsPial = 0;
    for jj = 1:length(lstArt)
        Anode1 = pialseg(lstArt(jj),4);
        Anode2 = pialseg(lstArt(jj),5);
        if ( ((Anode1 == segnode1) && (Anode2 == segnode2)) || ((Anode1 == segnode2) && (Anode2 == segnode1))  ),
            segBranchOrder(ii,2) = pialseg(lstArt(jj),2); % same branch order as pial artery
            %segBranchOrder(ii,3) = 0; % length from pial artery is 0 -  not good, it should look for '0' order arteries
            sameAsPial = 1;
            break;
        else            
            foo1 = dSPbr(segnode1,Anode1);
            foo2 = dSPbr(segnode1,Anode2);
            foo3 = dSPbr(segnode1,Anode1);
            foo4 = dSPbr(segnode1,Anode2);
            foo = min([foo1 foo2 foo3 foo4]);
            if foo<minBR
                minBR = foo;
                segBRord = pialseg(lstArt(jj),2);
            end
        end
    end
    if ~sameAsPial
        segBranchOrder(ii,2) = minBR+segBRord+1; 
    end
end


% DISTANCE FROM SURFACE PIAL ARTERIES WITH BRANCHING ORDER ZERO
maxDIS = 1000000; % some large number of micrones
for ii=1:numSegments
    segnode1 = find(lookuptable == im0.segEndNodes(segLst(ii),1));
    segnode2 = find(lookuptable == im0.segEndNodes(segLst(ii),2));
    % find minimum distance to each node of pial arterial segments
    lstArt = find( (pialseg(:,3)==1) & (pialseg(:,2) == 0) ); % find arterial pial vessels with branching order = 0
    minDIS = maxDIS+1;
    sameAsPial = 0;
    for jj = 1:length(lstArt)
        Anode1 = pialseg(lstArt(jj),4);
        Anode2 = pialseg(lstArt(jj),5);
        if ( ((Anode1 == segnode1) && (Anode2 == segnode2)) || ((Anode1 == segnode2) && (Anode2 == segnode1))  ),
            segBranchOrder(ii,3) = 0; % segment equal to pial atery with BR order zero, so length is 0
            sameAsPial = 1;
            break;
        else            
            foo1 = dSP(segnode1,Anode1);
            foo2 = dSP(segnode1,Anode2);
            foo3 = dSP(segnode1,Anode1);
            foo4 = dSP(segnode1,Anode2);
            foo = min([foo1 foo2 foo3 foo4]);
            if foo<minDIS
                minDIS = foo;
                segDIS = dSP(Anode1, Anode2);
            end
        end
    end
    if ~sameAsPial
        segBranchOrder(ii,3) = minDIS+segDIS/2.0+dSP(segnode1,segnode2)/2.0; % add 1/2 of lengths of segments which are compared
    end
end


% BRANCHING ORDER FROM PIAL VEINS
maxBR = max(dSPbr(:));
for ii=1:numSegments
    segnode1 = find(lookuptable == im0.segEndNodes(segLst(ii),1));
    segnode2 = find(lookuptable == im0.segEndNodes(segLst(ii),2));
    % find minimum branch distance to each node of pial arterial segments
    lstVein = find(pialseg(:,3)==3);
    minBR = maxBR+1;
    sameAsPial = 0;
    for jj = 1:length(lstVein)
        Vnode1 = pialseg(lstVein(jj),4);
        Vnode2 = pialseg(lstVein(jj),5);
        if ( ((Vnode1 == segnode1) && (Vnode2 == segnode2)) || ((Vnode1 == segnode2) && (Vnode2 == segnode1))  ),
            segBranchOrder(ii,4) = pialseg(lstVein(jj),2); % same branch order as pial vein
            sameAsPial = 1;
            break;
        else            
            foo1 = dSPbr(segnode1,Vnode1);
            foo2 = dSPbr(segnode1,Vnode2);
            foo3 = dSPbr(segnode1,Vnode1);
            foo4 = dSPbr(segnode1,Vnode2);
            foo = min([foo1 foo2 foo3 foo4]);
            if foo<minBR
                minBR = foo;
                segBRord = pialseg(lstVein(jj),2);
            end
        end
    end
    if ~sameAsPial
        segBranchOrder(ii,4) = minBR+segBRord+1; 
    end
end


% DISTANCE FROM SURFACE PIAL VEINS WITH BRANCHING ORDER ZERO

% this needs improvement! For example, there may not be zero branches of
% veins in FOV. It should look to min distance to given vessels in pialseg,
% but then it should know distances of vessels inside pialseg from surface
% pial vessels!

maxDIS = 1000000; % some large number of micrones
for ii=1:numSegments
    segnode1 = find(lookuptable == im0.segEndNodes(segLst(ii),1));
    segnode2 = find(lookuptable == im0.segEndNodes(segLst(ii),2));
    % find minimum distance to each node of pial venous segments
    lstVein = find( (pialseg(:,3)==3) & (pialseg(:,2) == 0) ); % find venous pial vessels with branching order = 0
    minDIS = maxDIS+1;
    sameAsPial = 0;
    for jj = 1:length(lstVein)
        Vnode1 = pialseg(lstVein(jj),4);
        Vnode2 = pialseg(lstVein(jj),5);
        if ( ((Vnode1 == segnode1) && (Vnode2 == segnode2)) || ((Vnode1 == segnode2) && (Vnode2 == segnode1))  ),
            segBranchOrder(ii,5) = 0; % segment equal to pial atery with BR order zero, so length is 0
            sameAsPial = 1;
            break;
        else            
            foo1 = dSP(segnode1,Vnode1);
            foo2 = dSP(segnode1,Vnode2);
            foo3 = dSP(segnode1,Vnode1);
            foo4 = dSP(segnode1,Vnode2);
            foo = min([foo1 foo2 foo3 foo4]);
            if foo<minDIS
                minDIS = foo;
                segDIS = dSP(Vnode1, Vnode2);
            end
        end
    end
    if ~sameAsPial
        segBranchOrder(ii,5) = minDIS+segDIS/2.0+dSP(segnode1,segnode2)/2.0; % add 1/2 of lengths of segments which are compared
    end
end


im0.grpStat.segBranchOrder = segBranchOrder;

% assign branching orders to all edges in the graph...
for ee=1:length(im0.edgeSegN)
    [Lia, Locb] = ismember(im0.edgeSegN(ee),im0.grpStat.segLst);
    if Lia % find if edge ee belongs to selected group
        im0.edgeBRorderVeins(ee) = im0.grpStat.segBranchOrder(Locb,4);
        im0.edgeBRorderArt(ee) = im0.grpStat.segBranchOrder(Locb,2);
    else % edge outside of desired group
        im0.edgeBRorderVeins(ee) = 0;
        im0.edgeBRorderArt(ee) = 0; % I may need better number for edges outside the group. This will make them branch order zero
    end
    %im0.edgeBRorderVeins(ee) = segBranchOrder(im0.edgeSegN(ee),4);
    %im0.edgeBRorderArt(ee) = segBranchOrder(im0.edgeSegN(ee),2);
end


sgmntIDs = 1: length(im0.segDiam);
im0.sgmntBrOrderArt = zeros(1, max(sgmntIDs(:)));
im0.sgmntBrOrderVein = zeros(1, max(sgmntIDs(:)));

for sg = sgmntIDs
    brOr = im0.edgeBRorderArt(im0.edgeSegN == sg);
    brOr = unique(brOr);
    im0.sgmntBrOrderArt(sg) = brOr;
    
    brOr = im0.edgeBRorderVeins(im0.edgeSegN == sg);
    brOr = unique(brOr);
    im0.sgmntBrOrderVein(sg) = brOr;
end
end

% --------------------------------------------------------------------
function menu_calcBrOrdr_Callback(hObject, eventdata, handles)
calcBranchOrdr()
msgbox('process complete')
% hObject    handle to menu_calcBrOrdr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
end

% --------------------------------------------------------------------
function menu_reassignBrOrdr_Callback(hObject, eventdata, handles)
global im0
% prompt = {'max capillary diameter (μm):'};
% dlgtitle = 'Input diameter limit';
% dims = [1 30];
% definput = {'8'};
% answer = cell2mat(inputdlg(prompt, dlgtitle, dims, definput));
% lmt = str2double(answer);
% 
% if isempty(lmt) || lmt < 0
%     errordlg({'process was aborted', 'due to improper diameter limit'}, 'input error')
%     return
% end



%--
% im0.grpStat.pialseg = im0.VsslLbeling.artsAndVens;
% im0.grpStat.pialseg = im0.VsslLbeling.pial;
% assignBranchOrdr()
%--

%- comented on 10 aug 2020
% labelSegments(handles, lmt)

% im0.segVesType(:) = 2;
% im0.segVesType(im0.VsslLbeling.pial(:,1)) = im0.VsslLbeling.pial(:,3);
% im0.segVesType(im0.VsslLbeling.artsAndVens(:,1)) = im0.VsslLbeling.artsAndVens(:,3);

sgmnts = find(im0.segVesType ~= 2)';
im0.grpStat.pialseg = [sgmnts, zeros(length(sgmnts), 1), im0.segVesType(sgmnts)'];
assignBranchOrdr()

%%-------------------------------------------------------------------------
% sgmnts = im0.selectedEdgesList(1, :);
% types = im0.selectedEdgesList(4, :);
% flows = im0.selectedEdgesList(3, :);
% artSg = sgmnts(types ~= 2 & flows < 0);
% veinSg = sgmnts(types ~= 2 & flows > 0);
% im0.segVesType(artSg) = 1;
% im0.segVesType(veinSg) = 3;
% 
% updateSelectedEdgesList(handles)
% 
% if ~isempty(im0.VsslLbeling.pial)
%     for sg = artSg
%         im0.VsslLbeling.pial(im0.VsslLbeling.pial(:,1) == sg, 3) = 1;
%     end
%     for sg = veinSg
%         im0.VsslLbeling.pial(im0.VsslLbeling.pial(:,1) == sg, 3) = 3;
%     end
% end
% if ~isempty(im0.VsslLbeling.artsAndVens)
%     for sg = artSg
%         im0.VsslLbeling.artsAndVens(im0.VsslLbeling.artsAndVens(:,1) == sg, 3) = 1;
%     end
%     for sg = veinSg
%         im0.VsslLbeling.artsAndVens(im0.VsslLbeling.artsAndVens(:,1) == sg, 3) = 3;
%     end
% end
%%-------------------------------------------------------------------------


%-------------14-may-2020
% sgmnts_I = find(im0.segDiam > lmt);
% sgmnts_II = im0.VsslLbeling.pial(:, 1)';
% sgmnts_III = im0.VsslLbeling.artsAndVens(:, 1)';
% sgmnts = unique([sgmnts_I, sgmnts_II, sgmnts_III]);
% types = im0.segVesType(sgmnts);
% order = zeros(length(sgmnts), 1);
% im0.grpStat.pialseg = [sgmnts', order, types'];
% 
% assignBranchOrdr()
% labelSegments(handles, lmt)
% msgbox('Operation Completed')


% hObject    handle to menu_reassignBrOrdr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

msgbox('Operation Completed')
end


function labelSegments(handles, thrshld)  
%%
global im0

im0.segVesType = ones(1,size(im0.segDiam, 2)) * 2;
im0.segVesType(im0.VsslLbeling.pial(:,1)') = im0.VsslLbeling.pial(:,3)';
im0.segVesType(im0.VsslLbeling.artsAndVens(:,1)') = im0.VsslLbeling.artsAndVens(:,3)';

%-- indx of capillaries
% CapInd = im0.segDiam <= lmt;
% im0.segVesType(CapInd) = 2;
% 
% remainingSgmnts = find(im0.segVesType == 0);
% 
% for sg = remainingSgmnts
%     ArtBrOr = im0.sgmntBrOrderArt(sg);
%     VBrOr = im0.sgmntBrOrderVein(sg);
%     if ArtBrOr < VBrOr
%         im0.segVesType(sg) = 1;
%     elseif ArtBrOr > VBrOr
%         im0.segVesType(sg) = 3;
%     end
% end

%---temp addition
% im0.segVesType(im0.segVesType == 0) = 2;
%---


labeledSegs = unique([im0.VsslLbeling.pial(:,1)', im0.VsslLbeling.artsAndVens(:,1)'])';
l = length(labeledSegs);
i = 1;

wtBr = waitbar(0, 'labeling');

while i <= l
    %%
    sg = labeledSegs(i);
    type = im0.segVesType(sg);
    waitbar(i/l, wtBr, ['labeling - sg ' num2str(sg)])
    %%
    endNodes = im0.segEndNodes(sg, :);
    cnctdSgs = find( (im0.segEndNodes(:, 1) == endNodes(1) |...
        im0.segEndNodes(:, 1) == endNodes(2) |...
        im0.segEndNodes(:, 2) == endNodes(1) |...
        im0.segEndNodes(:, 2) == endNodes(2)) &...
        ~(im0.segEndNodes(:, 1) == endNodes(1)&...
        im0.segEndNodes(:, 2) == endNodes(2)));
    %%
    alredyLbld_indx = [];
    for ii = 1:length(cnctdSgs)
        
        %%%
        endNodes_subSeg = im0.segEndNodes(cnctdSgs(ii), :);
        cnctdSubSgs = find( (im0.segEndNodes(:, 1) == endNodes_subSeg(1) |...
            im0.segEndNodes(:, 1) == endNodes_subSeg(2) |...
            im0.segEndNodes(:, 2) == endNodes_subSeg(1) |...
            im0.segEndNodes(:, 2) == endNodes_subSeg(2)) &...
            ~((im0.segEndNodes(:, 1) == endNodes_subSeg(1)&...
            im0.segEndNodes(:, 2) == endNodes_subSeg(2)) |...
            (im0.segEndNodes(:, 1) == endNodes(1)&...
            im0.segEndNodes(:, 2) == endNodes(2))) )';
        numberOfEdges = im0.segNedges(cnctdSgs(ii));
        subSegs_maxDiams = max(im0.segDiam(cnctdSubSgs));
        %%%
        
        
        if intersect(cnctdSgs(ii), labeledSegs(:,1))
            alredyLbld_indx = [alredyLbld_indx, ii];
        end
        if im0.segDiam(cnctdSgs(ii)) <= thrshld && numberOfEdges < 3 && subSegs_maxDiams >= thrshld
            im0.segDiam(cnctdSgs(ii)) = mean([subSegs_maxDiams, im0.segDiam(sg)]);
        end
        if im0.segDiam(cnctdSgs(ii)) <= thrshld
            alredyLbld_indx = [alredyLbld_indx, ii];
        end
        
        
    end
    %%
    cnctdSgs(alredyLbld_indx) = [];
    im0.segVesType(cnctdSgs) = type;
    
    
    
    
    labeledSegs = [labeledSegs; cnctdSgs];
    l = length(labeledSegs);
    i = i+1;
end


waitbar(1, wtBr, 'Done!')
close(wtBr)



updateSelectedEdgesList(handles)
% hObject    handle to menu_labelSegments (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
end


% --------------------------------------------------------------------
function menu_showCutoff_Callback(hObject, eventdata, handles)
switch hObject.Checked
    case 'on'
        hObject.Checked = 'off';
    case 'off'
        hObject.Checked = 'on';
end
% hObject    handle to menu_showCutoff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
end


% --------------------------------------------------------------------
function menu_showMedian_Callback(hObject, eventdata, handles)
switch hObject.Checked
    case 'on'
        hObject.Checked = 'off';
    case 'off'
        hObject.Checked = 'on';
end
% hObject    handle to menu_showMedian (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
end


% --------------------------------------------------------------------
function loadSession(handles, fullPath)
global im0
hwait = waitbar(0.2, 'Loading TPM Data');
im0 = load(fullPath);
waitbar(0.9, hwait, 'Loading TPM Data');
var = who('-file', fullPath);
eval(['im0 = im0.' var{1} ';']);

waitbar(1, hwait, 'loading complete');
pause(0.2)
close(hwait)

initialiseVar()
refreshGui(handles)

highlightEdges(handles, im0.gui.edgeID)
update_Colorbar(handles)
end

% --------------------------------------------------------------------
function loadTwofiles(handles, fPath_left, fPath_right)

global im0

%-- load TPM file
hwait = waitbar(0, 'Select TPM data');
im0 = load(fPath_left, '-mat');
waitbar(0.35, hwait, 'Select TPM data');
var = who('-file', fPath_left);
eval(['im0 = im0.' var{1} ';']);
waitbar(0.4, hwait, 'loading TPM data complete');
pause(0.1)

%-- load OCT flow data
waitbar(0.5, hwait, 'Select OCT Flow Data');
OCT_Flow = load(fPath_right); %#ok<NASGU>
var = who('-file', fPath_right);
eval(['im0.OCT_Flow = OCT_Flow.' var{1} ';']);
waitbar(1, hwait, 'loading complete');
pause(0.2)
close(hwait)

if isfield(im0, 'gui')
    im0.gui = [];
end

initialiseVar()
refreshGui(handles)

highlightEdges(handles, [])
update_Colorbar(handles)
end


%------------------------------------------------------------
function openingDlg(handles)
%%
[opsHist, ~, ~] = fileparts(handles.figure1.FileName);
opsHist = fullfile(opsHist, 'recentlyAccessed.txt');

if exist(opsHist, 'file')
    tblCells = table2cell(readtable(opsHist));
else
    tblCells = {};
    lbxCells = {};
end

l = size(tblCells, 1);
if l
    ind = 1;
    while ind <= l
        if ~exist(tblCells{ind, 2}, 'file')
            tblCells(ind,:) = [];
            l = l - 1;
        else
            tblCells{ind, 3} = datenum(tblCells{ind, 1});
            ind = ind + 1;
        end
    end

    if ~isempty(tblCells)
        tblCells = sortrows(tblCells,3,'Descend');
        tblCells(:,3) = [];
    end
    lbxCells = {};

    for ind = 1:l
        [~, fnm, ~] = fileparts(tblCells{ind, 2});
        lbxCells{ind} = [datestr(tblCells{ind, 1}), '     ', fnm];
    end
    lbxCells = squeeze(lbxCells);
end


fig = uifigure('Position',[700, 300, 518, 520]);
bg = uibuttongroup(fig,'Position',[0, 0, 520, 520],...
    'SelectionChangedFcn',@bselection);
p1 = uipanel(fig, 'Title','Recently Accessed Sessions','FontSize',13,...
    'Position',[40, 350, 460, 150]);
uilabel(fig, 'Position', [5, 260, 510, 15],'Text',...
    '____________________________________________________________________________');

% Create text area for recently accessed selection
fPathBox_1 = uitextarea(fig,'Position',[40, 280, 460, 40],'Value','', 'Editable', 'off',...
    'BackgroundColor', [0.9400 0.9400 0.9400]);
lbl1 = uilabel(fig, 'Position', [50, 320, 100, 15], 'Text', 'Path to file', 'FontSize', 10);
btn1 = uibutton(fig,'push','Position',[400, 320, 100, 22],'Text', 'load',...
    'ButtonPushedFcn', @loadRecent);


% Create list box
lbox = uilistbox(p1,'Items',lbxCells,'Position',[0, 0, 460, 130],...
    'ValueChangedFcn', @updateEditField);

%- create tab group
tabgp = uitabgroup(fig, 'Position',[40, 20, 460, 220]);
tab1 = uitab(tabgp,'Title','Brows for a Session File');
tab2 = uitab(tabgp,'Title','Brows for separate files');

%- create radio buttons
rb_recent = uiradiobutton(bg,'Position',[15, 480, 20, 20], 'Text', '');
rb2_brows = uiradiobutton(bg,'Position',[15, 220, 20, 20], 'Text', '', 'Value', 0);

%- create a text area for browsed session path
fPathBox_2 = uitextarea(tab1,'Position',[15, 60, 430, 60],'Value','', 'Editable', 'off',...
    'BackgroundColor', [0.9400 0.9400 0.9400]);
lbl2 = uilabel(tab1, 'Position', [25, 120, 100, 15], 'Text', 'Path to file', 'FontSize', 10);
btn2 = uibutton(tab1,'push','Position',[395, 120, 50, 22],'Text', 'brows',...
    'ButtonPushedFcn', @bowsSession);

fPathBox_3 = uitextarea(tab2,'Position',[15, 120, 430, 40],'Value','', 'Editable', 'off',...
    'BackgroundColor', [0.9400 0.9400 0.9400]);
lbl3 = uilabel(tab2, 'Position', [25, 160, 100, 15], 'Text', 'Path to angiogram file', 'FontSize', 10);
btn3 = uibutton(tab2,'push','Position',[395, 160, 50, 22],'Text', 'brows',...
    'ButtonPushedFcn', @bowsLeftVol);

fPathBox_4 = uitextarea(tab2,'Position',[15, 50, 430, 40],'Value','', 'Editable', 'off',...
    'BackgroundColor', [0.9400 0.9400 0.9400]);
lbl4 = uilabel(tab2, 'Position', [25, 90, 100, 15], 'Text', 'Path to OCT file', 'FontSize', 10);
btn4 = uibutton(tab2,'push','Position',[395, 90, 50, 22],'Text', 'brows',...
    'ButtonPushedFcn', @bowsRightVol);

btn5 = uibutton(tab1,'push','Position',[330, 15, 100, 22],'Text', 'load',...
    'ButtonPushedFcn', @load_tab1);
btn6 = uibutton(tab2,'push','Position',[330, 15, 100, 22],'Text', 'load',...
    'ButtonPushedFcn', @load_tab2);


if isempty(lbxCells)
    rb_recent.Value = 0;
    rb2_brows.Value = 1;
    rb_recent.Enable = 'off';
    bselection()
else
    rb_recent.Value = 1;
    rb2_brows.Value = 0;
    bselection()
    updateEditField(lbox)
end

uiwait(fig)


%functions
% ValueChangedFcn callback
    function updateEditField(src, ~)
        s = src.Value;
        for indx = 1:l
            if strcmp(s, src.Items(indx))
                fPathBox_1.Value = tblCells{indx, 2};
                continue
            end
        end
    end

% radio button callbak
    function bselection(~, ~)
        if rb_recent.Value
            lbox.Enable = 'on';
            fPathBox_1.Enable = 'on';
            btn1.Enable = 'on';
            lbl1.Enable = 'on';
            
            fPathBox_2.Enable = 'off';
            btn2.Enable = 'off';
            fPathBox_3.Enable = 'off';
            btn3.Enable = 'off';
            fPathBox_4.Enable = 'off';
            btn4.Enable = 'off';
            btn5.Enable = 'off';
            btn6.Enable = 'off';
            lbl2.Enable = 'off';
            lbl3.Enable = 'off';
            lbl4.Enable = 'off';
            
        else
            lbox.Enable = 'off';
            fPathBox_1.Enable = 'off';
            btn1.Enable = 'off';
            lbl1.Enable = 'off';
            
            fPathBox_2.Enable = 'on';
            btn2.Enable = 'on';
            fPathBox_3.Enable = 'on';
            btn3.Enable = 'on';
            fPathBox_4.Enable = 'on';
            btn4.Enable = 'on';
            btn5.Enable = 'on';
            btn6.Enable = 'on';
            lbl2.Enable = 'on';
            lbl3.Enable = 'on';
            lbl4.Enable = 'on';
        end
    end



    %------
    function bowsSession(~, ~)
        if ~isempty(fPathBox_2.Value{1})
            [intialP, ~, ~] = fileparts(fPathBox_2.Value{1});
        elseif ~isempty(fPathBox_1.Value{1})
            [intialP, ~, ~] = fileparts(fPathBox_1.Value{1});
        elseif ~isempty(fPathBox_3.Value{1})
            [intialP, ~, ~] = fileparts(fPathBox_3.Value{1});
        elseif ~isempty(fPathBox_4.Value{1})
            [intialP, ~, ~] = fileparts(fPathBox_4.Value{1});
        else
            intialP = pwd;
        end
        prev_pwd = pwd;
        cd(intialP)
        [fName, fPath] = uigetfile('*.mat','Select MAT file of processing session');
        cd(prev_pwd)
        if fName == 0; return; end
        fPathBox_2.Value = fullfile(fPath, fName);
    end

    %---------
    function bowsLeftVol(~, ~)
        if ~isempty(fPathBox_3.Value{1})
            [intialP, ~, ~] = fileparts(fPathBox_3.Value{1});
        elseif ~isempty(fPathBox_4.Value{1})
            [intialP, ~, ~] = fileparts(fPathBox_4.Value{1});
        elseif ~isempty(fPathBox_1.Value{1})
            [intialP, ~, ~] = fileparts(fPathBox_1.Value{1});
        elseif ~isempty(fPathBox_2.Value{1})
            [intialP, ~, ~] = fileparts(fPathBox_2.Value{1});
        else
            intialP = pwd;
        end
        prev_pwd = pwd;
        cd(intialP)
        [fName, fPath] = uigetfile('*.img;*.tif;*.tiff;*.mat;*.seed',...
            'Load TPM Graph file');
        cd(prev_pwd)
        if fName == 0; return; end
        fPathBox_3.Value = fullfile(fPath, fName);
    end

    %---------
    function bowsRightVol(~, ~)
        if ~isempty(fPathBox_4.Value{1})
            [intialP, ~, ~] = fileparts(fPathBox_4.Value{1});
        elseif ~isempty(fPathBox_3.Value{1})
            [intialP, ~, ~] = fileparts(fPathBox_3.Value{1});
        elseif ~isempty(fPathBox_1.Value{1})
            [intialP, ~, ~] = fileparts(fPathBox_1.Value{1});
        elseif ~isempty(fPathBox_2.Value{1})
            [intialP, ~, ~] = fileparts(fPathBox_2.Value{1});
        else
            intialP = pwd;
        end
        prev_pwd = pwd;
        cd(intialP)
        [fName, fPath] = uigetfile( '*.img;*.tif;*.tiff;*.mat;*.seed',...
            'Load OCT Flow Data');
        cd(prev_pwd)
        if fName == 0; return; end
        fPathBox_4.Value = fullfile(fPath, fName);
    end


%- push button callback
    function loadRecent(~, ~)
        fullPath = fPathBox_1.Value{1};
        if ~isempty(fullPath)
            i = 1;
            while i <= l
                if strcmp(fullPath, tblCells{i, 2})
                    tblCells(i,:) = [];
                    l = l - 1;
                else
                    i = i + 1;
                end
            end
            tblCells = cell2table([tblCells; {datetime('now'), fullPath}],...
                'VariableNames', {'Last_Accessed' 'Path'});
            writetable(tblCells, opsHist)
            delete(fig)
            loadSession(handles, fullPath)
        end
    end

    function load_tab1(~, ~)
        fullPath = fPathBox_2.Value{1};
        if ~isempty(fullPath)
            i = 1;
            while i <= l
                if strcmp(fullPath, tblCells{i, 2})
                    tblCells(i,:) = [];
                    l = l - 1;
                else
                    i = i + 1;
                end
            end
            tblCells = cell2table([tblCells; {datetime('now'), fullPath}],...
                'VariableNames', {'Last_Accessed' 'Path'});
            writetable(tblCells, opsHist)
            delete(fig)
            loadSession(handles, fullPath)
        end
    end

    function load_tab2(~, ~)
        fullPath_1 = fPathBox_3.Value{1}; fullPath_2 = fPathBox_4.Value{1};
        if ~isempty(fullPath_1) && ~isempty(fullPath_2)
            delete(fig)
            loadTwofiles(handles, fullPath_1, fullPath_2)
        end
    end
end


%---------------
function monitorsBtnDownFnc(hObject, ~, handles)
%%
global im0
slctType = get(gcf, 'SelectionType');
switch slctType
    case 'normal'
        p = get(gca, 'CurrentPoint');
        p = round([p(1,1), p(1,2)]);
        cylMid = ( im0.cylNodes(im0.gui.edgeID, 1:3) + im0.cylNodes(im0.gui.edgeID, 4:6) )/2;
        switch hObject.Parent.Tag
            case 'axes_OCT_Win_01.Tag'
                step = round(p - cylMid([1, 2]));
                im0.cylNodes(im0.gui.edgeID,[1,4, 2,5]) = ...
                    im0.cylNodes(im0.gui.edgeID,[1,4, 2,5]) + [step([1,1]), step([2,2])];
                
                
            case 'axes_OCT_Win_02.Tag'
                step = round(p - cylMid([2, 3]));
                im0.cylNodes(im0.gui.edgeID,[2,5, 3,6]) = ...
                    im0.cylNodes(im0.gui.edgeID,[5,2, 3,6]) + [step([1,1]),...
                    step([2,2])];
                im0.cylNodes(im0.gui.edgeID,[2,5]) = im0.szOCT(2) - im0.cylNodes(im0.gui.edgeID,[2,5]);
                
            case 'axes_OCT_Win_03.Tag'
                step = round(p - cylMid([1, 3]));
                im0.cylNodes(im0.gui.edgeID,[1,4, 3,6]) = ...
                    im0.cylNodes(im0.gui.edgeID,[4,1, 3,6]) + [step([1,1]),...
                    step([2,2])];
                im0.cylNodes(im0.gui.edgeID,[1,4]) = im0.szOCT(1) - im0.cylNodes(im0.gui.edgeID,[1,4]); 
        end
        
        update_Monitors(handles, im0.gui.edgeID, true)
        computeEdgeFlow(handles, im0.gui.edgeID)
    case 'alt'
        %- left click
end
end


% --- Executes on key press with focus on figure1 and none of its controls.
function figure1_KeyPressFcn(~, eventdata, handles)
global im0 

switch eventdata.Key
    case 'w'
        oldVal = strsplit(get(handles.cylDiam,'String'));
        if length(oldVal) > 1
            val = str2double(cell2mat(oldVal(1)));
        else
            val = str2double(get(handles.cylDiam,'String'));
        end
        dec = val - floor(val);
        if dec < 0.5 && dec > 0
            val = val - dec + 0.5;
        elseif dec > 0.5
            val = val - dec + 1;
        elseif dec == 0 || dec == 0.5
            val = val + 0.5;
        end
        
        im0.cylDia(im0.gui.edgeID) = val;
        set(handles.cylDiam, 'String', [num2str(val), ' μm'])
        computeEdgeFlow(handles, im0.gui.edgeID)
        update_Monitors(handles, im0.gui.edgeID, true)
        
    case 's'
        oldVal = strsplit(get(handles.cylDiam,'String'));
        if length(oldVal) > 1
            val = str2double(cell2mat(oldVal(1)));
        else
            val = str2double(get(handles.cylDiam,'String'));
        end
        dec = val - floor(val);
        if dec < 0.5 && dec > 0
            val = val - dec;
        elseif dec > 0.5
            val = val - dec + 0.5;
        elseif dec == 0 || dec == 0.5
            val = val - 0.5;
        end
        
        im0.cylDia(im0.gui.edgeID) = val;
        set(handles.cylDiam, 'String', [num2str(val), ' μm'])
        computeEdgeFlow(handles, im0.gui.edgeID)
        update_Monitors(handles, im0.gui.edgeID, true)
        
    case 'q'
        l = str2double(get(handles.colorbarUpperLim,'String'));
        set(handles.colorbarUpperLim,'String',num2str(l + 0.5));
        colorbarUpperLim_Callback(handles.colorbarUpperLim, [], handles)
    case 'a'
        l = str2double(get(handles.colorbarUpperLim,'String'));
        if l + 0.5 > 0
        set(handles.colorbarUpperLim,'String',num2str(l - 0.5));
        colorbarUpperLim_Callback(handles.colorbarUpperLim, [], handles)
        end
        
    case 'delete'
        rejectEdge_Callback([], [], handles)
        
    case 'return'
        %         moveListSelection(handles, false, false)
    case '3'
        set(handles.sgmntType, 'String', '3')
        sgmntType_Callback(handles.sgmntType, [], handles)
    case '1'
        set(handles.sgmntType, 'String', '1')
        sgmntType_Callback(handles.sgmntType, [], handles)
    case '2'
        set(handles.sgmntType, 'String', '2')
        sgmntType_Callback(handles.sgmntType, [], handles)
end


% hObject    handle to figure1 (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.FIGURE)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)
end


% --- Executes on key press with focus on win03Left and none of its controls.
function win03Left_KeyPressFcn(~, eventdata, handles)
figure1_KeyPressFcn([], eventdata, handles)
end

% --- Executes on key press with focus on win03Right and none of its controls.
function win03Right_KeyPressFcn(~, eventdata, handles)
figure1_KeyPressFcn([], eventdata, handles)
end

% --- Executes on key press with focus on win03Down and none of its controls.
function win03Down_KeyPressFcn(~, eventdata, handles)
figure1_KeyPressFcn([], eventdata, handles)
end

% --- Executes on key press with focus on win03Up and none of its controls.
function win03Up_KeyPressFcn(~, eventdata, handles)
figure1_KeyPressFcn([], eventdata, handles)
end

% --- Executes on key press with focus on win03RotAntiClock and none of its controls.
function win03RotAntiClock_KeyPressFcn(~, eventdata, handles)
figure1_KeyPressFcn([], eventdata, handles)
end

% --- Executes on key press with focus on win03RotClock and none of its controls.
function win03RotClock_KeyPressFcn(~, eventdata, handles)
figure1_KeyPressFcn([], eventdata, handles)
end

% --- Executes on key press with focus on win01Left and none of its controls.
function win01Left_KeyPressFcn(~, eventdata, handles)
figure1_KeyPressFcn([], eventdata, handles)
end

% --- Executes on key press with focus on win01Right and none of its controls.
function win01Right_KeyPressFcn(~, eventdata, handles)
figure1_KeyPressFcn([], eventdata, handles)
end

% --- Executes on key press with focus on win01Up and none of its controls.
function win01Up_KeyPressFcn(~, eventdata, handles)
figure1_KeyPressFcn([], eventdata, handles)
end

% --- Executes on key press with focus on win01Down and none of its controls.
function win01Down_KeyPressFcn(~, eventdata, handles)
figure1_KeyPressFcn([], eventdata, handles)
end

% --- Executes on key press with focus on win01RotClock and none of its controls.
function win01RotClock_KeyPressFcn(~, eventdata, handles)
figure1_KeyPressFcn([], eventdata, handles)
end

% --- Executes on key press with focus on win01RotAntiClock and none of its controls.
function win01RotAntiClock_KeyPressFcn(~, eventdata, handles)
figure1_KeyPressFcn([], eventdata, handles)
end

% --- Executes on key press with focus on win02Left and none of its controls.
function win02Left_KeyPressFcn(~, eventdata, handles)
figure1_KeyPressFcn([], eventdata, handles)
end

% --- Executes on key press with focus on win02Right and none of its controls.
function win02Right_KeyPressFcn(~, eventdata, handles)
figure1_KeyPressFcn([], eventdata, handles)
end

% --- Executes on key press with focus on win02Up and none of its controls.
function win02Up_KeyPressFcn(~, eventdata, handles)
figure1_KeyPressFcn([], eventdata, handles)
end

% --- Executes on key press with focus on win02Down and none of its controls.
function win02Down_KeyPressFcn(~, eventdata, handles)
figure1_KeyPressFcn([], eventdata, handles)
end

% --- Executes on key press with focus on win02RotClock and none of its controls.
function win02RotClock_KeyPressFcn(~, eventdata, handles)
figure1_KeyPressFcn([], eventdata, handles)
end

% --- Executes on key press with focus on win02RotAntiClock and none of its controls.
function win02RotAntiClock_KeyPressFcn(~, eventdata, handles)
figure1_KeyPressFcn([], eventdata, handles)
end

% --- Executes on key press with focus on cylDiam and none of its controls.
function cylDiam_KeyPressFcn(~, eventdata, handles)
figure1_KeyPressFcn([], eventdata, handles)
end

% --- Executes on key press with focus on selectedEdgesList and none of its controls.
function selectedEdgesList_KeyPressFcn(~, eventdata, handles)
figure1_KeyPressFcn([], eventdata, handles)
end

% --- Executes on key press with focus on colorbarUpperLim and none of its controls.
function colorbarUpperLim_KeyPressFcn(~, eventdata, handles)
figure1_KeyPressFcn([], eventdata, handles)
end

% --- Executes on key press with focus on colorbalLowerLim and none of its controls.
function colorbalLowerLim_KeyPressFcn(~, eventdata, handles)
figure1_KeyPressFcn([], eventdata, handles)
end

% --- Executes on key press with focus on rejectEdge and none of its controls.
function rejectEdge_KeyPressFcn(~, eventdata, handles)
figure1_KeyPressFcn([], eventdata, handles)
end

% --- Executes on key press with focus on rstCylPos and none of its controls.
function rstCylPos_KeyPressFcn(~, eventdata, handles)
figure1_KeyPressFcn([], eventdata, handles)
end

% --- Executes on key press with focus on radiobutton_XY and none of its controls.
function radiobutton_XY_KeyPressFcn(~, eventdata, handles)
figure1_KeyPressFcn([], eventdata, handles)
end

% --- Executes on key press with focus on radiobutton_YZ and none of its controls.
function radiobutton_YZ_KeyPressFcn(~, eventdata, handles)
figure1_KeyPressFcn([], eventdata, handles)
end

% --- Executes on key press with focus on radiobutton_XZ and none of its controls.
function radiobutton_XZ_KeyPressFcn(~, eventdata, handles)
figure1_KeyPressFcn([], eventdata, handles)
end
