function varargout = tuneSegmentationParams(varargin)
% TUNESEGMENTATIONPARAMS MATLAB code for tuneSegmentationParams.fig
%      TUNESEGMENTATIONPARAMS, by itself, creates a new TUNESEGMENTATIONPARAMS or raises the existing
%      singleton*.
%
%      H = TUNESEGMENTATIONPARAMS returns the handle to a new TUNESEGMENTATIONPARAMS or the handle to
%      the existing singleton*.
%
%      TUNESEGMENTATIONPARAMS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in TUNESEGMENTATIONPARAMS.M with the given input arguments.
%
%      TUNESEGMENTATIONPARAMS('Property','Value',...) creates a new TUNESEGMENTATIONPARAMS or raises
%      the existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before tuneSegmentationParams_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to tuneSegmentationParams_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help tuneSegmentationParams

% Last Modified by GUIDE v2.5 07-Jul-2016 12:32:55

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @tuneSegmentationParams_OpeningFcn, ...
                   'gui_OutputFcn',  @tuneSegmentationParams_OutputFcn, ...
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

% --- Executes just before tuneSegmentationParams is made visible.
function tuneSegmentationParams_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to tuneSegmentationParams (see VARARGIN)

% set default variables upon opening
defaults = struct();
defaults.startframe=1;
defaults.endframe=-1;
defaults.threshscl=0.5;
defaults.adjmin=0.01;
defaults.adjmax=0.99;

% input default variables
for index = 1:2:nargin-3
    switch(lower(varargin{index}))
        case 'storagedirpath'
             set(handles.storagedirbox,'String',varargin{index+1});            
        case 'storagedirname'
            set(handles.storageglobbox,'String',varargin{index+1});             
        case 'filename'
            set(handles.tifglobbox,'String',varargin{index+1});
        case 'metadata'
            set(handles.metadatabox,'String',varargin{index+1});
        case 'outputfile'
            set(handles.outputfilebox,'String',varargin{index+1});            
        case 'mnum'            
            set(handles.mnumbox,'String',varargin{index+1});            
        case 'startframe'
            defaults.startframe=varargin{index+1};
        case 'endframe'
            defaults.endframe=varargin{index+1};
        case 'threshscl'
            defaults.threshscl=varargin{index+1};
        case 'adjmin'
            defaults.adjmin=varargin{index+1};
        case 'adjmax'
            defaults.adjmax=varargin{index+1};
    end
end

handles.segoptions = struct();
handles.defaults = defaults;

% Choose default command line output for tuneSegmentationParams
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

setDefaults(hObject,handles);

%initialize_gui(hObject, handles, false);

function setDefaults(hObject,handles)
    set(handles.startframebox,'String',handles.defaults.startframe);
    
    % set endframe to last frame if movie is loaded and negative indexing
    if (isfield(handles,'cellmv') & handles.defaults.endframe<0)
        endframe = length(handles.cellmv)+handles.defaults.endframe+1;        
    else        
        endframe=handles.defaults.endframe;
    end
    set(handles.endframebox,'String',endframe);
    set(handles.threshsclbox,'String',handles.defaults.threshscl);
    set(handles.adjminbox,'String',handles.defaults.adjmin);
    set(handles.adjmaxbox,'String',handles.defaults.adjmax);
    
    guidata(hObject, handles);
    
% UIWAIT makes tuneSegmentationParams wait for user response (see UIRESUME)
% uiwait(handles.fullgui);

function loadCellmv(hObject,handles)
    % load in a information about the entire movie
    
    % which movie to load
    mnum =  str2num(get(handles.mnumbox,'String'));
    
    % get storage directory
    storagedirpath = get(handles.storagedirbox,'String');
    storagedirglob = get(handles.storageglobbox,'String');
    storagedirglob = fullfile(storagedirpath, sprintf(storagedirglob,mnum));
    storagedirname = dir(storagedirglob); 
    if isempty(storagedirname)
        msgbox(sprintf('Failed to find dir %s', storagedirglob),'Error')
        cla(handles.axes_first)
        cla(handles.axes_last)
        return
    end
    storagedirname = fullfile(storagedirname(1).folder, filesep)
    disp('Storagedirname')
    disp(storagedirname)
    
    tifglob = get(handles.tifglobbox,'String');
    metafileglob = get(handles.metadatabox,'String');
    
    % get metadata file
    metafile = dir(fullfile(storagedirname, metafileglob));
    disp(fullfile(storagedirname, metafileglob))
    metafile = metafile.name
    disp(sprintf('Processing data for Movie %d',mnum))
    
    % load in tif stacks
    files = dir(fullfile(storagedirname, '*.ome.tif'));
    files = {files.name}
    
    % read in movie data
    disp('Reading in data from file(s)...')
    cellmv = loadMovieData(storagedirname,files,metafile);
    
    set(handles.loadedfilebox,'String',files{end})
    handles.cellmv=cellmv;
    
    % update endframe if negative
    endframe = str2num(get(handles.endframebox,'String'));    
    if endframe<0
        endframe = length(cellmv) + endframe + 1;
        set(handles.endframebox,'String',endframe)        
    end
    
    fc1 = str2num(get(handles.startframebox,'String'));
    fc2 = endframe
    loadImg(hObject,handles,fc1,fc2)
    
    guidata(hObject, handles);
    
function loadImg(hObject,handles,fc1,fc2)
% Load in images into the axes

img1 = imread(handles.cellmv(fc1).fname,handles.cellmv(fc1).fpage);
axes(handles.axes_first);
imshow(img1,[],'initialmagnification','fit');
title(sprintf('Frame %d', fc1))

if (fc2>length(handles.cellmv))
    msgbox(sprintf('Could not find frame %d',fc2),'Error')
    cla(handles.axes_last)
else
    img2 = imread(handles.cellmv(fc2).fname,handles.cellmv(fc2).fpage);
    axes(handles.axes_last);
    imshow(img2,[],'initialmagnification','fit');
    title(sprintf('Frame %d', fc2))
end

% --- Outputs from this function are returned to the command line.
function varargout = tuneSegmentationParams_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in segmentationbutton.
function segmentationbutton_Callback(hObject, eventdata, handles)
% hObject    handle to segmentationbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if (~isfield(handles,'cellmv'))
    msgbox('No cell movie loaded!','Error')
    cla(handles.axes_first)
    cla(handles.axes_last)
    return
end

% get the segmentation parameters
handles.segoptions.threshscl = str2num(get(handles.threshsclbox,'String'))
handles.segoptions.dodisplay = 1;
handles.segoptions.adjmin = str2num(get(handles.adjminbox,'String'));
handles.segoptions.adjmax = str2num(get(handles.adjmaxbox,'String'));

startframe = str2num(get(handles.startframebox,'String'));
endframe = str2num(get(handles.endframebox,'String'));

% Run segmentation on first frame
axes(handles.axes_first)
[handles.cellmv, handles.segoptions] = segmentLysotracker(handles.cellmv,[startframe],'',handles.segoptions);
axes(handles.axes_last)
[handles.cellmv, handles.segoptions] = segmentLysotracker(handles.cellmv,[endframe],'',handles.segoptions);

guidata(hObject,handles)

% --- Executes on button press in reset.
function reset_Callback(hObject, eventdata, handles)
% hObject    handle to reset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

initialize_gui(gcbf, handles, true);

% --- Executes when selected object changed in unitgroup.
function unitgroup_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in unitgroup 
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if (hObject == handles.english)
    set(handles.text4, 'String', 'lb/cu.in');
    set(handles.text5, 'String', 'cu.in');
    set(handles.text6, 'String', 'lb');
else
    set(handles.text4, 'String', 'kg/cu.m');
    set(handles.text5, 'String', 'cu.m');
    set(handles.text6, 'String', 'kg');
end

% --------------------------------------------------------------------
% function initialize_gui(fig_handle, handles, isreset)
% % If the metricdata field is present and the reset flag is false, it means
% % we are we are just re-initializing a GUI by calling it from the cmd line
% % while it is up. So, bail out as we dont want to reset the data.
% if isfield(handles, 'metricdata') && ~isreset
%     return;
% end
% 
% handles.metricdata.density = 0;
% handles.metricdata.volume  = 0;
% 
% set(handles.density, 'String', handles.metricdata.density);
% set(handles.volume,  'String', handles.metricdata.volume);
% set(handles.mass, 'String', 0);
% 
% set(handles.unitgroup, 'SelectedObject', handles.english);
% 
% set(handles.text4, 'String', 'lb/cu.in');
% set(handles.text5, 'String', 'cu.in');
% set(handles.text6, 'String', 'lb');
% 
% % Update handles structure
% guidata(handles.fullgui, handles);


% --- Executes on button press in savebutton.
function savebutton_Callback(hObject, eventdata, handles)
% hObject    handle to savebutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

mnum = str2num(get(handles.mnumbox,'String'));
outfilename = get(handles.outputfilebox,'String');
outfilename = sprintf(outfilename,mnum);
disp(sprintf('Saving to file %s', outfilename))

cellmv = handles.cellmv;
segoptions = handles.segoptions;
save(outfilename,'cellmv','segoptions')

% --- Executes on button press in resetbutton.
function resetbutton_Callback(hObject, eventdata, handles)
% hObject    handle to resetbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

setDefaults(hObject,handles)

function storagedirbox_Callback(hObject, eventdata, handles)
% hObject    handle to storagedirbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of storagedirbox as text
%        str2double(get(hObject,'String')) returns contents of storagedirbox as a double


% --- Executes during object creation, after setting all properties.
function storagedirbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to storagedirbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function outputfilebox_Callback(hObject, eventdata, handles)
% hObject    handle to outputfilebox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of outputfilebox as text
%        str2double(get(hObject,'String')) returns contents of outputfilebox as a double


% --- Executes during object creation, after setting all properties.
function outputfilebox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to outputfilebox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function mnumbox_Callback(hObject, eventdata, handles)
% hObject    handle to mnumbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of mnumbox as text
%        str2double(get(hObject,'String')) returns contents of mnumbox as a double


% --- Executes during object creation, after setting all properties.
function mnumbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mnumbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in prevbutton.
function prevbutton_Callback(hObject, eventdata, handles)
% hObject    handle to prevbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

mnum = str2num(get(handles.mnumbox,'String'));
mnum=mnum-1;
set(handles.mnumbox,'String',mnum)
% update last frame if necessary
endframe = str2num(get(handles.endframebox,'String'));
if endframe>=length(handles.cellmv)
    set(handles.endframebox,'String',-1);
end
loadCellmv(hObject,handles)

% --- Executes on button press in nextbutton.
function nextbutton_Callback(hObject, eventdata, handles)
% hObject    handle to nextbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

mnum = str2num(get(handles.mnumbox,'String'));
mnum=mnum+1;

set(handles.mnumbox,'String',mnum)

% update last frame if necessary
endframe = str2num(get(handles.endframebox,'String'));
if endframe>=length(handles.cellmv)
    set(handles.endframebox,'String',-1);
end
loadCellmv(hObject,handles)

function startframebox_Callback(hObject, eventdata, handles)
% hObject    handle to startframebox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of startframebox as text
%        str2double(get(hObject,'String')) returns contents of startframebox as a double


% --- Executes during object creation, after setting all properties.
function startframebox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to startframebox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function threshsclbox_Callback(hObject, eventdata, handles)
% hObject    handle to threshsclbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of threshsclbox as text
%        str2double(get(hObject,'String')) returns contents of threshsclbox as a double


% --- Executes during object creation, after setting all properties.
function threshsclbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to threshsclbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function adjminbox_Callback(hObject, eventdata, handles)
% hObject    handle to adjminbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of adjminbox as text
%        str2double(get(hObject,'String')) returns contents of adjminbox as a double


% --- Executes during object creation, after setting all properties.
function adjminbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to adjminbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function adjmaxbox_Callback(hObject, eventdata, handles)
% hObject    handle to adjmaxbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of adjmaxbox as text
%        str2double(get(hObject,'String')) returns contents of adjmaxbox as a double


% --- Executes during object creation, after setting all properties.
function adjmaxbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to adjmaxbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in loadbutton.
function loadbutton_Callback(hObject, eventdata, handles)
% hObject    handle to loadbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

loadCellmv(hObject,handles);

function endframebox_Callback(hObject, eventdata, handles)
% hObject    handle to endframebox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of endframebox as text
%        str2double(get(hObject,'String')) returns contents of endframebox as a double


% --- Executes during object creation, after setting all properties.
function endframebox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to endframebox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function tifglobbox_Callback(hObject, eventdata, handles)
% hObject    handle to tifglobbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tifglobbox as text
%        str2double(get(hObject,'String')) returns contents of tifglobbox as a double


% --- Executes during object creation, after setting all properties.
function tifglobbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tifglobbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function metadatabox_Callback(hObject, eventdata, handles)
% hObject    handle to metadatabox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of metadatabox as text
%        str2double(get(hObject,'String')) returns contents of metadatabox as a double


% --- Executes during object creation, after setting all properties.
function metadatabox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to metadatabox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function storageglobbox_Callback(hObject, eventdata, handles)
% hObject    handle to storageglobbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of storageglobbox as text
%        str2double(get(hObject,'String')) returns contents of storageglobbox as a double


% --- Executes during object creation, after setting all properties.
function storageglobbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to storageglobbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in loadframesbutton.
function loadframesbutton_Callback(hObject, eventdata, handles)
% hObject    handle to loadframesbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if (~isfield(handles,'cellmv'))
    msgbox('No movie is loaded','Error')
    return
end

startframe= str2num(get(handles.startframebox,'String'));
endframe = str2num(get(handles.endframebox,'String'));
loadImg(hObject,handles,startframe,endframe)
