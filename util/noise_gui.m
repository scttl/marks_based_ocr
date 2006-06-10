function varargout = noise_gui(varargin)
%  NOISE_GUI  Lets one see and select image noise model parameters
%
%    [rem_pct, add_pct, thresh, fspec_args] = noise_gui img_file 
%
%    This GUI takes in a matrix version of an image file (or its filename),
%    and displays it in the top axis.  After selecting appropriate parameter
%    settings, a noisy version of the image is displayed in the bottom axis.
%    The user can then manipulate sliders and values to come up with an 
%    appropriate noise model, then close the GUI, which returns the parameter
%    values to the caller.
% 
%    img_file should either be the full path and name of a valid image file,
%    or it should be a binary matrix representing pixel values.  (Colour images
%    are converted to binary).
%
%    rem_pct, add_pct give the percentage of on (respectively off) pixels that
%    had there value randomly switched.
%
%    thresh determines the minimum pixel values (which are graryscale and lie 
%    between 0 and 1), to consider as being 'on'.
%
%    fspec_args is a cell array which can be passed as a sequence of arguments
%    to add_img_noise or fspecial to create an appropriate sized filter mask
%    that can be used to blur/sharpen an image.
%
%    See: add_img_noise and fspecial for more details.
%
%    If editing type "help guide" to see how this code was originally created
%

% CVS INFO %
%%%%%%%%%%%%
% $Id: noise_gui.m,v 1.1 2006-06-10 21:01:47 scottl Exp $
%
% REVISION HISTORY
% $Log: noise_gui.m,v $
% Revision 1.1  2006-06-10 21:01:47  scottl
% Initial revision.
%


% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @noise_gui_OpeningFcn, ...
                   'gui_OutputFcn',  @noise_gui_OutputFcn, ...
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


% --- Executes just before noise_gui is made visible.
function noise_gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to noise_gui (see VARARGIN)

if length(varargin) ~= 1
    error('must specify the image file, or its representation as an argument');
end
%initialize the output object
handles.output = {0.0, 0.0, .5, {'gaussian', [3 3], .5}};

%initialize the input object
if ischar(varargin{1})
    %open and read the image file
    handles.orig_img = imread(varargin{1});
else
    handles.orig_img = varargin{1};
end

% Update handles structure
guidata(hObject, handles);

movegui(hObject,'onscreen'); % To display application onscreen
movegui(hObject,'center');  % To display application in the center of screen

set(handles.orig_axs,'HandleVisibility','ON');
set(handles.noise_axs,'HandleVisibility','OFF');
imagesc(handles.orig_img);
colormap gray;
axis equal;
axis tight;
axis off;
set(handles.orig_axs,'XTickLabel',' ','YTickLabel',' ');
set(handles.noise_axs,'XTickLabel',' ','YTickLabel',' ');

%draw the noisy version
draw_noisy(handles);

%don't return right away
set(handles.figure1, 'Visible', 'ON');
waitfor(handles.figure1, 'Visible', 'OFF');


% --- Outputs from this function are returned to the command line.
function varargout = noise_gui_OutputFcn(hObject, eventdata, handles) 
% Pass the current output values to the command line
varargout = handles.output;
close(gcf); % to close GUI


% --- Executes on button press in apply_btn.
function apply_btn_Callback(hObject, eventdata, handles)
set(handles.figure1, 'Visible', 'OFF');


% --- Executes on button press in cncl_btn.
function cncl_btn_Callback(hObject, eventdata, handles)
handles.output{1} = -1;  %return -1 to the caller
guidata(hObject, handles);
set(handles.figure1, 'Visible', 'OFF');


% --- Executes on slider movement.
function rem_sld_Callback(hObject, eventdata, handles)
% update the display in the text box
val = get(hObject, 'Value');
set(handles.rem_val, 'String', num2str(val));
handles.output{1} = val;

guidata(hObject, handles);
draw_noisy(handles);


% --- Executes during object creation, after setting all properties.
function rem_sld_CreateFcn(hObject, eventdata, handles)
usewhitebg = 1;
if usewhitebg
    set(hObject,'BackgroundColor',[.9 .9 .9]);
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on slider movement.
function add_sld_Callback(hObject, eventdata, handles)
% update the display in the text box
val = get(hObject, 'Value');
set(handles.add_val, 'String', num2str(val));
handles.output{2} = val;

guidata(hObject, handles);
draw_noisy(handles);


% --- Executes during object creation, after setting all properties.
function add_sld_CreateFcn(hObject, eventdata, handles)
usewhitebg = 1;
if usewhitebg
    set(hObject,'BackgroundColor',[.9 .9 .9]);
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on slider movement.
function thresh_sld_Callback(hObject, eventdata, handles)
% update the display in the text box
val = get(hObject, 'Value');
set(handles.thresh_val, 'String', num2str(val));
handles.output{3} = val;

guidata(hObject, handles);
draw_noisy(handles);


% --- Executes during object creation, after setting all properties.
function thresh_sld_CreateFcn(hObject, eventdata, handles)
usewhitebg = 1;
if usewhitebg
    set(hObject,'BackgroundColor',[.9 .9 .9]);
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


function rem_val_Callback(hObject, eventdata, handles)
%update the slider position
valid_txt_pct(hObject);
set(handles.rem_sld, 'Value', str2double(get(hObject, 'String')));

guidata(hObject, handles);
draw_noisy(handles);


% --- Executes during object creation, after setting all properties.
function rem_val_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


function add_val_Callback(hObject, eventdata, handles)
%update the slider position
valid_txt_pct(hObject);
set(handles.add_sld, 'Value', str2double(get(hObject, 'String')));

guidata(hObject, handles);
draw_noisy(handles);


% --- Executes during object creation, after setting all properties.
function add_val_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


function thresh_val_Callback(hObject, eventdata, handles)
%update the slider position
valid_txt_pct(hObject);
set(handles.thresh_sld, 'Value', str2double(get(hObject, 'String')));

guidata(hObject, handles);
draw_noisy(handles);


% --- Executes during object creation, after setting all properties.
function thresh_val_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on selection change in model_selection.
function model_selection_Callback(hObject, eventdata, handles)
%update the display to reflect the new model, and update the output params
switch get(handles.model_selection, 'Value')
case 1
    type = 'gaussian';
    %requires mask size, par1 = sigma
    set(handles.mask_size_popup, 'Visible', 'ON');
    set(handles.mask_size_txt, 'Visible', 'ON');
    set(handles.par1_sld, 'Visible', 'ON');
    set(handles.par1_sld, 'Min', 0.1);
    set(handles.par1_sld, 'Max', 5.0);
    set(handles.par1_val, 'Visible', 'ON');
    set(handles.par1_txt, 'Visible', 'ON');
    set(handles.par1_txt, 'String', 'Sigma');
    sz = get(handles.mask_size_popup, 'Value');
    sigma = get(handles.par1_sld, 'Value');
    if sigma < 0.1
        sigma = 0.1;
        set(handles.par1_sld, 'Value', sigma);
        set(handles.par1_val, 'String', num2str(sigma));
    elseif sigma > 5.0
        sigma = 5.0;
        set(handles.par1_sld, 'Value', sigma);
        set(handles.par1_val, 'String', num2str(sigma));
    end
    params = {[sz sz], sigma};
case 2
    type = 'laplacian';
    %no mask size, par1 = alpha
    set(handles.mask_size_popup, 'Visible', 'OFF');
    set(handles.mask_size_txt, 'Visible', 'OFF');
    set(handles.par1_sld, 'Visible', 'ON');
    valid_txt_pct(handles.par1_val);
    valid_txt_pct(handles.par1_sld);
    set(handles.par1_sld, 'Min', 0.0);
    set(handles.par1_sld, 'Max', 1.0);
    set(handles.par1_val, 'Visible', 'ON');
    set(handles.par1_txt, 'Visible', 'ON');
    set(handles.par1_txt, 'String', 'Alpha');
    alpha = get(handles.par1_sld, 'Value');
    if alpha < 0.0
        alpha = 0.0;
        set(handles.par1_sld, 'Value', alpha);
        set(handles.par1_val, 'String', num2str(alpha));
    elseif alpha > 1.0
        alpha = 1.0;
        set(handles.par1_sld, 'Value', alpha);
        valid_txt_pct(handles.par1_val);
        set(handles.par1_val, 'String', num2str(alpha));
    end
    params = {alpha};
case 3
    type = 'unsharp';
    %no mask size, par1 = alpha
    set(handles.mask_size_popup, 'Visible', 'OFF');
    set(handles.mask_size_txt, 'Visible', 'OFF');
    set(handles.par1_sld, 'Visible', 'ON');
    valid_txt_pct(handles.par1_val);
    set(handles.par1_sld, 'Min', 0.0);
    set(handles.par1_sld, 'Max', 1.0);
    set(handles.par1_val, 'Visible', 'ON');
    set(handles.par1_txt, 'Visible', 'ON');
    set(handles.par1_txt, 'String', 'Alpha');
    alpha = get(handles.par1_sld, 'Value');
    if alpha < 0.0
        alpha = 0.0;
        set(handles.par1_sld, 'Value', alpha);
        set(handles.par1_val, 'String', num2str(alpha));
    elseif alpha > 1.0
        alpha = 1.0;
        set(handles.par1_sld, 'Value', alpha);
        valid_txt_pct(handles.par1_val);
        set(handles.par1_val, 'String', num2str(alpha));
    end
    params = {alpha};
case 4
    type = 'sobel';
    % no mask size, no par1
    set(handles.mask_size_popup, 'Visible', 'OFF');
    set(handles.mask_size_txt, 'Visible', 'OFF');
    set(handles.par1_sld, 'Visible', 'OFF');
    set(handles.par1_val, 'Visible', 'OFF');
    set(handles.par1_txt, 'Visible', 'OFF');
    params = {};
case 5
    type = 'prewitt';
    % no mask size, no par1
    set(handles.mask_size_popup, 'Visible', 'OFF');
    set(handles.mask_size_txt, 'Visible', 'OFF');
    set(handles.par1_sld, 'Visible', 'OFF');
    set(handles.par1_val, 'Visible', 'OFF');
    set(handles.par1_txt, 'Visible', 'OFF');
    params = {};
case 6
    type = 'average';
    % no par1
    set(handles.mask_size_popup, 'Visible', 'ON');
    set(handles.mask_size_txt, 'Visible', 'ON');
    set(handles.par1_sld, 'Visible', 'OFF');
    set(handles.par1_val, 'Visible', 'OFF');
    set(handles.par1_txt, 'Visible', 'OFF');
    sz = get(handles.mask_size_popup, 'Value');
    params = {[sz sz]};
end
handles.output{4} = {type, params{:}};

guidata(hObject, handles);
draw_noisy(handles);

% --- Executes during object creation, after setting all properties.
function model_selection_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on slider movement.
function par1_sld_Callback(hObject, eventdata, handles)
%update text value above slider
val = get(hObject, 'Value');
set(handles.par1_val, 'String', num2str(val));
switch get(handles.model_selection, 'Value')
case 1
    %gaussian -> has mask size parameter before
    handles.output{4}{3} = val;
case {2, 3}
    %laplacian or unsharp -> has no mask parameter
    handles.output{4}{2} = val;
otherwise
    error('not sure how you got here, but something is messed!');
end

guidata(hObject, handles);
draw_noisy(handles);


% --- Executes during object creation, after setting all properties.
function par1_sld_CreateFcn(hObject, eventdata, handles)
usewhitebg = 1;
if usewhitebg
    set(hObject,'BackgroundColor',[.9 .9 .9]);
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function par1_val_Callback(hObject, eventdata, handles)
%update the slider position
switch get(handles.model_selection, 'Value')
case 1
    %gaussian -> sigma must be positive, and less than 5 (or maximum value)
    valid_sigma(hObject);
case {2, 3}
    %laplacian or unsharp -> alpha must be between 0 and 1
    valid_txt_pct(hObject);
end
set(handles.par1_sld, 'Value', str2double(get(hObject, 'String')));

guidata(hObject, handles);
draw_noisy(handles);


% --- Executes during object creation, after setting all properties.
function par1_val_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on selection change in mask_size_popup.
function mask_size_popup_Callback(hObject, eventdata, handles)
%update the size of the mask matrix 
sz = get(hObject, 'Value');
handles.output{4}{2} = [sz sz];

guidata(hObject, handles);
draw_noisy(handles);



% --- Executes during object creation, after setting all properties.
function mask_size_popup_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function valid_txt_pct(hObject)
% ensure that the text value is set to a valid percentage value i.e. between
% 0 and 1
val = str2double(get(hObject, 'String'));
if val < 0
    set(hObject, 'String', '0.0');
elseif val > 1
    set(hObject, 'String', '1.0');
end


function valid_sigma(hObject)
% ensure that the text value is set to a valid positive value i.e. larger
% than 0 and less than or equal to 5
val = str2double(get(hObject, 'String'));
if val <= 0
    set(hObject, 'String', '0.1');
elseif val > 5
    set(hObject, 'String', '5.0');
end


function draw_noisy(handles)
% re-draw the noisy version of the image on the lower axis
set(handles.noise_axs,'HandleVisibility','ON');
noise_img = add_img_noise(handles.orig_img, handles.output{1:3}, ...
            handles.output{4}{:});
set(handles.orig_axs,'HandleVisibility','OFF');
set(handles.noise_axs,'HandleVisibility','ON');
imagesc(noise_img);
colormap gray;
axis equal;
axis tight;
axis off;
set(handles.orig_axs,'XTickLabel',' ','YTickLabel',' ');
set(handles.noise_axs,'XTickLabel',' ','YTickLabel',' ');
