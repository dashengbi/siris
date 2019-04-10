function varargout = Choose_Spindles(varargin)
% CHOOSE_SPINDLES MATLAB code for Choose_Spindles.fig
%      CHOOSE_SPINDLES, by itself, creates a new CHOOSE_SPINDLES or raises the existing
%      singleton*.
%
%      H = CHOOSE_SPINDLES returns the handle to a new CHOOSE_SPINDLES or the handle to
%      the existing singleton*.
%
%      CHOOSE_SPINDLES('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CHOOSE_SPINDLES.M with the given input arguments.
%
%      CHOOSE_SPINDLES('Property','Value',...) creates a new CHOOSE_SPINDLES or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Choose_Spindles_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Choose_Spindles_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Choose_Spindles

% Last Modified by GUIDE v2.5 17-Dec-2018 07:13:55

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Choose_Spindles_OpeningFcn, ...
                   'gui_OutputFcn',  @Choose_Spindles_OutputFcn, ...
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

% --- Executes just before Choose_Spindles is made visible.
function Choose_Spindles_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Choose_Spindles (see VARARGIN)

% Choose default command line output for Choose_Spindles
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);
clear global;
global loaded_data loaded_scores loaded_spindles plotted num_comp color_plot name_comp;
loaded_data = false;
loaded_scores = false;
loaded_spindles = false;
plotted = false;
num_comp = 0;
color_plot = [];
name_comp = [];

% UIWAIT makes Choose_Spindles wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = Choose_Spindles_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

%% ==========Start helper functions==========
function y = pre_process(data)

global f fs nfft win;

noise_band = find(f>90,1):find(f<100,1,'last');
temp_pow = pwelch(data,hanning(win),[],nfft,fs);
noise_lev = mean(temp_pow(noise_band));

data_smooth = smooth(data,10);
data_base = smooth(data_smooth,2000);
data_base = smooth(data_base,2000);
data_s = data_smooth - data_base;

y = data_s./sqrt(noise_lev);

function data_all = find_seg_2(ind)
data_seg = [];
if length(ind) > 1
    diff_ind = diff(ind); diff_ind2 = [diff_ind(2:end);0];
    if diff_ind(1) == 1
        data_seg(1,1) = ind(1);
    end

    count = 1;
    for i=1:length(diff_ind)
        if (diff_ind(i) ~= 1) && (diff_ind2(i) == 1)
            data_seg(count,1) = ind(i+1);
        elseif (diff_ind(i) == 1) && (diff_ind2(i) ~= 1) 
            data_seg(count,2) = ind(i+1);
            count = count + 1;
        end
    end

    if diff_ind(end) == 1
        data_seg(end,2) = ind(end);
    end
    
    data_segs = [];
    for iseg = 1 : size(data_seg, 1)
        data_segs = [data_segs, data_seg(iseg,1): data_seg(iseg,2)];
    end
    data_individuals = setdiff(ind, data_segs);
    data_individual = [];
    count_individual = 1;
    for j = 1: numel(data_individuals)
        data_individual(count_individual, 1) = data_individuals(j);
        data_individual(count_individual, 2) = data_individuals(j);
        count_individual = count_individual + 1;
    end
    data_all = sort([data_seg; data_individual], 1); 
end

function cleared = clear_window()
global temp_points temp_plots_n temp_plots_f;
delete(temp_plots_n(:));
temp_plots_n = [];
delete(temp_plots_f(:));
temp_plots_f = [];
temp_points = [];
cleared = 1;

function [spindle_in, count] = count_win_spindles(num_comp, spindle_comp, window_ind, window_len, fs)
spindle_in = cell(num_comp, 1);
count = 0;
for i=1:num_comp
    spindle_in{i} = [];
    for j=1:size(spindle_comp{i}, 1)
        start = spindle_comp{i}(j, 1);
        stop = spindle_comp{i}(j, 2);
        if ((start >= window_ind) && (start <= window_ind + window_len * fs)) || ((stop >= window_ind) && (stop <= window_ind + window_len * fs))
            spindle_in{i} = cat(1, spindle_in{i}, [start stop]);
        end
    end
    count = count + size(spindle_in{i}, 1);
end

function changed = change_view(handles, window_len, window_ind, fs)
axis_labels = cell(1, window_len + 1);
for i = 0:window_len
    axis_labels{i + 1} = (window_ind / fs) + i;
end
handles.normalized_signal.XLim = [window_ind, window_ind + window_len * fs];
handles.normalized_signal.XTick = window_ind:fs:(window_ind + window_len * fs);
handles.normalized_signal.XTickLabel= axis_labels;
handles.filtered_signal.XTick = window_ind:fs:(window_ind + window_len * fs);   
handles.filtered_signal.XTickLabel= axis_labels;


function graphed = graph_function(handles, ask_questions, graph_data, config_axes, graph_spindles, graph_labels, onoff_prim, onoff_revise)
global data_n data_filt data_plot data_filt_p window_ind window_len fs plotted recorded spindle_points label_name duration_trace start_time num_comp name_comp spindle_comp color_plot spindle_plot_n spindle_plot_f comp_handles fname_edf;
if ask_questions
    def_answers = {'SpindleScorer', fname_edf, '3600', '0', '20'};
    temp = false;
    while (temp ~= true)
        temp_answers = inputdlg({'Name/Identification Information for Scorer:', 'Name of File:', 'Duration of Trace Examined (seconds):', 'Start timestamp (seconds):', 'Length of One Window (seconds):'}, 'InformationInput', 1, def_answers);
        temp = true;
        if (isempty(temp_answers) == false)
            for i=3:5
                if isnan(str2double(temp_answers{i})) == true
                    temp = false;
                end
            end
        end
    end
    if isempty(temp_answers) == true
        temp_answers = def_answers;
    end
    window_ind = str2double(temp_answers{4}) * fs;
    start_time = str2double(temp_answers{4});
    window_len = str2double(temp_answers{5});
    duration_trace = min(str2double(temp_answers{3}), max(size(data_n)) / fs);
    label_name = temp_answers{1};
    fname_edf = temp_answers{2};
end

if graph_data
    data_plot = data_n(start_time * fs + 1: (start_time + duration_trace) * fs);
    data_filt_p = data_filt(start_time * fs + 1: (start_time + duration_trace) * fs);
    data_x = start_time * fs + 1:1:(start_time + duration_trace) * fs;

    axes(handles.normalized_signal);
    set(gcf, 'NextPlot','add');
    set(gca, 'NextPlot', 'replacechildren');
    norplot = plot(data_x, data_plot, 'r');
    set(norplot, 'HitTest', 'off');
    tempax = gca;
    tempax.YLabel.String = 'Amplitude/uv';
    title(tempax, 'Normalized Signal');
    drawnow;

    axes(handles.filtered_signal);
    set(gca, 'NextPlot', 'replacechildren');
    filplot = plot(data_x, data_filt_p, 'b');
    set(filplot, 'HitTest', 'off');
    tempax = gca;
    tempax.XLabel.String = 'Time/s';
    tempax.YLabel.String = 'Amplitude/uv';
    title(tempax, 'Filtered Signal');
    drawnow;
    plotted = true;
end

if config_axes
    set(handles.load_info, 'String', strcat(num2str(window_ind), '-', num2str(window_ind + window_len * fs), ' frames'));
    drawnow;
    recorded = true;

    linkaxes([handles.normalized_signal,handles.filtered_signal],'x');
    handles.normalized_signal.XLim = [window_ind, window_ind + window_len * fs];
    handles.normalized_signal.YLim = [min(data_plot), max(data_plot)];
    handles.filtered_signal.YLim = [min(data_filt_p), max(data_filt_p)];
    handles.normalized_signal.XTick = window_ind:fs:(window_ind + window_len * fs);
    handles.filtered_signal.XTick = window_ind:fs:(window_ind + window_len * fs);     
    axis_labels = cell(1, window_len + 1);
    for i = 0:window_len
        axis_labels{i + 1} = (window_ind / fs) + i;
    end
    handles.normalized_signal.XTickLabel= axis_labels;
    handles.filtered_signal.XTickLabel= axis_labels;
    handles.normalized_signal.XLimMode = 'manual';
    handles.normalized_signal.XTickMode = 'manual';
    handles.normalized_signal.YLimMode = 'manual';
    handles.filtered_signal.YLimMode = 'manual';
end

if graph_spindles && size(spindle_points, 1) > 0
    spindle_plot_n = [];
    spindle_plot_f = [];

    axes(handles.normalized_signal);
    set(gca, 'NextPlot', 'add');
    for i = 1:size(spindle_points, 1)
        spindle_plot_n = cat(1, spindle_plot_n, plot([spindle_points(i, 1) spindle_points(i, 1)],[handles.normalized_signal.YLim(1) handles.normalized_signal.YLim(2)],'k-'));
        spindle_plot_n = cat(1, spindle_plot_n, plot([spindle_points(i, 2) spindle_points(i, 2)],[handles.normalized_signal.YLim(1) handles.normalized_signal.YLim(2)],'k:'));
    end
    axes(handles.filtered_signal);
    set(gca, 'NextPlot', 'add');
    for i = 1:size(spindle_points, 1)
        spindle_plot_f = cat(1, spindle_plot_f, plot([spindle_points(i, 1) spindle_points(i, 1)],[handles.filtered_signal.YLim(1) handles.filtered_signal.YLim(2)],'k-'));
        spindle_plot_f = cat(1, spindle_plot_f, plot([spindle_points(i, 2) spindle_points(i, 2)],[handles.filtered_signal.YLim(1) handles.filtered_signal.YLim(2)],'k:'));
    end
end

if graph_labels && num_comp > 0
    comp_handles = zeros(1, num_comp);
    for i = 1:num_comp
        axes(handles.normalized_signal);

        % New visualization
        spacing_height_lines = 20;
        height_line = handles.normalized_signal.YLim(1) + (spacing_height_lines * i);
        % End new visualization

        set(gca, 'NextPlot', 'add');
        for j = 1:size(spindle_comp{i}, 1)
            if j == 1
                comp_handles(i) = plot([spindle_comp{i}(j, 1) spindle_comp{i}(j, 2)],[height_line height_line], '-', 'LineWidth', 1, 'Color', color_plot{i});
            else
                plot([spindle_comp{i}(j, 1) spindle_comp{i}(j, 2)],[height_line height_line], '-', 'LineWidth', 1, 'Color', color_plot{i});
            end
        end
        axes(handles.filtered_signal);
        set(gca, 'NextPlot', 'add');
        for j = 1:size(spindle_comp{i}, 1)
            plot([spindle_comp{i}(j, 1) spindle_comp{i}(j, 1)],[handles.filtered_signal.YLim(1) handles.filtered_signal.YLim(2)], '-', 'LineWidth', 1, 'Color', color_plot{i});
            plot([spindle_comp{i}(j, 2) spindle_comp{i}(j, 2)],[handles.filtered_signal.YLim(1) handles.filtered_signal.YLim(2)], ':', 'LineWidth', 1, 'Color', color_plot{i});
        end
    end
    axes(handles.normalized_signal);
    legend(handles.legend_axes, comp_handles, name_comp);
    set(handles.legend_axes, 'Visible', 'on');
    
    [~, count] = count_win_spindles(num_comp, spindle_comp, window_ind, window_len, fs);
    set(handles.win_spindle_ct, 'String', strcat('Window Spindle Count: ', num2str(count)));
    drawnow;
end

set(handles.record_spindle, 'Visible', onoff_prim);
set(handles.del_spindle, 'Visible', onoff_prim);
set(handles.prev_win, 'Visible', onoff_prim);
set(handles.next_win, 'Visible', onoff_prim);
set(handles.text8, 'Visible', onoff_prim);
set(handles.text7, 'Visible', onoff_prim);
set(handles.goto_sec, 'Visible', onoff_prim);
set(handles.axes_scales, 'Visible', onoff_prim);
set(handles.slide_ctrl, 'Visible', onoff_prim);
set(handles.slide_ctrl, 'Max', (max(size(data_plot)) - window_len * fs));
set(handles.slide_ctrl, 'SliderStep', [2 * fs / max(size(data_plot)), 5 * fs / max(size(data_plot))]);

set(handles.accept_spindle, 'Visible', onoff_revise);
set(handles.reject_spindle, 'Visible', onoff_revise);
set(handles.update_spindle, 'Visible', onoff_revise);
set(handles.break_spindle, 'Visible', onoff_revise);
graphed = true;
%% ===========End helper functions===========

% --- Executes on button press in load_data.
function load_data_Callback(hObject, eventdata, handles)
% hObject    handle to load_data (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%% Initialization
% clear global;
global epoch f fs nfft win ch loaded_data sig_band loaded_scores loaded_spindles plotted record fname_edf fname_edf_full fpath_edf spindle_points num_comp spindle_comp spindle_plot_n spindle_plot_f name_comp color_plot comp_handles;
%% Get files
set(handles.rec_status, 'String', '');
drawnow;
set(handles.load_info, 'String', 'Reading edf files');
drawnow;
[fname_edf_full, fpath_edf] = uigetfile('*.edf');
[hdr, record] = edfread([fpath_edf fname_edf_full]);
fname_edf = fname_edf_full(1:7);
fs = hdr.frequency(1); % 1000 (Hz) = frequency of data sampling
win = fs*2;
nfft = 2^nextpow2(win); % 2048
f = fs/2*linspace(0,1,nfft/2+1);
sig_band = find(f>10,1):find(f<16,1,'last');
epoch = 10;
ch = 1;
cla(handles.normalized_signal);
cla(handles.filtered_signal);
loaded_data = true;
loaded_scores = false;
loaded_spindles = false;
plotted = false;
spindle_points = [];
spindle_comp = [];
num_comp = 0;
name_comp = [];
color_plot = [];
comp_handles = [];
spindle_plot_n = [];
spindle_plot_f = [];

cla(handles.legend_axes);
set(handles.legend_axes, 'Visible', 'off');
set(handles.record_spindle, 'Visible', 'off');
set(handles.del_spindle, 'Visible', 'off');
set(handles.prev_win, 'Visible', 'off');
set(handles.next_win, 'Visible', 'off');
set(handles.text8, 'Visible', 'off');
set(handles.text7, 'Visible', 'off');
set(handles.goto_sec, 'Visible', 'off');
set(handles.slide_ctrl, 'Visible', 'off');
set(handles.axes_scales, 'Visible', 'off');
set(handles.load_info, 'String', 'Raw data read. Load scores/spindles to continue.');
set(handles.spindle_count, 'String', '');
drawnow;

% --- Executes on button press in load_scores.
function load_scores_Callback(hObject, eventdata, handles)
% hObject    handle to load_scores (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global data_n data_filt data_n_segmented ind_NREM epoch fs loaded_data fname_s fpath_s loaded_scores record ch spindle_points num_comp spindle_comp;
load BP8_18Hz;
if loaded_data ~= true
    errordlg('Load data first.', 'NullPointerException');
else
    set(handles.load_info, 'String', 'Reading txt files');
    drawnow;
    [fname_s, fpath_s] = uigetfile('*.txt');
    tmp_fname = [fpath_s fname_s];
    fid = fopen(tmp_fname);
    fgets(fid);
    tline = fgets(fid); 
    count = 0;
    while tline ~= -1
        tmp = strsplit(tline,',');
        count = count + 1;
        tmp_score(count) = str2num(tmp{end-1});
        tline = fgets(fid);
    end
    fclose(fid);

    %% Read scores
    ind = find(tmp_score == 2);     % NREM
    epoch_NREM = find_seg_2(ind');
    ind_NREM(:,1) = (epoch_NREM(:,1)-1)*epoch*fs+1;
    ind_NREM(:,2) = (epoch_NREM(:,2))*epoch*fs;

    %% Stitch NREM, save data
    set(handles.load_info, 'String', 'Normalizing NREM segments');
    drawnow;
    data_n = zeros(1,1);
    data_n_segmented = cell(1, size(ind_NREM, 1));
    seg_start = 1;
    for i=1:size(ind_NREM,1)
        
        ind = ind_NREM(i,1):ind_NREM(i,2);
        data = record(ch,ind);
        processed_data = pre_process(data);

        data_n_segmented{i} = processed_data;

        data_n(seg_start: (seg_start + length(data) - 1)) = processed_data; % normalize impendance level
        seg_start = seg_start + length(data);
        set(handles.load_info, 'String', strcat('Normalizing NREM segments: ', num2str(i), '/', num2str(size(ind_NREM, 1))));
        drawnow;
    end
    set(handles.load_info, 'String', 'Filtering EEG data');
    drawnow;
    data_filt = filtfilt(BP8_18Hz.tf.num, BP8_18Hz.tf.den, data_n);
    
    loaded_scores = true;
    set(handles.load_info, 'String', 'Plot or Load Spindles to continue.');
    drawnow;
    
    spindle_points = [];
    spindle_comp = [];
    num_comp = 0;
    set(handles.spindle_count, 'String', '');
    drawnow;
end

% --- Executes on button press in load_spindles.
function load_spindles_Callback(hObject, eventdata, handles)
% hObject    handle to load_spindles (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global data_n data_filt loaded_data loaded_scores loaded_spindles record ch spindle_points spindle_plot_n spindle_plot_f;
load BP8_18Hz;
if loaded_data ~= true
    errordlg('Load data first.', 'NullPointerException');
else
    spindle_plot_n = [];
    spindle_plot_f = [];
    spindle_points = [];
    if (loaded_scores ~= true && (strcmp(questdlg('Process unscored data trace? The timestamp may be different from that of scored data.', 'ConfirmProcess', 'No', 'Yes', 'No'), 'Yes') == true))
        %% Preprocess data
        set(handles.load_info, 'String', 'Normalizing data segments');
        drawnow;
        data_n = pre_process(record(ch, :)); % normalize impendance level
        set(handles.load_info, 'String', 'Filtering EEG data');
        drawnow;
        data_filt = filtfilt(BP8_18Hz.tf.num, BP8_18Hz.tf.den, data_n);

        loaded_scores = true;
    end
    if loaded_scores == true
        set(handles.load_info, 'String', 'Reading spindle files');
        drawnow;
        [fname_sp, fpath_sp] = uigetfile('*.mat');
        load(strcat(fpath_sp, fname_sp));
        set(handles.load_info, 'String', 'Spindles loaded');
        drawnow;
    end
    set(handles.spindle_count, 'String', strcat('Spindle Count: ', num2str(size(spindle_points, 1))));
    drawnow;
    loaded_spindles = true;
end

% --- Executes on button press in load_comp.
function load_comp_Callback(hObject, eventdata, handles)
% hObject    handle to load_comp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global loaded_data spindle_comp num_comp color_plot name_comp;
if loaded_data ~= true
    errordlg('Load data first.', 'NullPointerException');
else
    set(handles.load_info, 'String', 'Reading spindle files');
    drawnow;
    [fname_sp, fpath_sp] = uigetfile('*.mat');
    temp = load(strcat(fpath_sp, fname_sp));
    num_comp = num_comp + 1;
    spindle_comp{num_comp} = temp.spindle_points;
    temp_cell{1} = spindle_comp{num_comp};
    spindle_comp{num_comp} = d2s2(temp_cell, ceil(max(max(temp_cell{1}))), 1, 1);
    set(handles.load_info, 'String', strcat('Spindles loaded; Total: ', num2str(size(spindle_comp{num_comp}, 1))));
    drawnow;
    
%     list = {'b-Blue','g-Green','r-Red','c-Cyan','m-Magenta','y-Yellow'};
%     tf = false;
%     while (tf == false)
%     [indx,tf] = listdlg('PromptString','Select a color for this set of labels:',...
%                            'SelectionMode','single',...
%                            'ListString',list);
%     end
%     color_plot(num_comp) = list{indx}(1);
    color_plot{num_comp} = uisetcolor('Choose Color');
    if (isempty(color_plot{num_comp})) || (max(color_plot{num_comp}) == 0)
        color_plot{num_comp} = [0 0 0];
    end
    name_comp{num_comp} = string(inputdlg('Enter a name for this set of labels:', 'Input', 1, {fname_sp(1:end-4)}, 'on'));
end

% --- Executes on button press in clear_comp.
function clear_comp_Callback(hObject, eventdata, handles)
% hObject    handle to clear_comp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global loaded_data spindle_comp num_comp color_plot name_comp window_ind window_len fs recorded norplot data_n filplot data_filt spindle_points spindle_plot_n spindle_plot_f plotted comp_handles;
if loaded_data ~= true
    errordlg('Load data first.', 'NullPointerException');
else
    set(handles.load_info, 'String', 'Clearing spindle files');
    drawnow;
    
    spindle_comp = [];
    color_plot = [];
    num_comp = 0;
    name_comp = [];
    comp_handles = [];
    cla(handles.normalized_signal);
    cla(handles.filtered_signal);
    cla(handles.legend_axes);
    set(handles.legend_axes, 'Visible', 'off');
    
    set(handles.load_info, 'String', 'Spindles cleared');
    drawnow;
    if(plotted)
        
        graph_function(handles, false, true, true, false, false, 'on', 'off');

        if size(spindle_points, 1) > 0
            graph_function(handles, false, false, false, true, false, 'on', 'off');
        end
    end
end

% --- Executes on button press in plot_graph.
function plot_graph_Callback(hObject, eventdata, handles)
% hObject    handle to plot_graph (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global data_n data_n_segmented data_filt ind_NREM window_ind window_len fs loaded_data loaded_scores loaded_spindles plotted recorded record ch spindle_points label_name duration_trace start_time num_comp name_comp spindle_comp color_plot spindle_plot_n spindle_plot_f comp_handles;
if (loaded_data ~= true)
    errordlg('Load data first.', 'NullPointerException');
else
    if (loaded_scores ~= true) && (loaded_spindles ~= true)
        if(strcmp(questdlg('Plot unscored data trace?', 'ConfirmPlot', 'No', 'Yes', 'No'), 'No') == true)
            return;
        end
        load BP8_18Hz;
        set(handles.load_info, 'String', 'Normalizing data segments');
        drawnow;
        data_n = pre_process(record(ch, :)); % normalize impendance level
        set(handles.load_info, 'String', 'Filtering EEG data');
        drawnow;
        data_filt = filtfilt(BP8_18Hz.tf.num, BP8_18Hz.tf.den, data_n);
    end
    
    graph_function(handles, true, true, true, true, true, 'on', 'off');
    %% Plot data
%     def_answers = {'SpindleScorer', '3600', '0', '20'};
%     temp = false;
%     while (temp ~= true)
%         temp_answers = inputdlg({'Name/Identification Information for Scorer:', 'Duration of Trace Examined (seconds):', 'Start timestamp (seconds):', 'Length of One Window (seconds):'}, 'InformationInput', 1, def_answers);
%         temp = true;
%         if (isempty(temp_answers) == false)
%             for i=2:4
%                 if isnan(str2double(temp_answers{i})) == true
%                     temp = false;
%                 end
%             end
%         end
%     end
%     if isempty(temp_answers) == true
%         temp_answers = def_answers;
%     end
%     window_ind = str2double(temp_answers{3}) * fs;
%     start_time = str2double(temp_answers{3});
%     window_len = str2double(temp_answers{4});
%     duration_trace = min(str2double(temp_answers{2}), max(size(data_n)) / fs);
%     label_name = temp_answers{1};
%     
%     data_n = data_n(window_ind + 1: duration_trace * fs);
%     data_filt = data_filt(window_ind + 1: duration_trace * fs);

    %% =================For Testing of STFT algorithm only================%
%     disp('Entering STFT.')
%     lt = 0.49;
%     ut = 0.62;
%     for i=1:size(ind_NREM, 1)
%     data_f{i} = spectrogram(data_n_segmented{i}, hamming(300), 200, [], 1000); % Check parameters for this function
%     end
%     for i=1:size(ind_NREM, 1)
%         % For each 0.5s window
%         for j=1:size(data_f{i}, 2)
%             data_f_abs{i}{j} = abs(data_f{i}(:, j));
%         end
%     end
% 
%     rel_pwr = cell(1,size(ind_NREM,1));
%     data_th = cell(1, size(ind_NREM, 1));
% 
%     total_pwr = 0;
%     count = 0;
%     for i=1:size(ind_NREM, 1)
%         total_pwr = total_pwr + sum(sum(cell2mat(data_f_abs{i})));
%         count = count + size(cell2mat(data_f_abs{i}), 2);
%     end
%     total_pwr_mean = total_pwr / count;
% 
%     for i=1:size(ind_NREM, 1)
%         for j=1:size(data_f_abs{i}, 2)
%             rel_pwr{i}(j) = sum(data_f_abs{i}{j}(6:16)) / total_pwr_mean;
%             data_th{i}((j-1)*100+1:(j-1)*100+300) = rel_pwr{i}(j);
%         end
%     end
%     data_th_avg = cell(1, size(ind_NREM, 1));
%     for i=1:size(ind_NREM, 1)
%         for j=1:size(data_f_abs{i}, 2)
%             data_th_avg{i}((j-1)*100+1:(j-1)*100+300) = mean(rel_pwr{i}(max(j-2, 1):j));
%         end
%     end
%     data_filt = cell2mat(data_th_avg);
%     data_filt = data_filt(1:3600 * fs);
%     disp('End STFT.')
    %% ================================End================================%
    
%     
%     axes(handles.normalized_signal);
%     set(gcf, 'NextPlot','add');
%     set(gca, 'NextPlot', 'replacechildren');
%     norplot = plot(data_n, 'r');
%     set(norplot, 'HitTest', 'off');
%     tempax = gca;
%     tempax.YLabel.String = 'Amplitude/uv';
%     title(tempax, 'Normalized Signal');
%     drawnow;
% 
%     axes(handles.filtered_signal);
%     set(gca, 'NextPlot', 'replacechildren');
%     filplot = plot(data_filt, 'b');
%     set(filplot, 'HitTest', 'off');
%     tempax = gca;
%     tempax.XLabel.String = 'Time/s';
%     tempax.YLabel.String = 'Amplitude/uv';
%     title(tempax, 'Filtered Signal');
%     drawnow;
%     
%     %% ========================Start STFT related=========================%
% %     set(gca, 'NextPlot', 'add');
% %     data_th_mat = cell2mat(data_th);
% %     data_th_mat = data_th_mat(1:3600 * fs);
% %     plot(data_th_mat, 'm');
% %     plot([1 length(data_th_mat)], [lt lt], 'k');
% %     plot([1 length(data_th_mat)], [ut ut], 'k');
%     %% =========================End STFT related===========================%
%     
%     set(handles.load_info, 'String', strcat(num2str(window_ind), '-', num2str(window_ind + window_len * fs), ' frames'));
%     drawnow;
%     recorded = true;
% 
%     %% Configure Axes Scales
%     linkaxes([handles.normalized_signal,handles.filtered_signal],'x');
%     handles.normalized_signal.XLim = [window_ind, window_ind + window_len * fs];
%     handles.normalized_signal.YLim = [min(data_n), max(data_n)];
%     handles.filtered_signal.YLim = [min(data_filt), max(data_filt)];
%     handles.normalized_signal.XTick = window_ind:fs:(window_ind + window_len * fs);
%     handles.filtered_signal.XTick = window_ind:fs:(window_ind + window_len * fs);     
%     axis_labels = cell(1, window_len + 1);
%     for i = 0:window_len
%         axis_labels{i + 1} = (window_ind / fs) + i;
%     end
%     handles.normalized_signal.XTickLabel= axis_labels;
%     handles.filtered_signal.XTickLabel= axis_labels;
%     handles.normalized_signal.XLimMode = 'manual';
%     handles.normalized_signal.XTickMode = 'manual';
%     handles.normalized_signal.YLimMode = 'manual';
%     handles.filtered_signal.YLimMode = 'manual';
%     plotted = true;
%     
%     if size(spindle_points, 1) > 0
%         %% Plot spindle points
%         spindle_plot_n = [];
%         spindle_plot_f = [];
%         
%         axes(handles.normalized_signal);
%         set(gca, 'NextPlot', 'add');
%         for i = 1:size(spindle_points, 1)
%             spindle_plot_n = cat(1, spindle_plot_n, plot([spindle_points(i, 1) spindle_points(i, 1)],[handles.normalized_signal.YLim(1) handles.normalized_signal.YLim(2)],'k-'));
%             spindle_plot_n = cat(1, spindle_plot_n, plot([spindle_points(i, 2) spindle_points(i, 2)],[handles.normalized_signal.YLim(1) handles.normalized_signal.YLim(2)],'k:'));
%         end
%         axes(handles.filtered_signal);
%         set(gca, 'NextPlot', 'add');
%         for i = 1:size(spindle_points, 1)
%             spindle_plot_f = cat(1, spindle_plot_f, plot([spindle_points(i, 1) spindle_points(i, 1)],[handles.filtered_signal.YLim(1) handles.filtered_signal.YLim(2)],'k-'));
%             spindle_plot_f = cat(1, spindle_plot_f, plot([spindle_points(i, 2) spindle_points(i, 2)],[handles.filtered_signal.YLim(1) handles.filtered_signal.YLim(2)],'k:'));
%         end
%     end
%     
%     if num_comp > 0
%         %% Plot Comp Points
%         comp_handles = zeros(1, num_comp);
%         for i = 1:num_comp
%             axes(handles.normalized_signal);
%             
%             % New visualization
%             spacing_height_lines = 20;
%             height_line = handles.normalized_signal.YLim(1) + (spacing_height_lines * i);
%             % End new visualization
%             
%             set(gca, 'NextPlot', 'add');
%             for j = 1:size(spindle_comp{i}, 1)
%                 if j == 1
%                     comp_handles(i) = plot([spindle_comp{i}(j, 1) spindle_comp{i}(j, 2)],[height_line height_line], strcat(color_plot(i), '-'), 'LineWidth', 1);                    
%                 else
%                     plot([spindle_comp{i}(j, 1) spindle_comp{i}(j, 2)],[height_line height_line], strcat(color_plot(i), '-'), 'LineWidth', 1);
%                 end
%             end
%             axes(handles.filtered_signal);
%             set(gca, 'NextPlot', 'add');
%             for j = 1:size(spindle_comp{i}, 1)
%                 plot([spindle_comp{i}(j, 1) spindle_comp{i}(j, 1)],[handles.filtered_signal.YLim(1) handles.filtered_signal.YLim(2)],strcat(color_plot(i), '-'), 'LineWidth', 1);
%                 plot([spindle_comp{i}(j, 2) spindle_comp{i}(j, 2)],[handles.filtered_signal.YLim(1) handles.filtered_signal.YLim(2)],strcat(color_plot(i), ':'), 'LineWidth', 1);
%             end
%         end
%         axes(handles.normalized_signal);
%         legend(handles.legend_axes, comp_handles, name_comp);
%         set(handles.legend_axes, 'Visible', 'on');
%     end
%     %% Adjust other interactive elements
%     set(handles.record_spindle, 'Visible', 'on');
%     set(handles.del_spindle, 'Visible', 'on');
%     set(handles.prev_win, 'Visible', 'on');
%     set(handles.next_win, 'Visible', 'on');
%     set(handles.text8, 'Visible', 'on');
%     set(handles.text7, 'Visible', 'on');
%     set(handles.goto_sec, 'Visible', 'on');
%     set(handles.axes_scales, 'Visible', 'on');
%     set(handles.slide_ctrl, 'Visible', 'on');
%     set(handles.slide_ctrl, 'Max', (max(size(data_n)) - window_len * fs));
%     set(handles.slide_ctrl, 'SliderStep', [2 * fs / max(size(data_n)), 5 * fs / max(size(data_n))]);
%     
%     [~, count] = count_win_spindles(num_comp, spindle_comp, window_ind, window_len, fs);
%     set(handles.win_spindle_ct, 'String', strcat('Window Spindle Count: ', num2str(count)));
%     drawnow;
%     
end

% --- Executes on button press in axes_scales.
function axes_scales_Callback(hObject, eventdata, handles)
% hObject    handle to axes_scales (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
def_answers = {num2str(handles.normalized_signal.YLim(2)), num2str(handles.normalized_signal.YLim(1)), num2str(handles.filtered_signal.YLim(2)), num2str(handles.filtered_signal.YLim(1))};
temp = false;
while (temp ~= true)
    temp_answers = inputdlg({'Normalized Signal: Y-Axis Maximum', 'Normalized Signal: Y-Axis Minimum', 'Filtered Signal: Y-Axis Maximum', 'Filtered Signal: Y-Axis Minimum'}, 'InformationInput', 1, def_answers);
    temp = true;
    if (isempty(temp_answers) == false)
        for i=1:4
            if isnan(str2double(temp_answers{i})) == true
                temp = false;
            end
        end
    end
end
if isempty(temp_answers) == true
    temp_answers = def_answers;
end
ny2 = str2double(temp_answers{1});
ny1 = str2double(temp_answers{2});
fy2 = str2double(temp_answers{3});
fy1 = str2double(temp_answers{4});
handles.normalized_signal.YLim = [ny1, ny2];
handles.filtered_signal.YLim = [fy1, fy2];

% --- Executes on button press in prev_win.
function prev_win_Callback(hObject, eventdata, handles)
% hObject    handle to prev_win (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global window_ind window_len fs plotted recorded num_comp spindle_comp;
if plotted ~= true
    errordlg('No data loaded.', 'WindowOutOfBoundsException');
% elseif window_ind == 0
%     errordlg('This is the first window.', 'WindowOutOfBoundsException');
else
    if (recorded ~= true && strcmp(questdlg('You have not recorded the spindle in this window. Continue to next window?', 'ConfirmWindowSwitch', 'No', 'Yes', 'No'), 'Yes') == true) || (recorded == true)
        set(handles.rec_status, 'String', '');
        drawnow;
        clear_window();
        window_ind = max(window_ind - window_len * fs, 0);
        axis_labels = cell(1, window_len + 1);
        
        change_view(handles, window_len, window_ind, fs);
        
        set(handles.slide_ctrl, 'Value', window_ind);
        
        set(handles.load_info, 'String', strcat(num2str(window_ind), '-', num2str(window_ind + window_len * fs), ' frames'));
        drawnow;
        recorded = true;
        [~, count] = count_win_spindles(num_comp, spindle_comp, window_ind, window_len, fs);
        set(handles.win_spindle_ct, 'String', strcat('Window Spindle Count: ', num2str(count)));
        drawnow;
    end
end

% --- Executes on button press in next_win.
function next_win_Callback(hObject, eventdata, handles)
% hObject    handle to next_win (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global data_n num_comp spindle_comp window_ind window_len fs plotted recorded;
if plotted ~= true
    errordlg('No data loaded.', 'WindowOutOfBoundsException');
% elseif window_ind + window_len * fs == 3600*fs
%     errordlg('This is the last window.', 'WindowOutOfBoundsException');
else
    if (recorded ~= true && strcmp(questdlg('You have not recorded the spindle in this window. Continue to next window?', 'ConfirmWindowSwitch', 'No', 'Yes', 'No'), 'Yes') == true) || (recorded == true)
        set(handles.rec_status, 'String', '');
        drawnow;
        clear_window();
        window_ind = min(window_ind + window_len * fs, (max(size(data_n)) - window_len * fs));

        change_view(handles, window_len, window_ind, fs);

        set(handles.slide_ctrl, 'Value', window_ind);

        set(handles.load_info, 'String', strcat(num2str(window_ind), '-', num2str(window_ind + window_len * fs), ' frames'));
        drawnow;
        recorded = true;
        [~, count] = count_win_spindles(num_comp, spindle_comp, window_ind, window_len, fs);
        set(handles.win_spindle_ct, 'String', strcat('Window Spindle Count: ', num2str(count)));
        drawnow;
    end
end

% --- Executes on change in goto_sec text.
function goto_sec_Callback(hObject, eventdata, handles)
% hObject    handle to goto_sec (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of goto_sec as text
%        str2double(get(hObject,'String')) returns contents of goto_sec as a double
global data_n num_comp spindle_comp window_ind window_len fs loaded recorded;
timestamp_goto = str2double(get(hObject, 'String'));
if loaded ~= true
    errordlg('No data loaded.', 'WindowOutOfBoundsException');
elseif isnan(timestamp_goto) == true
    errordlg('Enter a number.', 'InputNaNException');
elseif (timestamp_goto > (max(size(data_n))/fs - window_len)) || (timestamp_goto < 0)
    errordlg('Enter a number in the interval of the signal.', 'WindowOutOfBoundsException'); %Customize # seconds loaded.
else
    if (recorded ~= true && strcmp(questdlg('You have not recorded the spindle in this window. Continue to specified window?', 'ConfirmWindowSwitch', 'No', 'Yes', 'No'), 'Yes') == true) || (recorded == true)
        set(handles.rec_status, 'String', '');
        drawnow;
        clear_window();
        window_ind = timestamp_goto * fs;

        change_view(handles, window_len, window_ind, fs);

        set(handles.slide_ctrl, 'Value', window_ind);

        set(handles.load_info, 'String', strcat(num2str(window_ind), '-', num2str(window_ind + window_len * fs), ' frames'));
        drawnow;
        recorded = true;
        [~, count] = count_win_spindles(num_comp, spindle_comp, window_ind, window_len, fs);
        set(handles.win_spindle_ct, 'String', strcat('Window Spindle Count: ', num2str(count)));
        drawnow;
    end
end

% --- Executes during goto_sec creation, after setting all properties.
function goto_sec_CreateFcn(hObject, eventdata, handles)
% hObject    handle to goto_sec (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on slider movement.
function slide_ctrl_Callback(hObject, eventdata, handles)
% hObject    handle to slide_ctrl (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
global num_comp spindle_comp window_ind window_len fs loaded recorded;
if loaded ~= true
    errordlg('No data loaded.', 'WindowOutOfBoundsException');
else
    if (recorded ~= true && strcmp(questdlg('You have not recorded the spindle in this window. Continue to next window?', 'ConfirmWindowSwitch', 'No', 'Yes', 'No'), 'Yes') == true) || (recorded == true)
        set(handles.rec_status, 'String', '');
        drawnow;
        clear_window();
        window_ind = round(get(hObject, 'Value'), -3);
        hObject.Value = window_ind;

        change_view(handles, window_len, window_ind, fs);

        set(handles.load_info, 'String', strcat(num2str(window_ind), '-', num2str(window_ind + window_len * fs), ' frames'));
        drawnow;
        recorded = true;
        [~, count] = count_win_spindles(num_comp, spindle_comp, window_ind, window_len, fs);
        set(handles.win_spindle_ct, 'String', strcat('Window Spindle Count: ', num2str(count)));
        drawnow;
    end
end

% --- Executes during object creation, after setting all properties.
function slide_ctrl_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slide_ctrl (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

% --- Executes on mouse press over axes background.
function normalized_signal_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to normalized_signal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global temp_points temp_plots_n temp_plots_f recorded plotted;
if plotted ~= true
    errordlg('No data loaded.', 'WindowOutOfBoundsException');
else
    recorded = false;
    set(handles.rec_status, 'String', '');
    drawnow;
    new_point = [handles.normalized_signal.CurrentPoint(1, 1), handles.normalized_signal.CurrentPoint(1, 2)];
    temp_points = cat(1, temp_points, new_point);

    handles.normalized_signal.XLimMode = 'manual';
    handles.normalized_signal.YLimMode = 'manual';
    handles.filtered_signal.YLimMode = 'manual';

    axes(handles.normalized_signal);
    set(gca, 'NextPlot', 'add');
%     temp_plots_n = cat(1, temp_plots_n, plot([temp_points(size(temp_points, 1), 1) temp_points(size(temp_points, 1), 1)],[handles.normalized_signal.YLim(1) handles.normalized_signal.YLim(2)],'g'));
    temp_n_line = line([temp_points(size(temp_points, 1), 1) temp_points(size(temp_points, 1), 1)], [handles.normalized_signal.YLim(1) handles.normalized_signal.YLim(2)], ...
                                            'color', 'g', ...
                                            'linewidth', 1, ...
                                            'ButtonDownFcn', @startDragFcn);
    temp_plots_n = cat(1, temp_plots_n, temp_n_line);
    set(gcf, 'WindowButtonUpFcn', @stopDragFcn);

    axes(handles.filtered_signal);
    set(gca, 'NextPlot', 'add');
%     temp_plots_f = cat(1, temp_plots_f, plot([temp_points(size(temp_points, 1), 1) temp_points(size(temp_points, 1), 1)],[handles.filtered_signal.YLim(1) handles.filtered_signal.YLim(2)],'g'));
    temp_f_line = line([temp_points(size(temp_points, 1), 1) temp_points(size(temp_points, 1), 1)], [handles.filtered_signal.YLim(1) handles.filtered_signal.YLim(2)], ...
                                            'color', 'g', ...
                                            'linewidth', 1, ...
                                            'ButtonDownFcn', @startDragFcn);
    temp_plots_f = cat(1, temp_plots_f, temp_f_line);

    if size(temp_points, 1) > 2
        delete(temp_plots_n(size(temp_plots_n, 1) - 2));
        temp_plots_n = temp_plots_n(2:end, :);
        delete(temp_plots_f(size(temp_plots_f, 1) - 2));
        temp_plots_f = temp_plots_f(2:end, :);
        temp_points = temp_points(2:end, :);
    end
end

function startDragFcn(varargin)
global temp_plots_n current_line;
tp = get(gca, 'CurrentPoint');
if(size(temp_plots_n) < 2)
    current_line = 1;
elseif (abs(temp_plots_n(1).XData(1) - tp(1)) < abs(temp_plots_n(2).XData(1) - tp(1)))
    current_line = 1;
else
    current_line = 2;
end
set(gcf, 'WindowButtonMotionFcn', @draggingFcn);

function stopDragFcn(varargin)
set(gcf, 'WindowButtonMotionFcn', '');

function draggingFcn(varargin)
global temp_plots_n temp_plots_f temp_points current_line;
tp = get(gca, 'CurrentPoint');
% if(abs(temp_plots_n(1).XData(1) - tp(1)) < 200)
%     set(temp_plots_n(1), 'XData', tp(1) * [1 1]);
%     set(temp_plots_f(1), 'XData', tp(1) * [1 1]);
%     temp_points(1, 1) = tp(1);
% elseif(size(temp_plots_n) < 2)
%     
% elseif(abs(temp_plots_n(2).XData(1) - tp(1)) < 200)
%     set(temp_plots_n(2), 'XData', tp(1) * [1 1]);
%     set(temp_plots_f(2), 'XData', tp(1) * [1 1]);
%     temp_points(2, 1) = tp(1);
% else
% end
set(temp_plots_n(current_line), 'XData', tp(1) * [1 1]);
set(temp_plots_f(current_line), 'XData', tp(1) * [1 1]);
temp_points(current_line, 1) = tp(1);

% --- Executes on button press in record_spindle.
function record_spindle_Callback(hObject, eventdata, handles)
% hObject    handle to record_spindle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.rec_status, 'String', '');
drawnow;
global temp_points spindle_points recorded label_name duration_trace fname_edf start_time loaded_scores spindle_plot_n spindle_plot_f;
if size(temp_points, 1) < 2
    errordlg('Please choose the start/end points first.', 'NullPointerException');
else
    temp_cat = [temp_points(size(temp_points, 1) - 1, 1), temp_points(size(temp_points, 1), 1)];
    if temp_cat(1) > temp_cat(2)
        temp_cat = fliplr(temp_cat);
    end
    spindle_points = cat(1, spindle_points, temp_cat);
    axes(handles.normalized_signal);
    spindle_plot_n = cat(1, spindle_plot_n, plot([temp_cat(1) temp_cat(1)],[handles.normalized_signal.YLim(1) handles.normalized_signal.YLim(2)],'k-'));
    spindle_plot_n = cat(1, spindle_plot_n, plot([temp_cat(2) temp_cat(2)],[handles.normalized_signal.YLim(1) handles.normalized_signal.YLim(2)],'k:'));
    axes(handles.filtered_signal);
    spindle_plot_f = cat(1, spindle_plot_f, plot([temp_cat(1) temp_cat(1)],[handles.filtered_signal.YLim(1) handles.filtered_signal.YLim(2)],'k-'));
    spindle_plot_f = cat(1, spindle_plot_f, plot([temp_cat(2) temp_cat(2)],[handles.filtered_signal.YLim(1) handles.filtered_signal.YLim(2)],'k:'));
    
    if loaded_scores == true
        score_status = 'scored';
    else
        score_status = 'unscored';
    end
    save(strcat(fname_edf, '_', label_name, '_', num2str(duration_trace), 's_', num2str(start_time), '_start_', score_status, '_labels', '.mat'), 'spindle_points');
    set(handles.rec_status, 'String', 'Spindle recorded successully.');
    set(handles.spindle_count, 'String', strcat('Spindle Count: ', num2str(size(spindle_points, 1))));
    drawnow;
    clear_window();
    recorded = true;
end

% --- Executes on button press in del_spindle. (For future release)
function del_spindle_Callback(hObject, eventdata, handles)
% hObject    handle to del_spindle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.rec_status, 'String', '');
drawnow;
global temp_points spindle_points recorded label_name duration_trace fname_edf start_time loaded_scores spindle_plot_n spindle_plot_f;
if size(temp_points, 1) < 2
    errordlg('Please choose the start/end points first.', 'NullPointerException');
else
    temp_del = [temp_points(size(temp_points, 1) - 1, 1), temp_points(size(temp_points, 1), 1)];
    if temp_del(1) > temp_del(2)
        temp_del = fliplr(temp_del);
    end
    j = size(spindle_points, 1);
    i = 1;
   	while (i <= j)
        if spindle_points(i, 1) >= temp_del(1) && spindle_points(i, 2) <= temp_del(2)
            spindle_points(i, :) = [];
            delete(spindle_plot_n(2 * i - 1));
            delete(spindle_plot_n(2 * i));
            delete(spindle_plot_f(2 * i - 1));
            delete(spindle_plot_f(2 * i));
            spindle_plot_n(2 * i - 1, :) = [];
            spindle_plot_f(2 * i - 1, :) = [];
            spindle_plot_n(2 * i - 1, :) = [];
            spindle_plot_f(2 * i - 1, :) = [];
%             temp_plots_f = temp_plots_f(2:end, :);
            i = i - 1;
            j = j - 1;
        end
        i = i + 1;
    end
    
    if loaded_scores == true
        score_status = 'scored';
    else
        score_status = 'unscored';
    end
    save(strcat(fname_edf, '_', label_name, '_', num2str(duration_trace), 's_', num2str(start_time), '_start_', score_status, '_labels', '.mat'), 'spindle_points');
    set(handles.rec_status, 'String', 'Spindle(s) deleted successully.');
    set(handles.spindle_count, 'String', strcat('Spindle Count: ', num2str(size(spindle_points, 1))));
    drawnow;
    clear_window();
    recorded = true;
end

% --- Executes on button press in statistics.
function statistics_Callback(hObject, eventdata, handles)
% hObject    handle to statistics (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global spindle_comp name_comp spindle_events fname_edf duration_trace start_time score_status fs;

%% Overlap Calculation + TP/FP/FN/F1/Histogram
tf = false;
while (tf == false)
[stnd_ind,tf] = listdlg('PromptString', 'Choose the standard label set',...
                       'SelectionMode', 'single',...
                       'ListString', name_comp, ...
                       'ListSize', [200 100]);
end
tf = false;
while (tf == false)
[algo_ind,tf] = listdlg('PromptString', 'Choose the algorithm label set',...
                       'SelectionMode', 'single',...
                       'ListString', name_comp, ...
                       'ListSize', [200 100]);
end

% Analyze spindles using external function

if ~isempty(duration_trace)
    totduration = duration_trace * 1000;
    stats = Spindle_Analytics_Core_2(spindle_comp{stnd_ind}, spindle_comp{algo_ind}, totduration, start_time * fs);
else
    stats = Spindle_Analytics_Core_2(spindle_comp{stnd_ind}, spindle_comp{algo_ind}, ceil(max([max(max(spindle_comp{stnd_ind})), max(max(spindle_comp{algo_ind}))])), 0);
end

stats_text = {strcat("TP: ", num2str(stats.ntp), ", rate = ", num2str(stats.tpr));
        strcat("FP: ", num2str(stats.nfp), ", rate = ", num2str(stats.fpr));
        strcat("FN: ", num2str(stats.nfn), ", rate = ", num2str(stats.fnr));
        strcat("Recall: ", num2str(stats.recall));
        strcat("Precision: ", num2str(stats.precision));
        strcat("F1 Score: ", num2str(stats.f1))};

figure;
ax1 = axes('Position',[0 0 1 1],'Visible','off');
ax2 = axes('Position',[.3 .1 .6 .8]);
histogram(stats.overlap(:, 3), [0:10:100]);
axes(ax1);
text(.025, 0.6, stats_text);

% save(strcat(fname_edf, num2str(duration_trace), 's_', num2str(start_time), '_start_', score_status, '_total_labels', '.mat'), 'spindle_events')

%Saving disabled%
% save(strcat(fname_edf, num2str(duration_trace), 's_', num2str(start_time), '_start_', score_status, '_overlap_', name_comp{algo_ind}, '_wrt_', name_comp{stnd_ind}, '.mat'))

% --- Executes on button press in curr_win_stats.
function curr_win_stats_Callback(hObject, eventdata, handles)
% hObject    handle to curr_win_stats (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global spindle_comp name_comp num_comp window_ind window_len fs plotted;
if ~plotted
    errordlg('No data plotted', 'NullPointerException');
    return;
end
[spindle_in, count] = count_win_spindles(num_comp, spindle_comp, window_ind, window_len, fs);
set(handles.win_spindle_ct, 'String', strcat('Window Spindle Count: ', num2str(count)));
drawnow;
stats = cell(num_comp + 1, 1);
stats{1} = "-----Start Times-----";
for i=1:num_comp
    stats{i + 1} = "";
    stats{i + 1} = strcat(stats{i + 1}, name_comp{i}, ": ");
    for j=1:size(spindle_in{i}, 1)
        stats{i + 1} = strcat(stats{i + 1}, num2str(spindle_in{i}(j, 1) / fs), " ");
    end
end
msgbox(stats, 'Window Spindle Count');


% --- Executes on button press in review_labels.
function review_labels_Callback(hObject, eventdata, handles)
% hObject    handle to review_labels (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global spindle_comp name_comp duration_trace fs window_len start_time window_ind plotted sp_int review_labels data_n data_filt break_review idinfo spindle_revised;
if ~plotted
    errordlg('Plot data first.', 'NullPointerException');
    return;
end
spindle_revised = [];
tf = false;
while (tf == false)
[inds,tf] = listdlg('PromptString', 'Choose the label sets to combine',...
                       'ListString', name_comp, ...
                       'ListSize', [200 100]);
end
for i=1:length(inds)
    spindle_cell{i} = spindle_comp{inds(i)};
end
review_labels = d2s3(spindle_cell, duration_trace * fs, start_time * fs, 1, length(inds));

cla(handles.normalized_signal);
cla(handles.filtered_signal);
cla(handles.legend_axes);
set(handles.legend_axes, 'Visible', 'off');

graph_function(handles, false, true, true, false, false, 'off', 'on');

set(handles.load_info, 'String', '');
set(handles.spindle_count, 'String', '');
set(handles.win_spindle_ct, 'String', '');

break_review = false;
set(handles.normalized_signal, 'NextPlot', 'add');
set(handles.filtered_signal, 'NextPlot', 'add');
def_answers = {'SpindleReviewer'};
temp_answers = inputdlg({'Name/Identification Information for Labeller:'}, 'InformationInput', 1, def_answers);
if isempty(temp_answers)
    temp_answers = def_answers;
end
idinfo = temp_answers{1};
for sp_int = 1:size(review_labels, 1)
    axes(handles.normalized_signal);
    line_height = handles.normalized_signal.YLim(2);
    avg_point = mean([review_labels(sp_int, 1) review_labels(sp_int, 2)]);
    clear_window();
    handles.normalized_signal.XLim = [avg_point - window_len * fs / 2, avg_point + window_len * fs / 2];
    window_ind = handles.normalized_signal.XLim(1);
    change_view(handles, window_len, window_ind, fs);
%     graph_function(handles, false, true, true, false, false, 'off', 'on');
%     axis_labels = cell(1, window_len + 1);
%     for i = 0:window_len
%         axis_labels{i + 1} = (window_ind / fs) + i;
%     end
%     handles.normalized_signal.XTick = window_ind:fs:(window_ind + window_len * fs);
%     handles.normalized_signal.XTickLabel= axis_labels;
%     handles.filtered_signal.XTick = window_ind:fs:(window_ind + window_len * fs);
%     handles.filtered_signal.XTickLabel= axis_labels;
%     set(handles.slide_ctrl, 'Value', window_ind);
    axes(handles.normalized_signal);
    set(handles.normalized_signal, 'NextPlot', 'add');
    set(handles.filtered_signal, 'NextPlot', 'add');
    height_line = handles.normalized_signal.YLim(1)+50;
    temp_handle_2 = plot([review_labels(sp_int, 1) review_labels(sp_int, 2)],[height_line height_line], 'k-', 'LineWidth', 1);
%     temp_handle_2 = area([review_labels(sp_int, 1) review_labels(sp_int, 2)], [line_height line_height], handles.normalized_signal.YLim(1), 'FaceColor', 'y', 'FaceAlpha', 0.3, 'EdgeAlpha', 0.3, 'LineStyle', 'none');
    set(temp_handle_2, 'HitTest', 'off');
    set(handles.rec_status, 'String', strcat([num2str(sp_int) ' of ' num2str(size(review_labels, 1))]));
    uiwait(gcf);
    delete(temp_handle_2);
    if break_review
        break;
    end
end
graph_function(handles, true, true, true, false, true, 'on', 'off');
% set(handles.record_spindle, 'Visible', 'on');
% set(handles.del_spindle, 'Visible', 'on');
% set(handles.prev_win, 'Visible', 'on');
% set(handles.next_win, 'Visible', 'on');
% set(handles.text8, 'Visible', 'on');
% set(handles.text7, 'Visible', 'on');
% set(handles.goto_sec, 'Visible', 'on');
% set(handles.slide_ctrl, 'Visible', 'on');
% set(handles.axes_scales, 'Visible', 'on');
% 
% set(handles.accept_spindle, 'Visible', 'off');
% set(handles.reject_spindle, 'Visible', 'off');
% set(handles.update_spindle, 'Visible', 'off');
% set(handles.break_spindle, 'Visible', 'off');

% --- Executes on button press in accept_spindle.
function accept_spindle_Callback(hObject, eventdata, handles)
% hObject    handle to accept_spindle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global spindle_revised sp_int review_labels duration_trace start_time score_status fname_edf idinfo;
spindle_revised = cat(1, spindle_revised, [review_labels(sp_int, 1) review_labels(sp_int, 2)]);
spindle_points = spindle_revised;
save(strcat(fname_edf, '_', idinfo, '_', num2str(duration_trace), 's_', num2str(start_time), '_start_', score_status, '_labels_revised', '.mat'), 'spindle_points');
set(handles.rec_status, 'String', 'Accepted.');
uiresume(gcf);

% --- Executes on button press in reject_spindle.
function reject_spindle_Callback(hObject, eventdata, handles)
% hObject    handle to reject_spindle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.rec_status, 'String', 'Rejected.');
uiresume(gcf);

% --- Executes on button press in update_spindle.
function update_spindle_Callback(hObject, eventdata, handles)
% hObject    handle to update_spindle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global spindle_revised temp_points duration_trace start_time score_status fname_edf idinfo;
if size(temp_points, 1) < 2
    errordlg('Choose start/end points first.', 'NullPointerException');
    return;
end
temp_cat = [temp_points(size(temp_points, 1) - 1, 1), temp_points(size(temp_points, 1), 1)];
if temp_cat(1) > temp_cat(2)
    temp_cat = fliplr(temp_cat);
end
spindle_revised = cat(1, spindle_revised, temp_cat);
spindle_points = spindle_revised;
save(strcat(fname_edf, '_', idinfo, '_', num2str(duration_trace), 's_', num2str(start_time), '_start_', score_status, '_labels_revised', '.mat'), 'spindle_points');
set(handles.rec_status, 'String', 'Updated.');
uiresume(gcf);

% --- Executes on button press in break_spindle.
function break_spindle_Callback(hObject, eventdata, handles)
% hObject    handle to break_spindle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global break_review;
break_review = true;
uiresume(gcf);


% --- Executes on button press in print_axes.
function print_axes_Callback(hObject, eventdata, handles)
% hObject    handle to print_axes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
fignew = figure;
newAxes = copyobj(handles.normalized_signal, fignew);
set(newAxes, 'Position', get(groot, 'DefaultAxesPosition'));
% savefig(fignew, input('Name:'));
print('-painters', '-dmeta', strcat(input('Name:'), '.emf'))
delete(fignew);
msgbox('Current axes saved.');


% --- Executes on button press in runstft.
function runstft_Callback(hObject, eventdata, handles)
% hObject    handle to runstft (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global ch epoch fname_edf fpath_edf fname_s fpath_s;
msgbox('Select .edf recording file');
uiwait(gcf);
[fname_edf, fpath_edf] = uigetfile('*.edf');
msgbox('Select .txt scoring file');
uiwait(gcf);
[fname_s, fpath_s] = uigetfile('*.txt');
def_answers = {'8', '16', '300', '250', '0.1', '0.8', '0.01', '0.01'};
temp = false;
while (temp~= true)
    temp_answers = inputdlg({'Lower Frequency Band Range (Hz)', 'Upper Frequency Band Range (Hz)', 'Window Length (frames)', 'Window Overlap (frames)', 'Lower Bound of Thresholds (0-1)', 'Upper Bound of Thresholds (0-1)', 'Higher Threshold Step Size', 'Lower Threshold Step Size'}, 'InformationInput', 1, def_answers);
    temp = true;
    if(isempty(temp_answers) == false)
        for i=1:8
            if isnan(str2double(temp_answers{i})) == true
                temp = false;
            else
                temp_answers{i} = str2double(temp_answers{i});
            end
        end
    end
end
if (isempty(temp_answers) == true)
    for i=1:8
        temp_answers{i} = str2double(def_answers{i});
    end
end
b1 = temp_answers{1};
b2 = temp_answers{2};
wl = temp_answers{3};
wo = temp_answers{4};
ut1 = temp_answers{5};
ut2 = temp_answers{6};
utstep = temp_answers{7};
ltstep = temp_answers{8};

if (isempty(ch) || isempty(epoch))
def_answers = {'1', '10'};
temp = false;
while (temp~= true)
    temp_answers = inputdlg({'Channel', 'Epoch size'}, 'InformationInput', 1, def_answers);
    temp = true;
    if(isempty(temp_answers) == false)
        for i=1:2
            if isnan(str2double(temp_answers{i})) == true
                temp = false;
            else
                temp_answers{i} = str2double(temp_answers{i});
            end
        end
    end
end
if (isempty(temp_answers) == true)
    for i=1:2
        temp_answers{i} = str2double(def_answers{i});
    end
end
ch = temp_answers{1};
epoch = temp_answers{2};
end
Spindle_STFT_Revision_Modularized_new(ch, epoch, fpath_edf, fname_edf, fpath_s, fname_s, b1, b2, wl, wo, ut1, ut2, utstep, ltstep);


% --- Executes on button press in stftanalytics.
function stftanalytics_Callback(hObject, eventdata, handles)
% hObject    handle to stftanalytics (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
msgbox('Choose Labelling Set to Evaluate Against');
uiwait(gcf);
[fname_sp, fpath_sp] = uigetfile('*.mat');
msgbox('Choose edf file for data');
uiwait(gcf);
[fname_edf, ~] = uigetfile('*.edf');
msgbox('Choose folder for STFT results (folder directly containing STFT labels)');
uiwait(gcf);
fpath_fold = uigetdir;
def_answers = {'8', '16', '300', '250', '0.1', '0.8', '0.01', '0.01', '0.2'};
temp = false;
while (temp~= true)
    temp_answers = inputdlg({'Lower Frequency Band Range (Hz)', 'Upper Frequency Band Range (Hz)', 'Window Length (frames)', 'Window Overlap (frames)', 'Lower Bound of Threshold (0-1)', 'Upper Bound of Threshold (0-1)', 'Higher Threshold Step Size', 'Lower Threshold Step Size', 'False Negative Tolerance (0-1)'}, 'InformationInput', 1, def_answers);
    temp = true;
    if(isempty(temp_answers) == false)
        for i=1:9
            if isnan(str2double(temp_answers{i})) == true
                temp = false;
            else
                temp_answers{i} = str2double(temp_answers{i});
            end
        end
    end
end
if (isempty(temp_answers) == true)
    for i=1:9
        temp_answers{i} = str2double(def_answers{i});
    end
end
b1 = temp_answers{1};
b2 = temp_answers{2};
wl = temp_answers{3};
wo = temp_answers{4};
ut1 = temp_answers{5};
ut2 = temp_answers{6};
utstep = temp_answers{7};
ltstep = temp_answers{8};
maxfn = temp_answers{9};
Spindle_Analytics_STFT_Revision_Modularized_new(fname_sp, fpath_sp, fpath_fold, fname_edf, b1, b2, wl, wo, ut1, ut2, utstep, ltstep, maxfn);
