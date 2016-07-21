% Ivan Ny Hanitra - Master thesis
%       -- Algorithm to discriminate PPG/ECG signal peaks and recover HR --

% Requirements
%   - discriminate peaks in function of notes
%       note 1(pentagram): peak amplitude of sampled signal (smax)
%       note 2 (triangles): local maxima/minima around smax
%   - discriminate peaks in function of timing
%       selected peaks have the same frequency (HR)

%   - Init. - DO NOT EDIT -
function varargout = discrimination(varargin)
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @discrimination_OpeningFcn, ...
    'gui_OutputFcn',  @discrimination_OutputFcn, ...
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

%   - Executes just before discrimination is made visible -
function discrimination_OpeningFcn(hObject, ~, h, varargin)
scrsz = get(groot,'ScreenSize');
set(gcf,'Units','pixels','Position',[1 1 scrsz(3)*0.8 scrsz(4)*0.9]);
set(0,'ScreenPixelsPerInch', 95);

% scrsz = get(groot,'ScreenSize');
% set(gcf,'Units','pixels','Position',[1 .5 scrsz(3)*0.8 scrsz(4)*0.9]);
% set(0,'ScreenPixelsPerInch', 70);

fl_ecg   = {}; ppg_ecg   = {}; ecg_ecg   = {}; dt0_ecg   = {};
fl_noecg = {}; ppg_noecg = {}; ecg_noecg = {}; dt0_noecg = {};
fl_meas = {}; ppg_meas = {}; dt0_meas = {};
fl = what;
for f = fl.mat'
    f = f{1}; %#ok<FXSET>
    if exist([f(1:end-3) 'info'],'file')
        fid = fopen([f(1:end-3) 'info'], 'rt');
        fgetl(fid); fgetl(fid); fgetl(fid);
        dt0 = sscanf(fgetl(fid), 'Sampling frequency: %f Hz  Sampling interval: %f sec');
        dt0 = dt0(2);
        fgetl(fid); s = fgetl(fid);
        ppg = 0; ecg = 0; meas = 0;
        while s
            [k,signal,~,~,~] = strread(s,'%d%s%f%f%s','delimiter','\t'); %#ok<DSTRRD>
            if strcmp(signal,'II');    ecg = k; end
            if strcmp(signal,'PLETH'); ppg = k; end
            if strcmp(signal,'MEAS'); meas = k; end
            s = fgetl(fid);
        end
        fclose(fid);
        if ppg
            if ecg
                fl_ecg    = [fl_ecg    f(1:end-4)]; %#ok<AGROW>
                ppg_ecg   = [ppg_ecg   ppg];        %#ok<AGROW>
                ecg_ecg   = [ecg_ecg   ecg];        %#ok<AGROW>
                dt0_ecg   = [dt0_ecg dt0];          %#ok<AGROW>
            else
                fl_noecg  = [fl_noecg  f(1:end-4)]; %#ok<AGROW>
                ppg_noecg = [ppg_noecg ppg];        %#ok<AGROW>
                ecg_noecg = [ecg_noecg 0];          %#ok<AGROW>
                dt0_noecg = [dt0_noecg dt0];        %#ok<AGROW>
            end
        elseif meas
            fl_meas = [fl_meas f(1:end-4)];
            ppg_meas = meas;
            dt0_meas = [dt0_meas dt0];
        end
    end
end


h.f_list_database    = [fl_ecg  fl_noecg];
h.f_ppg_row = [ppg_ecg ppg_noecg];
h.f_ecg_row = [ecg_ecg ecg_noecg ];
h.f_dt0_database     = [dt0_ecg dt0_noecg];
h.nf_ecg    = length(fl_ecg);

h.f_list_meas = [fl_meas];
h.f_dt0_meas = [dt0_meas];
h.nf_meas = length(fl_meas);

update_filelist(h);
h.popupmenu_files.Value = 1;
h.ft = []; h.ft_ = []; h.fs = []; h.fs_ = [];
h.output = hObject;
h = update_infile(h);
guidata(hObject, h);

%   - Outputs from this function are returned to the command line -
function varargout = discrimination_OutputFcn(~, ~, handles)
varargout{1} = handles.output;

function file_chgd = update_filelist(h)
file_chgd = 0;
if h.checkbox_ecg.Value && ~h.checkbox_ppg_meas.Value
    if h.popupmenu_files.Value > h.nf_ecg
        h.popupmenu_files.Value = 1;
        file_chgd = 1;
    end
    h.popupmenu_files.String = h.f_list_database(1:h.nf_ecg);
    
elseif h.checkbox_ppg_meas.Value && ~h.checkbox_ecg.Value
    if h.popupmenu_files.Value > h.nf_meas
        h.popupmenu_files.Value = 1;
        file_chgd = 1;
    end
    h.popupmenu_files.String = h.f_list_meas(1:h.nf_meas);
elseif ~h.checkbox_ppg_meas.Value && ~h.checkbox_ecg.Value
    h.popupmenu_files.String = h.f_list_database;
    
else
    uiwait(msgbox('Select only one database or none by default','Warning','warn'));
end

function h = update_infile(h)
n = h.popupmenu_files.Value;
% Read data
if h.checkbox_ppg_meas.Value
    h.dt0 = h.f_dt0_meas{n};
    load([h.f_list_meas{n} '.mat']);
    Vout(isnan(Vout)) = [];
    h.s0  = Vout';
    h.s0  = (h.s0  - mean(h.s0 ))/sqrt(var(h.s0 ));
    
else
    h.dt0 = h.f_dt0_database{n};                         % timeline with the initial sampling time
    load([h.f_list_database{n} '.mat']);
    h.s0  = val(h.f_ppg_row{n}, :); %#ok<NODEF>
    %k = find(h.s0 ~= h.s0(end),1,'last');       % index k: last non zero value
    %h.s0  = h.s0(1:k);
    h.s0  = (h.s0  - mean(h.s0 ))/sqrt(var(h.s0 ));
    if h.checkbox_ecg.Value
        h.ecg = val(h.f_ecg_row{n}, :);
        h.ecg = h.ecg(1:end);
        h.ecg = (h.ecg - mean(h.ecg))/sqrt(var(h.ecg));
        h.s0 = h.ecg;
    else
        h.ecg = [];
    end
    
end

guidata(h.output, h);
h = quantize_input(h);

%   - Timeline, noise, integration, quantization -
function h = quantize_input(h)
if isempty(h.edit_f_sample.String)
    uiwait(msgbox('Fill f_sample.','Warning','warn'));
else
    h.dt   = 1/str2double(h.edit_f_sample.String);                   % apply t_sample
end

h.t_int = h.dt * h.slider_t_int.Value;

if isempty(h.edit_dNdS.String)
    uiwait(msgbox('Fill dN/dS.','Warning','warn'));
else
    h.dNds = str2double(h.edit_dNdS.String);                         % apply vertical precision (delta_samp)
end

% Define timeline
h.t0 = (1:length(h.s0)) * h.dt0;

% Integration of all signal
[h.t,h.s] = integration(h.t0,h.s0,h.dt0, h.dt,h.t_int,h.dNds,0);

% % Divide timeline
% if isempty(h.edit_frame_length.String)
%     uiwait(msgbox('Fill Frame length.','Warning','warn'));
% elseif isempty(h.edit_frame.String)
%     uiwait(msgbox('Fill Frame.','Warning','warn'));
% else
% h.frame_length = str2double(h.edit_frame_length.String);
% h.frame = str2double(h.edit_frame.String);
% end
%
% if h.frame > floor(h.t0(end)/h.frame_length)
%     h.frame = floor(h.t0(end)/h.frame_length);
% end
%
% [h.t0_ h.s0_ h.t_ h.s_] = time_div(h.t0,h.s0,h.dt0, h.t,h.s,h.dt, h.frame_length,h.frame);

% % Timeline grids
% if h.toggle_FF.Value == 0
%     h.axes.XLim = [h.t0_(1) h.t0_(end)];
%     h.axes.YLim = [min(h.s0_)-1 max(h.s0_)+1];
%
% elseif h.toggle_FF.Value == 1
%     h.axes.XLim = [h.t0(1) h.t0(end)];
%     h.axes.YLim = [min(h.s0) max(h.s0)];
% end

if isempty(h.edit_frame_length.String)
    uiwait(msgbox('Fill Frame length.','Warning','warn'));
else
    frame_length = str2double(h.edit_frame_length.String);
end

h.slider_frame.Max = floor( h.t0(end)/frame_length );
h.slider_frame.SliderStep = [1/(h.slider_frame.Max - h.slider_frame.Min) 1];
k = h.slider_frame.Value;

h.frame_init =  (k-1) * frame_length;
h.frame_end =  k * frame_length;

index0 = find(h.t0 >= h.frame_init & h.t0 <= h.frame_end);
h.t0_frame = h.t0(index0);
h.s0_frame = h.s0(index0);

% Timeline grids
if h.toggle_FF.Value == 0
    h.axes.XLim = [h.t0_frame(1) h.t0_frame(end)];
    h.axes.YLim = [min(h.s0_frame)-1 max(h.s0_frame)+1];
    
elseif h.toggle_FF.Value == 1
    h.axes.XLim = [h.t0(1) h.t0(end)];
    h.axes.YLim = [min(h.s0) max(h.s0)];
end

h = grids(h);

if h.t_int ~= 0
    if h.checkbox_clustering.Value
        h = process_clustering(h);
        plot_cluster(h);
        plot_cluster_distribution(h);
    else
        if isempty(h.edit_eps.String)
            uiwait(msgbox('Fill eps.','Warning','warn'));
        else
            h.eps = str2double(h.edit_eps.String);
        end
        
        h = process_sig(h);
        if h.warning == 0
            plot_(h);
        else
            uiwait(msgbox('No peaks detected.','Warning','warn'));
        end
    end
else
    plot_(h);
end

function [t,s] = integration(t0,s0,dt0,dt,t_int,quant,add_noise)
t = t0(1):dt:t0(end);                                   % timeline with new sampling frequency

% Noise
for k = 1:length(t)-1
    frameNoise (:,k) = [ floor(t(k)/dt0) :  floor(t(k)/dt0) + floor(dt/dt0) ];
end
noise= random('Normal',mean(s0(frameNoise)),std(s0(frameNoise)),1,length(frameNoise));                     % Gaussian distribution (model thermal noise of finite BW)

% Integration
if t_int ~=0
    
    for k = 2:length(t)
        index = min (length(    [floor( (t(k-1)+ dt - t_int)/dt0 ): floor( t(k)/ dt0 ) ]    ));
    end
    index = index-1;
    for k = 2:length(t)
        frameInteg(:,k-1) = [ floor( t(k)/ dt0 ) - index : floor( t(k)/ dt0 ) ];
        frameInteg_(:,k-1)= s0(frameInteg(:,k-1)) ;
    end
    
    if add_noise == 'noise'                               % add Gaussian noise before integration
        frameInteg_ = vertcat(frameInteg_,noise);
    elseif add_noise == 0
        frameInteg_ = frameInteg_;
    end
    
    s(1) = s0(1);
    for k = 2:length(t)
        s(k) = mean(frameInteg_(:,k-1));            % sampled signal = average of Nint last values + noise during dt
    end
    
else
    s = zeros(1,length(t));
end

s = quant * floor( s / quant );                % quantization

function h = grids(h)
h.xgrid = [ h.t0(1) h.t0(end) ];               % grid in dotted lines
h.ygrid = [ 0       0         ];

if h.checkbox_f_sample.Value
    h.xgrid = [ h.xgrid  nan  kron(h.t,[1 1 1]) ];
    h.ygrid = [ h.ygrid  nan  repmat(10*[h.axes.YLim nan],1,length(h.t)) ];
end
if h.checkbox_dNdS.Value
    y = min(h.s)-h.dNds : h.dNds : max(h.s)+h.dNds;
    h.xgrid = [ h.xgrid  nan  repmat([min(h.t) max(h.t) nan],1,length(y)) ];
    h.ygrid = [ h.ygrid  nan  kron(y,[1 1 1]) ];
end

function [ft,fs] = apply_filter_(t,s,f)
f = str2num(f); %#ok<ST2NM>
if length(f) > 0
    fs = conv( s , fliplr(f)     , 'valid' ) ;
    ft = conv( t , ones(size(f)) , 'valid' ) / length(f) ;
else
    fs = []; ft = [];
end

% function h = process_sig(h) %#ok<DEFNU>
% [h.ft ,h.fs ] = apply_filter_( h.t  , h.s , h.edit_F1.String );
% if isempty(h.ft)
%
%     [h.kx,h.tx,h.sx, h.dhi,h.dlo, h.td, h.d, h.kx_n,h.tx_N,h.sx_N, h.note_x] = signal_peaks(h.t, h.s);
% else
%     [h.kx,h.tx,h.sx, h.dhi,h.dlo, h.td, h.d, h.kx_n,h.tx_N,h.sx_N, h.note_x] =  signal_peaks(h.ft, h.fs);      %   filter applied before derivative
% end
% [h.ft_,h.fs_] = apply_filter_( h.td , h.d , h.edit_F2.String );
% if isempty(h.ft_)
%
%     [h.ky,h.ty,h.sy,h.d2hi,h.d2lo,h.td2,h.d2,h.ky_n,h.ty_N,h.sy_N,~] = signal_peaks(h.td, h.d  );
% else
%     [h.ky,h.ty,h.sy,h.d2hi,h.d2lo,h.td2,h.d2,h.ky_n,h.ty_N,h.sy_N,~] = signal_peaks(h.ft_, h.fs_);
% end

function h = process_sig(h) %#ok<DEFNU>
[h.ft ,h.fs ] = apply_filter_( h.t  , h.s , h.edit_F1.String );

if isempty(h.ft)
    
    %   - Peaks identification -
    [h.kx,h.tx,h.sx, h.dhi,h.dlo, h.td,h.d, h.kx_n,h.tx_N,h.sx_N, h.note_x] = signal_peaks(h.t,h.s);
    
    if ~h.checkbox_ecg.Value
        % Select events in the frame
        [kx_frame,tx_frame,sx_frame,note_x_frame] = frame_select(h.kx,h.tx,h.sx,h.note_x, h.frame_init,h.frame_end);
        
        %   - Minimum variance algorithm -
        [h.kx_major,h.tx_major,h.sx_major, h.T, h.warning] = min_variance(kx_frame,tx_frame,sx_frame, note_x_frame, h.eps);
    else
        h.tx_major = nan; h.sx_major = nan;             % no ecg algorithm yet
    end
    
else
    [h.kx,h.tx,h.sx, h.dhi,h.dlo, h.td, h.d, h.kx_n,h.tx_N,h.sx_N, h.note_x] =  signal_peaks(h.ft, h.fs);      %   filter applied before derivative
end
[h.ft_,h.fs_] = apply_filter_( h.td , h.d , h.edit_F2.String );
if isempty(h.ft_)
    
    [h.ky,h.ty,h.sy,h.d2hi,h.d2lo,h.td2,h.d2,h.ky_n,h.ty_N,h.sy_N,~] = signal_peaks(h.td, h.d  );
else
    [h.ky,h.ty,h.sy,h.d2hi,h.d2lo,h.td2,h.d2,h.ky_n,h.ty_N,h.sy_N,~] = signal_peaks(h.ft_, h.fs_);
end


function h = process_clustering(h)
[h.ft ,h.fs ] = apply_filter_( h.t  , h.s , h.edit_F1.String );
if isempty(h.ft)
    
    [h.tx,h.sx, h.dhi,h.dlo, h.td, h.d, h.tx_N,h.sx_N, h.note_x, h.clust_note_x, h.clust_tx, h.clust_periodicity, h.kmax, h.tx_major, h.sx_major] = events_clustering(h.t, h.s);
else
    [h.tx,h.sx, h.dhi,h.dlo, h.td, h.d, h.tx_N,h.sx_N, h.note_x, h.clust_note_x, h.clust_tx, h.clust_periodicity, h.kmax, h.tx_major, h.sx_major] =  events_clustering(h.ft, h.fs);      %   filter applied before derivative
end
[h.ft_,h.fs_] = apply_filter_( h.td , h.d , h.edit_F2.String );
if isempty(h.ft_)
    
    [h.ty,h.sy,h.d2hi,h.d2lo,h.td2,h.d2,h.ty_N,h.sy_N,~,~,~,~,~,~,~] = events_clustering(h.td, h.d  );
else
    [h.ty,h.sy,h.d2hi,h.d2lo,h.td2,h.d2,h.ty_N,h.sy_N,~,~,~,~,~,~,~] = events_clustering(h.ft_, h.fs_);
end


function [kx,tx,sx, dhi,dlo, td,d, kx_n,tx_N,sx_N, note_x] = signal_peaks(t,s)
%   - Derivative -
d = s(2:end) -  s(1:end-1);
%td = (  t(2:end) +  t(1:end-1) ) / 2;      % timeline of derivative shifted by t_sample/2
td = t(2:end);

kx = d > 0;
kx = find(kx(1:end-1) & ~kx(2:end));       % k_{x}:index where d > 0; d( k_{x} + 1 ) <= 0

%   - Local maxima sx, maximum slope around sx -
[tx,sx, dhi,dlo, kx_n,tx_N,sx_N, note_x] = peaks_processing(t,s,kx);


%   - Hierarchical clustering (agglomerative) -
function [tx,sx, dhi,dlo, td,d, tx_N,sx_N, note_x, clust_note_x, clust_tx, clust_periodicity, kmax, tx_major,sx_major ] = events_clustering(t,s)
%   - Derivative -
d = s(2:end) -  s(1:end-1);
%td = (  t(2:end) +  t(1:end-1) ) / 2;      % timeline of derivative shifted by t_sample/2
td = t(2:end);

kx = d > 0;
kx = find(kx(1:end-1) & ~kx(2:end));       % k_{x}:index where d > 0; d( k_{x} + 1 ) <= 0

%   - Local maxima sx, maximum slope around sx -
[tx,sx, dhi,dlo, kx_n,tx_N,sx_N, note_x] = peaks_processing(t,s,kx);

% Initialization
kmax_init = 6;
[clust_index,  ~,~,  ~,~,  kmax, diff] = agglo_clustering(note_x, tx, kmax_init);

% Remove oultiers
kx = outlier(kx,clust_index, floor (0.05*length(kx)));      % remove cluster containing population <= 5% length(kx)

[tx,sx, dhi,dlo, kx_n,tx_N,sx_N, note_x] = peaks_processing(t,s,kx);

if diff(2,2) >= 1    % EMPIRICAL: no clustering if 2-clustering clusters are too close
    
    % Initialization with outliers removed
    [clust_index,  clust_note_x,mean_clust,  clust_tx,clust_periodicity,  kmax, diff] = agglo_clustering(note_x, tx, kmax_init);
    
    % Search for best number of clusters
    div = 2;            % EMPIRICAL: merge clusters that are too close (    min(mean_clust difference) <= (max(mean_clust) - min(mean_clust)) / div
    while min( diff(2:end,kmax) ) <= ( max(mean_clust(:,kmax)) - min(mean_clust(:,kmax)) )/div && kmax >= 3
        
        kmax = kmax - 1;
        [clust_index,  clust_note_x,mean_clust,  clust_tx,clust_periodicity,  kmax, diff] = agglo_clustering(note_x, tx, kmax);
        
    end
    
    [~,clust_major_index] = max(mean_clust(:,kmax));
    kx_major = clust_index{clust_major_index,kmax};
    tx_major = tx(kx_major);
    sx_major = sx(kx_major);
    
else
    clust_note_x = nan;
    clust_tx = nan;
    [clust_periodicity(1),clust_periodicity(2),clust_periodicity(3)] = periodicity(tx);
    kmax = 1;
    tx_major = tx;
    sx_major = sx;
end

function plot_(h)
xl = h.axes.XLim; yl = h.axes.YLim;
hold off

if h.t_int == 0
    plot( h.axes ...
        , h.t , h.s , '-k','LineWidth',.5);
    legend('Signal');
else
    
    if isempty(h.ft)
        if isempty(h.ft_)
            if h.checkbox_signal.Value
                if h.checkbox_detect.Value == 0 && h.checkbox_detect_.Value == 0        % plot signal
                    
                    plot( h.axes ...
                        ,h.t0 , h.s0 , '-k','LineWidth',.5);
                    hold on
                    plot( h.t  , h.s  , 'ok');
                    plot( h.td , h.d , 'x:b');
                    
                    plot(h.xgrid,h.ygrid , ':k');
                    legend({'Signal','Sampled signal','First derivative'},'FontSize',8,'Orientation','Horizontal');
                    hold off
                    
                elseif h.checkbox_detect.Value == 1 && h.checkbox_detect_.Value == 0        % plot events
                    plot( h.axes ...
                        ,kron(h.tx,[1 1 1]) , kron(h.dlo,[1 0 nan]) + kron(h.dhi,[0 1 nan]), '-b');
                    hold on
                    plot( h.tx , h.sx   , 'dc','MarkerSize',10);
                    plot( h.tx,h.sx_N, 'dc','MarkerSize',10);
                    plot(kron(h.tx,[1 1 1]), kron(h.sx_N,[1 0 nan]) + kron(h.sx,[0 1 nan]),'-c');
                    plot(h.tx_major,h.sx_major, 'pr','MarkerSize',20);
                    
                    plot( h.tx , h.dhi  , '^b');
                    plot( h.tx , h.dlo  , 'vb');
                    
                    plot(h.xgrid,h.ygrid , ':k');
                    legend({'Maximum slope difference','Maxima','Minima','Peak to peak amplitude','Major peaks'},'FontSize',8,'Orientation','Horizontal');
                    hold off
                    
                elseif h.checkbox_detect_.Value == 1                                        % plot signal + events
                    plot( h.axes ...
                        ,h.t0 , h.s0 , '-k','LineWidth',.5);
                    hold on
                    plot( h.t  , h.s  , 'ok');
                    plot( h.td , h.d , 'x:b');
                    plot(kron(h.tx,[1 1 1]) , kron(h.dlo,[1 0 nan]) + kron(h.dhi,[0 1 nan]), '-b');
                    plot( h.tx , h.sx   , 'dc','MarkerSize',10);
                    plot( h.tx,h.sx_N, 'dc','MarkerSize',10);
                    plot(kron(h.tx,[1 1 1]), kron(h.sx_N,[1 0 nan]) + kron(h.sx,[0 1 nan]),'-c');
                    plot(h.tx_major,h.sx_major, 'pr','MarkerSize',20);
                    
                    plot( h.tx , h.dhi  , '^b');
                    plot( h.tx , h.dlo  , 'vb');
                    
                    plot(h.xgrid,h.ygrid , ':k');
                    legend({'Signal','Sampled signal','First derivative','Maximum slope difference','Maxima','Minima','Peak to peak amplitude','Major peaks'},'FontSize',8,'Orientation','Horizontal');
                    hold off
                    
                end
            else
                if h.checkbox_detect.Value == 0 && h.checkbox_detect_.Value == 0        % plot signal
                    
                    plot( h.axes ...
                        ,h.t  , h.s  , 'ok');
                    hold on
                    plot( h.td , h.d , 'x:b');
                    plot(h.xgrid,h.ygrid , ':k');
                    legend({'Sampled signal','First derivative'},'FontSize',8,'Orientation','Horizontal');
                    hold off
                    
                elseif h.checkbox_detect.Value == 1 && h.checkbox_detect_.Value == 0        % plot events
                    plot( h.axes ...
                        ,kron(h.tx,[1 1 1]) , kron(h.dlo,[1 0 nan]) + kron(h.dhi,[0 1 nan]), '-b');
                    hold on
                    plot( h.tx , h.sx   , 'dc','MarkerSize',10);
                    plot( h.tx,h.sx_N, 'dc','MarkerSize',10);
                    plot(kron(h.tx,[1 1 1]), kron(h.sx_N,[1 0 nan]) + kron(h.sx,[0 1 nan]),'-c');
                    plot(h.tx_major,h.sx_major, 'pr','MarkerSize',20);
                    
                    plot( h.tx , h.dhi  , '^b');
                    plot( h.tx , h.dlo  , 'vb');
                    
                    plot(h.xgrid,h.ygrid , ':k');
                    legend({'Maximum slope difference','Maxima','Minima','Peak to peak amplitude','Major peaks'},'FontSize',8,'Orientation','Horizontal');
                    hold off
                    
                elseif h.checkbox_detect_.Value == 1                                        % plot signal + events
                    plot( h.axes ...
                        ,h.t  , h.s  , 'ok');
                    hold on
                    plot( h.td , h.d , 'x:b');
                    plot(kron(h.tx,[1 1 1]) , kron(h.dlo,[1 0 nan]) + kron(h.dhi,[0 1 nan]), '-b');
                    plot( h.tx , h.sx   , 'dc','MarkerSize',10);
                    plot( h.tx,h.sx_N, 'dc','MarkerSize',10);
                    plot(kron(h.tx,[1 1 1]), kron(h.sx_N,[1 0 nan]) + kron(h.sx,[0 1 nan]),'-c');
                    plot(h.tx_major,h.sx_major, 'pr','MarkerSize',20);
                    
                    plot( h.tx , h.dhi  , '^b');
                    plot( h.tx , h.dlo  , 'vb');
                    
                    plot(h.xgrid,h.ygrid , ':k');
                    legend({'Sampled signal','First derivative','Maximum slope difference','Maxima','Minima','Peak to peak amplitude','Major peaks'},'FontSize',8,'Orientation','Horizontal');
                    hold off
                    
                end
                
            end
            
        else
            plot( h.axes ...
                ,h.t0 , h.s0 , '-k' ...  %TODO  %
                , h.t  , h.s  , 'ok'  ...
                , h.td , h.d , 'x:b' ...
                , h.ft_, h.fs_ , 'd:r' ...  % filter
                , h.td2, h.d2, 'x:r' ...
                , h.xgrid,h.ygrid , ':k' ...
                , h.tx , h.sx   , 'pb' ...
                , h.tx , h.dhi  , '^b' ...
                , h.tx , h.dlo  , 'vb' ...
                , kron(h.tx,[1 1 1]) , kron(h.sx,[0 1 nan]) , '--b' ...
                , kron(h.tx,[1 1 1]) , kron(h.dlo,[1 0 nan]) + kron(h.dhi,[0 1 nan]) , '-b' ...
                , h.ty , h.sy   , 'pr' ...
                , h.ty , h.d2hi , '^r' ...
                , h.ty , h.d2lo , 'vr' ...
                , kron(h.ty,[1 1 1]) , kron(h.sy,[0 1 nan]) , '--r' ...
                , kron(h.ty,[1 1 1]) , kron(h.d2lo,[1 0 nan]) + kron(h.d2hi,[0 1 nan]) , '-r');
            legend({'Signal','Sampled signal','D1: First derivative','F1: D1 filtered','D2: F1 derivative'});           %*
            
        end
    else
        if isempty(h.ft_)
            plot( h.axes ...
                ,h.t0 , h.s0 , '-k' ...  %TODO  %
                , h.t  , h.s  , 'ok'  ...
                , h.ft , h.fs , 'd:b'  ...
                , h.td , h.d , 'x:b' ...
                , h.td2, h.d2, 'x:r' ...
                , h.xgrid,h.ygrid , ':k' ...
                , h.tx , h.sx   , 'pb' ...
                , h.tx , h.dhi  , '^b' ...
                , h.tx , h.dlo  , 'vb' ...
                , kron(h.tx,[1 1 1]) , kron(h.sx,[0 1 nan]) , '--b' ...
                , kron(h.tx,[1 1 1]) , kron(h.dlo,[1 0 nan]) + kron(h.dhi,[0 1 nan]) , '-b' ...
                , h.ty , h.sy   , 'pr' ...
                , h.ty , h.d2hi , '^r' ...
                , h.ty , h.d2lo , 'vr' ...
                , kron(h.ty,[1 1 1]) , kron(h.sy,[0 1 nan]) , '--r' ...
                , kron(h.ty,[1 1 1]) , kron(h.d2lo,[1 0 nan]) + kron(h.d2hi,[0 1 nan]) , '-r');
            legend({'Signal','Sampled signal','F1: sampled signal filtered','D1: F1 derivative','D2: D1 derivative'});                       %*
            
        else
            plot( h.axes ...
                ,h.t0 , h.s0 , '-k' ...  %TODO  %
                , h.t  , h.s  , 'ok'  ...
                , h.ft , h.fs , 'd:b'  ...
                , h.td , h.d , 'x:b' ...
                , h.ft_, h.fs_ , 'd:r' ...
                , h.td2, h.d2, 'x:r' ...
                , h.xgrid,h.ygrid , ':k' ...
                , h.tx , h.sx   , 'pb' ...
                , h.tx , h.dhi  , '^b' ...
                , h.tx , h.dlo  , 'vb' ...
                , kron(h.tx,[1 1 1]) , kron(h.sx,[0 1 nan]) , '--b' ...
                , kron(h.tx,[1 1 1]) , kron(h.dlo,[1 0 nan]) + kron(h.dhi,[0 1 nan]) , '-b' ...
                , h.ty , h.sy   , 'pr' ...
                , h.ty , h.d2hi , '^r' ...
                , h.ty , h.d2lo , 'vr' ...
                , kron(h.ty,[1 1 1]) , kron(h.sy,[0 1 nan]) , '--r' ...
                , kron(h.ty,[1 1 1]) , kron(h.d2lo,[1 0 nan]) + kron(h.d2hi,[0 1 nan]) , '-r');
            legend({'Signal','Sampled signal','F1: sampled signal filtered','D1: F1 derivative','F2: D1 filtered','F2: F2 derivative'});                %*
            
        end
    end
end
h.axes.XLim = xl; h.axes.YLim = yl;

function plot_cluster(h)
xl = h.axes.XLim; yl = h.axes.YLim;
hold off

if isempty(h.ft)
    if isempty(h.ft_)
        if h.checkbox_detect.Value == 0 && h.checkbox_detect_.Value == 0        % plot signal
            plot( h.axes ...
                ,h.t0 , h.s0 , '-k','LineWidth',.5);
            hold on
            plot( h.t  , h.s  , 'ok');
            plot( h.td , h.d , 'x:b');
            legend({'Signal','Sampled signal','First derivative'},'FontSize',8,'Orientation','Horizontal');
            hold off
            
        elseif h.checkbox_detect.Value == 1 && h.checkbox_detect_.Value == 0  % plot events
            
            plot( kron(h.tx,[1 1 1]) , kron(h.dlo,[1 0 nan]) + kron(h.dhi,[0 1 nan]), '-c');       % link note_2
            hold on
            plot( h.tx_major , h.sx_major   , 'pk','MarkerSize',15);
            plot( h.tx , h.sx   , 'dr','MarkerSize',12);
            %plot( h.tx_N, h.sx_N , 'dm','MarkerSize',8);
            plot( h.tx,h.sx_N, 'dr','MarkerSize',12);
            plot(kron(h.tx,[1 1 1]), kron(h.sx_N,[1 0 nan]) + kron(h.sx,[0 1 nan]),'r-');
            
            
            plot( h.tx , h.dhi  , '^c');
            plot( h.tx , h.dlo  , 'vc');
            
            plot(h.xgrid,h.ygrid , ':k');
            legend({'Maximum slope difference','Major peaks','Maxima','Minima','Peak to peak amplitude'},'FontSize',8,'Orientation','Horizontal');
            hold off
            
        elseif h.checkbox_detect_.Value == 1                                % plot signal + events
            
            plot( h.axes ...
                ,h.t0 , h.s0 , '-k','LineWidth',.5);
            hold on
            plot( h.t  , h.s  , 'ok');
            plot( h.td , h.d , 'x:b');
            plot( kron(h.tx,[1 1 1]) , kron(h.dlo,[1 0 nan]) + kron(h.dhi,[0 1 nan]), '-c');       % link note_2
            
            plot( h.tx_major , h.sx_major   , 'pk','MarkerSize',15);
            
            plot( h.tx , h.sx   , 'dr','MarkerSize',12);
            %plot( h.tx_N, h.sx_N , 'dm','MarkerSize',8);
            plot( h.tx,h.sx_N, 'dr','MarkerSize',12);
            plot(kron(h.tx,[1 1 1]), kron(h.sx_N,[1 0 nan]) + kron(h.sx,[0 1 nan]),'r-');
            
            plot( h.tx , h.dhi  , '^c');
            plot( h.tx , h.dlo  , 'vc');
            
            plot(h.xgrid,h.ygrid , ':k');
            
            legend({'Signal','Sampled signal','First derivative','Maximum slope difference','Major peaks','Maxima','Minima','Peak to peak amplitude'},'FontSize',8,'Orientation','Horizontal');
            hold off
            
        end
        
    else
        plot( h.axes ...
            ,h.t0 , h.s0 , '-k' ...  %TODO  %
            , h.t  , h.s  , 'ok'  ...
            , h.td , h.d , 'x:b' ...
            , h.ft_, h.fs_ , 'd:r' ...  % filter
            , h.td2, h.d2, 'x:r' ...
            , h.xgrid,h.ygrid , ':k' ...
            , h.tx , h.sx   , 'pb' ...
            , h.tx , h.dhi  , '^b' ...
            , h.tx , h.dlo  , 'vb' ...
            , kron(h.tx,[1 1 1]) , kron(h.sx,[0 1 nan]) , '--b' ...
            , kron(h.tx,[1 1 1]) , kron(h.dlo,[1 0 nan]) + kron(h.dhi,[0 1 nan]) , '-b' ...
            , h.ty , h.sy   , 'pr' ...
            , h.ty , h.d2hi , '^r' ...
            , h.ty , h.d2lo , 'vr' ...
            , kron(h.ty,[1 1 1]) , kron(h.sy,[0 1 nan]) , '--r' ...
            , kron(h.ty,[1 1 1]) , kron(h.d2lo,[1 0 nan]) + kron(h.d2hi,[0 1 nan]) , '-r');
        legend({'Signal','Sampled signal','D1: First derivative','F1: D1 filtered','D2: F1 derivative'});           %*
        
    end
else
    if isempty(h.ft_)
        plot( h.axes ...
            ,h.t0 , h.s0 , '-k' ...  %TODO  %
            , h.t  , h.s  , 'ok'  ...
            , h.ft , h.fs , 'd:b'  ...
            , h.td , h.d , 'x:b' ...
            , h.td2, h.d2, 'x:r' ...
            , h.xgrid,h.ygrid , ':k' ...
            , h.tx , h.sx   , 'pb' ...
            , h.tx , h.dhi  , '^b' ...
            , h.tx , h.dlo  , 'vb' ...
            , kron(h.tx,[1 1 1]) , kron(h.sx,[0 1 nan]) , '--b' ...
            , kron(h.tx,[1 1 1]) , kron(h.dlo,[1 0 nan]) + kron(h.dhi,[0 1 nan]) , '-b' ...
            , h.ty , h.sy   , 'pr' ...
            , h.ty , h.d2hi , '^r' ...
            , h.ty , h.d2lo , 'vr' ...
            , kron(h.ty,[1 1 1]) , kron(h.sy,[0 1 nan]) , '--r' ...
            , kron(h.ty,[1 1 1]) , kron(h.d2lo,[1 0 nan]) + kron(h.d2hi,[0 1 nan]) , '-r');
        legend({'Signal','Sampled signal','F1: sampled signal filtered','D1: F1 derivative','D2: D1 derivative'});                       %*
        
    else
        plot( h.axes ...
            ,h.t0 , h.s0 , '-k' ...  %TODO  %
            , h.t  , h.s  , 'ok'  ...
            , h.ft , h.fs , 'd:b'  ...
            , h.td , h.d , 'x:b' ...
            , h.ft_, h.fs_ , 'd:r' ...
            , h.td2, h.d2, 'x:r' ...
            , h.xgrid,h.ygrid , ':k' ...
            , h.tx , h.sx   , 'pb' ...
            , h.tx , h.dhi  , '^b' ...
            , h.tx , h.dlo  , 'vb' ...
            , kron(h.tx,[1 1 1]) , kron(h.sx,[0 1 nan]) , '--b' ...
            , kron(h.tx,[1 1 1]) , kron(h.dlo,[1 0 nan]) + kron(h.dhi,[0 1 nan]) , '-b' ...
            , h.ty , h.sy   , 'pr' ...
            , h.ty , h.d2hi , '^r' ...
            , h.ty , h.d2lo , 'vr' ...
            , kron(h.ty,[1 1 1]) , kron(h.sy,[0 1 nan]) , '--r' ...
            , kron(h.ty,[1 1 1]) , kron(h.d2lo,[1 0 nan]) + kron(h.d2hi,[0 1 nan]) , '-r');
        legend({'Signal','Sampled signal','F1: sampled signal filtered','D1: F1 derivative','F2: D1 filtered','F2: F2 derivative'});                %*
        
    end
end
h.axes.XLim = xl; h.axes.YLim = yl;


function plot_cluster_distribution(h)
if h.t_int > 0
    if h.kmax >= 2
        
        % plot clust_note_x
        for i = 1 : h.kmax
            figure(2);
            subplot(2,1,1);
            plot(h.clust_note_x{i,h.kmax} , '.');
            hold on
        end
        Legend=cell(h.kmax,1);
        for iter=1:h.kmax
            Legend{iter}=strcat('cluster ', num2str(iter));
        end
        legend(Legend);
        title('Global note clustering');
        xlabel('k');
        ylabel('note_{x,k}, a.u');
        
        hold off
        subplot(2,1,2);
        plot( h.note_x, '.');
        
        title('Global note distribution');
        xlabel('k');
        ylabel('note_{x,k}, a.u');
        
        % plot clust_periodicity
        base_array = cellfun(@length,h.clust_tx);
        base_max = max (base_array(:,h.kmax));
        base = [1:base_max];
        
        data = nan(base_max,h.kmax);
        
        for i = 1 : h.kmax
            
            data(1:base_array(i,h.kmax),i) = h.clust_tx{i,h.kmax};
            figure(3);
            tx_disp(i) = plot(base,data(:,i),'.');
            hold on
        end
        
        Legend=cell(h.kmax,1);
        
        for iter=1:h.kmax
            Legend{iter}=strcat('cluster ', num2str(iter),': T = ', num2str(h.clust_periodicity{iter,h.kmax}(1)), '; eps = ', num2str(h.clust_periodicity{iter,h.kmax}(2)), '; R = ', num2str(h.clust_periodicity{iter,h.kmax}(3)));
        end
        legend(Legend);
        
        title('Linear regression of t_{x,k}');
        xlabel('k');
        ylabel('t_{x,k}, s');
        hold off
        
    elseif h.kmax == 1
        figure(4);
        plot( h.note_x, '.');
        
        title('Global note distribution');
        xlabel('k');
        ylabel('note_{x,k}, a.u');
        
        figure(3);
        base = [1:length(h.tx)];
        plot(base,h.tx,'.');
        Legend = strcat('T = ', num2str(h.clust_periodicity(1)), '; eps = ', num2str(h.clust_periodicity(2)), '; R = ', num2str(h.clust_periodicity(3)));
        legend(Legend);
        
        title('Linear regression of t_{x,k}');
        xlabel('k');
        ylabel('t_{x,k}, s');
    end
end


function callback_infile(h)  %#ok<DEFNU>
h = update_infile(h);
guidata(h.output, h);

function callback_grids(h) %#ok<DEFNU>
if h.checkbox_clustering.Value
    h = grids(h);
    plot_cluster(h);
    guidata(h.output, h);
else
    h = grids(h);
    plot_(h);
    guidata(h.output, h);
end

function callback_sampling(h) %#ok<DEFNU>
h = quantize_input(h);
guidata(h.output, h);

function callback_t_int(h) %#ok<DEFNU>
h = quantize_input(h);
h.value_t_int.String = h.slider_t_int.Value;
guidata(h.output, h);

function callback_frame(h) %#ok<DEFNU>
h = quantize_input(h);
h.value_frame.String = h.slider_frame.Value;
guidata(h.output, h);

function callback_filters(h) %#ok<DEFNU>
if h.checkbox_clustering.Value
    h = process_clustering(h);
    plot_cluster(h);
    guidata(h.output, h);
else
    h = process_sig(h);
    plot_(h);
    guidata(h.output, h);
end

function callback_detect(h) %#ok<DEFNU>
if h.checkbox_clustering.Value
    plot_cluster(h);
else
    plot_(h);
end

function callback_clustering(h)
if h.checkbox_clustering.Value
    h = process_clustering(h);
    plot_cluster(h);
    plot_cluster_distribution(h);
    guidata(h.output, h);
end


function callback_PPG_meas(h) %#ok<DEFNU>
if update_filelist(h)
    h = update_infile(h);
    guidata(h.output, h);
end

function callback_ECG(h) %#ok<DEFNU>
if update_filelist(h)
    h = update_infile(h);
    guidata(h.output, h);
end


function detect_points(h) %#ok<DEFNU>
d  = h.s(2:end) - h.s(1:end-1);  td  = ( h.t(2:end) + h.t(1:end-1) ) / 2;                   % first derivative
d2 =   d(2:end) -   d(1:end-1);  td2 = (  td(2:end) +  td(1:end-1) ) / 2;                   % second derivative
[kx,tx,sx, dhi, dlo,~,~,~,~,~,~] = signal_peaks(h.t,h.s);                                        % detect peaks of signal
[ky,ty,sy,d2hi,d2lo,~,~,~,~,~,~] = signal_peaks( td,  d);                                        % detect peaks of first derivative
xl = h.axes.XLim;
yl = h.axes.YLim;
hold off
plot( h.axes ...
    , h.t0 , h.s0 , '-k' ...  %TODO  %
    , h.t  , h.s  , 'ok'  ...
    , td , d , 'x:' ...
    , td2, d2, 'x:' ...
    , [h.t0(1) h.t0(end)] , [0 0] , ':k' ...
    , tx , sx   , 'pb' ...
    , tx , dhi  , '^b' ...
    , tx , dlo  , 'vb' ...
    , kron(tx,[1 1 1]) , kron(sx,[0 1 nan]) , '--b' ...
    , kron(tx,[1 1 1]) , kron( dlo,[1 0 nan]) + kron( dhi,[0 1 nan]) , '-b' ...
    , ty , sy   , 'pr' ...
    , ty , d2hi , '^r' ...
    , ty , d2lo , 'vr' ...
    , kron(ty,[1 1 1]) , kron(sy,[0 1 nan]) , '--r' ...
    , kron(ty,[1 1 1]) , kron(d2lo,[1 0 nan]) + kron(d2hi,[0 1 nan]) , '-r' ...
    );
legend({'Signal','Sampled','D1','D2'});

function ECG_analysis(h)
%     if sum(h.ecg .^ 3) < 0; h.ecg = -h.ecg; end
%     [tx,sx, dhi, dlo] = signal_peaks(h.t0,h.ecg);
%     l = sort(sx); [~,k] = max( l(2:end) - l(1:end-1) ); l = (l(k)+l(k+1))/2;
%     k = find(sx > l);
%     hold off
%     plot( h.axes ...
%         , h.t0 , h.ecg , 'x-k' ...
%         , tx(k), sx(k), 'or' ...
%         , [h.t0(1) h.t0(end)], [1 1]*l, '-b' ...
%         )
if h.edit_F1.String; s = str2double(h.edit_F1.String); else s = 1; end
s = h.dt0/s;
f = s:s:sqrt(-2*log(.01));
k = length(f);
f = exp(-.5*f.^2); sum(f)
f = [ fliplr(f) 1 f ];
if h.edit_G1.String
    s = str2double(h.edit_G1.String);
    s = h.dt0/s;
    g = s*(1:length(f)); g = g-mean(g);
    g = exp(-.5*g.^2);
    %         g = g/(1+2*sum(g));
    f = f - g;
end
f = f/sum(f);
f = [zeros(1,k) 1 zeros(1,k)] - f;
ecg_hf = conv(h.ecg,f,'valid'); thf = h.t0(k+1:end-k);
if sum(ecg_hf .^ 3) < 0; ecg_hf = -ecg_hf; secg = -1; else secg = 1; end
[kx,tx,sx, dhi, dlo,~,~,~,~,~,~] = signal_peaks(thf,ecg_hf);
l = sort(sx); [~,k] = max( l(2:end) - l(1:end-1) ); l = (l(k)+l(k+1))/2;
k = find(sx > l);
hold off
%     plot( sort(sx) ,'x'); return

plot( h.axes ...
    , h.t0 , secg*h.ecg , 'x-k' ...
    , thf , ecg_hf, '-r' ...
    , tx(k), sx(k), 'or' ...
    , [h.t0(1) h.t0(end)], [1 1]*l, '-b' ...
    , [h.t0(1) h.t0(end)], [1 1]*sqrt(var(sx)), '--b' ...
    );
