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
fl_ecg   = {}; ppg_ecg   = {}; ecg_ecg   = {}; dt0_ecg   = {};
fl_noecg = {}; ppg_noecg = {}; ecg_noecg = {}; dt0_noecg = {};
fl = what;
for f = fl.mat'
    f = f{1}; %#ok<FXSET>
    if exist([f(1:end-3) 'info'],'file')
        fid = fopen([f(1:end-3) 'info'], 'rt');
        fgetl(fid); fgetl(fid); fgetl(fid);
        dt0 = sscanf(fgetl(fid), 'Sampling frequency: %f Hz  Sampling interval: %f sec');
        dt0 = dt0(2);
        fgetl(fid); s = fgetl(fid);
        ppg = 0; ecg = 0;
        while s
            [k,signal,~,~,~] = strread(s,'%d%s%f%f%s','delimiter','\t'); %#ok<DSTRRD>
            if strcmp(signal,'II');    ecg = k; end
            if strcmp(signal,'PLETH'); ppg = k; end
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
        end
    end
end
h.f_list    = [fl_ecg  fl_noecg];
h.f_ppg_row = [ppg_ecg ppg_noecg];
h.f_ecg_row = [ecg_ecg ecg_noecg];
h.f_dt0     = [dt0_ecg dt0_noecg];
h.nf_ecg    = length(fl_ecg);
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
if h.checkbox_ecg.Value
    if h.popupmenu_files.Value > h.nf_ecg
        h.popupmenu_files.Value = 1;
        file_chgd = 1;
    end
    h.popupmenu_files.String = h.f_list(1:h.nf_ecg);
else
    h.popupmenu_files.String = h.f_list;
end

function h = update_infile(h)
n = h.popupmenu_files.Value;
% Read data
h.dt0 = h.f_dt0{n};                         % timeline with the initial sampling time
load([h.f_list{n} '.mat']);
h.s0  = val(h.f_ppg_row{n}, :); %#ok<NODEF>
k = find(h.s0 ~= h.s0(end),1,'last');       % index k: last non zero value
h.s0  = h.s0(1:k);
h.s0  = (h.s0  - mean(h.s0 ))/sqrt(var(h.s0 ));
if h.checkbox_ecg.Value
    h.ecg = val(h.f_ecg_row{n}, :);
    h.ecg = h.ecg(1:k);
    h.ecg = (h.ecg - mean(h.ecg))/sqrt(var(h.ecg));
    h.s0 = h.ecg;
else
    h.ecg = [];
end

% Define timeline
h.t0 = (1:length(h.s0)) * h.dt0;
h.axes.XLim = [h.t0(1) h.t0(end)];

% Define y-axis
h.axes.YLim = [min(h.s0) max(h.s0)];
guidata(h.output, h);
h = quantize_input(h);


%   - timeline, noise, integration, quantization -
function h = quantize_input(h)
h.dt   = 1/str2double(h.edit_f_sample.String);                    % apply t_sample
h.t_int = h.dt * h.slider_t_int.Value;
h.dNds = str2double(h.edit_dNdS.String);                        % apply vertical precision (delta_samp)

% h.t = h.t0(1):h.dt:h.t0(end);                                   % timeline with new sampling frequency
% h.s = h.dNds * floor( interp1(h.t0,h.s0,h.t) / h.dNds );

[h.t,h.s] = integration(h.t0,h.s0,h.dt0,h.dt,h.t_int,h.dNds,0);

if h.t_int ~0
    h = grids(h);
    h = process_sig(h);
    plot_(h);
else
    h = grids(h);
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
if t_int ~0
    
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
h.xgrid = [ h.t0(1) h.t0(end) ];
h.ygrid = [ 0       0         ];

if h.checkbox_f_sample.Value                                       % ?
    h.xgrid = [ h.xgrid  nan  kron(h.t,[1 1 1]) ];
    h.ygrid = [ h.ygrid  nan  repmat(10*[h.axes.YLim nan],1,length(h.t)) ];
end
if h.checkbox_dNdS.Value                                        % ?
    y = min(h.s)-h.dNds : h.dNds : max(h.s)+h.dNds;
    h.xgrid = [ h.xgrid  nan  repmat([min(h.t) max(h.t) nan],1,length(y)) ];
    h.ygrid = [ h.ygrid  nan  kron(y,[1 1 1]) ];
end

function [ft,fs] = apply_filter_(t,s,f)         %? values entered
f = str2num(f); %#ok<ST2NM>
if length(f) > 0
    fs = conv( s , fliplr(f)     , 'valid' ) ;
    ft = conv( t , ones(size(f)) , 'valid' ) / length(f) ;
else
    fs = []; ft = [];
end

function h = process_sig(h) %#ok<DEFNU>
[h.ft ,h.fs ] = apply_filter_( h.t  , h.s , h.edit_F1.String );
if isempty(h.ft)
    %         h.ft  = [nan nan]; h.fs  = [nan nan];
    [h.tx, h.sx, h.dhi, h.dlo, h.td , h.d, h.tx_N, h.sx_N, h.note_1, h.note_2, h.delta, h.note_x, h.clust,h.F] = signal_peaks(h.t, h.s);
else
    [h.tx, h.sx, h.dhi, h.dlo, h.td , h.d,h.tx_N, h.sx_N, h.note_1, h.note_2, h.delta, h.note_x, h.clust,h.F] = signal_peaks(h.ft, h.fs );         %
    
    filter applied before derivative
end
[h.ft_,h.fs_] = apply_filter_( h.td , h.d , h.edit_F2.String );
if isempty(h.ft_)
    %         h.ft_ = [nan nan]; h.fs_ = [nan nan];
    [h.ty,h.sy,h.d2hi,h.d2lo,h.td2,h.d2,h.ty_N,h.sy_N,~,~,~,~,~,~] = signal_peaks(h.td, h.d  );
else
    [h.ty,h.sy,h.d2hi,h.d2lo,h.td2,h.d2,h.ty_N,h.sy_N,~,~,~,~,~,~] = signal_peaks(h.ft_, h.fs_);
end

function plot_(h)
xl = h.axes.XLim; yl = h.axes.YLim;
hold off

if h.t_int == 0
    plot( h.axes ...
        ,h.t , h.s , '-k','LineWidth',.5);
    legend('Signal');
else
    if isempty(h.ft)
        if isempty(h.ft_)
            plot( h.axes ...
                ,h.t0 , h.s0 , '-k','LineWidth',.5);
            hold on
            plot( h.t  , h.s  , 'ok');
            plot( h.td , h.d , 'x:b');
            plot( h.tx , h.sx   , 'pr','MarkerSize',15,'LineWidth',2);
            plot( h.tx_N, h.sx_N , 'pm','MarkerSize',15,'LineWidth',2);
            plot( kron(h.tx,[1 1 1]) , kron(h.dlo,[1 0 nan]) + kron(h.dhi,[0 1 nan]), '-c');       % link note_2
            plot(kron(h.tx,[1 1 1]), kron(h.sx_N,[1 0 nan]) + kron(h.sx,[0 1 nan]),'r-');
            
            plot( h.tx , h.dhi  , '^c');
            plot( h.tx , h.dlo  , 'vc');
            plot(h.xgrid,h.ygrid , ':k');
            
            legend({'Signal','Sampled signal','D1','Major peaks','Minima','Note_2','Peak to peak amplitude'});
            hold off
            
            %   - plot  note_x clustering -
            for i = 1 : 2
                figure(2);
                subplot(2,1,1);
                plot(h.clust{i,2} , '.');
                hold on
            end
            hold off
            subplot(2,1,2);
            plot( abs(normlist(h.note_x)), '.');
            
            
            %             figure(6);
            %             plot (h.F,'b-');
            %
            %   - plot peaks distribution -
            %             figure(1);
            %             subplot(1,1,3);
            %             plot(h.note_1);
            %             subplot(2,1,3);
            %             plot(h.note_2);
            %             subplot(3,1,3)
            %             plot(h.delta,'r.');
            
            %             xlabel('kx'); ylabel('note_x');
            %             title('Peaks note distribution');
            
        else
            plot( h.axes ...
                , h.t0 , h.s0 , '-k' ...  %TODO  %
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
                , kron(h.ty,[1 1 1]) , kron(h.d2lo,[1 0 nan]) + kron(h.d2hi,[0 1 nan]) , '-r' ...
                );
            legend({'Signal','Sampled','D1','FD2','D2'});           %*
        end
    else
        if isempty(h.ft_)
            plot( h.axes ...
                , h.t0 , h.s0 , '-k' ...  %TODO  %
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
                , kron(h.ty,[1 1 1]) , kron(h.d2lo,[1 0 nan]) + kron(h.d2hi,[0 1 nan]) , '-r' ...
                );
            legend({'Signal','Sampled','FD1','D1','D2'});                       %*
        else
            plot( h.axes ...
                , h.t0 , h.s0 , '-k' ...  %TODO  %
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
                , kron(h.ty,[1 1 1]) , kron(h.d2lo,[1 0 nan]) + kron(h.d2hi,[0 1 nan]) , '-r' ...
                );
            legend({'Signal','Sampled','FD1','D1','FD2','D2'});                %*
        end
    end
end
h.axes.XLim = xl; h.axes.YLim = yl;

function [tx,sx,dhi,dlo,td,d,tx_N,sx_N,note_1,note_2,delta,note_x,clust,F] = signal_peaks(t,s)
%   - Derivative, local maxima sx, maximum slope around sx -
d = s(2:end) -  s(1:end-1);
td = (  t(2:end) +  t(1:end-1) ) / 2;

kx = d > 0;
kx = find(kx(1:end-1) & ~kx(2:end));       % k_{x}:index where d > 0; d( k_{x} + 1 ) <= 0

sx = s(kx+1);                          % local maxima
tx = td(kx) + (td(kx+1)-td(kx)) .* d(kx)./(d(kx)-d(kx+1));      % linear interpolation of dhi and dho to get tx (@zero crossing)

dhi = d(kx);
dlo = d(kx+1);

for k = 1:length(kx)
    i = kx(k)-1;   while i > 0             && d(i) >= dhi(k); dhi(k) = d(i); i = i-1; end    % search for maximum positive slope at kx-
    i = kx(k)+2;   while i < length(d) && d(i) <= dlo(k); dlo(k) = d(i); i = i+1; end    % search for maximmum negative slope at kx+
end

kx_n = d < 0;                               % search for local minima
kx_n = find(kx_n(1:end-1) & ~kx_n(2:end));

if kx_n(1) < kx(1)
    for k = 1:length(kx)
        kx_index(k) = max( find( kx_n < kx(k) ) );
    end
    sx_N = s(kx_n( kx_index ) + 1);
    tx_N = td(kx_n( kx_index )) + (td(kx_n( kx_index )+1)-td(kx_n( kx_index ))) .* d(kx_n( kx_index ))./(d(kx_n( kx_index ))-d(kx_n( kx_index )+1));
else
    kx_index(1) = nan;
    sx_N(1) = nan;
    tx_N(1) = nan;
    
    for k = 2:length(kx)
        kx_index(k) = max( find( kx_n < kx(k) ) );
        sx_N(k) = s(kx_n( kx_index(k)) + 1);
        tx_N(k) = td(kx_n( kx_index(k) )) + (td(kx_n( kx_index(k) )+1)-td(kx_n( kx_index(k) ))) .* d(kx_n( kx_index(k) ))./(d(kx_n( kx_index(k) ))-d(kx_n( kx_index(k) )+1));
    end
    
end

%   - Peaks notation
note_1 = sx;
for k = 2:length(kx)-1
    note_1(k) = 2*sx(k) - sx(k+1) - sx(k-1);  % average peak value (doubled)
end

note_2 = dhi - dlo;                           % maximum slope difference around peak

delta = sx - sx_N;

note_3 = delta;
for k = 2:length(kx)-1
    note_3(k) = 2*delta(k) - delta(k+1) - delta(k-1);  % average peak to peak value (doubled)
end

note_x = note_3;                                % global note

%   - Hierarchical clustering according to Ward's criterion and F-statistics to evaluate best number of cluster 
k_max = 5;
clust{1} = note_x;

for k = 2:k_max
    c = clusterdata(note_x','linkage','ward','savememory','on','maxclust',k);
    
    for i = 1 : k     % inter clust
        clust_index{i,k} = find(c == i);
        clust{i,k} = note_x (clust_index{i,k});   % clust partition
        
        n = cellfun(@length,clust);
        
        %         num_F_(i) =( n(i,k) * (distance(mean(clust{i,k}), mean(note_x), 2))^2 ) / (k - 1);            % distance INTER - clust
        %
        %         for j = 1 : n(i,k)     % intra clust
        %             den_F_d(j) = distance( clust{i,k}(j), mean(clust{i,k}), 2)^2 / (length(kx) - k);        % distance INTRA - clust j
        %         end
        %
        %         den_F_(i) = sum(den_F_d);
        %         clearvars den_F_d;
        
        num_F_(i) =( n(i,k) * (distance(mean(clust{i,k}), mean(note_x), 2))^2 );
        
        for j = 1 : n(i,k)     % intra clust
            den_F_d(j) = distance( clust{i,k}(j), mean(clust{i,k}), 2)^2;        % distance INTRA - clust j
        end
        
        den_F_(i) = sum(den_F_d);
        clearvars den_F_d
        
    end
    
    num_F(k) = sum(num_F_);
    den_F(k) = sum(den_F_);
    F(k) = num_F(k) / den_F(k);             % F-statistics notation
    warning('off','all');
end


function callback_infile(h)  %#ok<DEFNU>
h = update_infile(h);
guidata(h.output, h);

function callback_sampling(h) %#ok<DEFNU>
h = quantize_input(h);
guidata(h.output, h);

function callback_filters(h) %#ok<DEFNU>
h = process_sig(h);
plot_(h);
guidata(h.output, h);

% function callback_detect(h) %#ok<DEFNU>
%     h = apply_filters(h);
%     plot_(h);
%     guidata(h.output, h);

function callback_grids(h) %#ok<DEFNU>
h = grids(h);
plot_(h);
guidata(h.output, h);

function callback_ECG(h) %#ok<DEFNU>
if update_filelist(h)
    h = update_infile(h);
    guidata(h.output, h);
end

function callback_t_int(h) %#ok<DEFNU>
h = quantize_input(h);
h.value_t_int.String = h.slider_t_int.Value;
guidata(h.output, h);

function detect_points(h) %#ok<DEFNU>
d  = h.s(2:end) - h.s(1:end-1);  td  = ( h.t(2:end) + h.t(1:end-1) ) / 2;                   % first derivative
d2 =   d(2:end) -   d(1:end-1);  td2 = (  td(2:end) +  td(1:end-1) ) / 2;                   % second derivative
[tx,sx, dhi, dlo,~,~,~,~,~,~,~,~,~,~] = signal_peaks(h.t,h.s);                                      % detect peaks of signal
[ty,sy,d2hi,d2lo,~,~,~,~,~,~,~,~,~,~] = signal_peaks( td,  d);                                      % detect peaks of first derivative
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
[tx,sx, dhi, dlo,~,~,~,~,~,~,~,~,~,~] = signal_peaks(thf,ecg_hf);
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


function edit13_Callback(hObject, eventdata, handles)
% hObject    handle to edit13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit13 as text
%        str2double(get(hObject,'String')) returns contents of edit13 as a double


% --- Executes during object creation, after setting all properties.
function edit13_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox5.
function checkbox5_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%
% Hint: get(hObject,'Value') returns toggle state of checkbox5

