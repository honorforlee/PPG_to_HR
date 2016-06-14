% Ivan Ny Hanitra - Master thesis
% Algorithm to discriminate PPG/ECG signal peaks and recover HR

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
update_filelist(h)
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
h.t0 = (0:length(h.s0)-1) * h.dt0;
h.axes.XLim = [h.t0(1) h.t0(end)];

%Define y-axis
h.axes.YLim = [min(h.s0) max(h.s0)];
guidata(h.output, h);
h = quantisize_input(h);

function h = quantisize_input(h)        
h.dt = 1/str2double(h.edit_f_sample.String);                    % t_samp
h.t_int = str2double(h.edit_t_int.String)*h.dt;                 % integration time
h.dNds = str2double(h.edit_dNdS.String);                        % apply vertical precision (delta_samp)

subels = (1:floor(h.dt/h.dt0):length(h.t0));                    % timeline subdivision in terms of samples   
h.t = h.t0(subels);                                             % timeline with new sampling frequency                                   

% Noise
frameNoise = (0:round(h.dt/h.dt0))';
frameNoise = bsxfun(@minus, subels, frameNoise);
frameNoise_zero = find (frameNoise <= 0);
frameNoise(frameNoise_zero) = 1;

noise= random('Normal',mean(h.s0(frameNoise)),std(h.s0(frameNoise)),1,length(subels));         % Gaussian distribution (model thermal noise of finite BW)

% Integration
frameInteg = (0:floor(h.t_int/h.dt0))';
frameInteg = bsxfun(@minus, subels, frameInteg);
frameInteg_zero = find (frameInteg <= 0);
frameInteg(frameInteg_zero) = 1;                            % t_int < dt

%h.s = mean( vertcat(h.s0(frameInteg),noise) );             % sampled signal = average of Nint last values + noise during dt
h.s = mean( h.s0(frameInteg) );                             % sampled signal = average of Nint last values

h.s = h.dNds * floor( h.s / h.dNds );                       % quantization

h = grids(h);
h = process_sig(h);
plot_(h);


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


% function h = apply_filters(h)
%     [h.ft1,h.fs1] = apply_filter_( h.t , h.s , h.edit_F1.String );
%     if h.ft1
%         [h.ft2,h.fs2] = apply_filter_( h.ft1 , h.fs1 , h.edit_F2.String );
%         if h.ft2
%             [h.ft3,h.fs3] = apply_filter_( h.ft2 , h.fs2 , h.edit_F3.String );
%         end
%     end
%     [h.gt1,h.gs1] = apply_filter_( h.t , h.s , h.edit_G1.String );
%     if h.gt1
%         [h.gt2,h.gs2] = apply_filter_( h.gt1 , h.gs1 , h.edit_G2.String );
%         if h.ft2
%             [h.gt3,h.gs3] = apply_filter_( h.gt2 , h.gs2 , h.edit_G3.String );
%         end
%     end

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
    [h.tx,h.sx, h.dhi, h.dlo,h.td ,h.d ] = signal_peaks(h.t  ,h.s  );
else
    [h.tx,h.sx, h.dhi, h.dlo,h.td ,h.d ] = signal_peaks(h.ft ,h.fs );         % filter applied before derivative
end
[h.ft_,h.fs_] = apply_filter_( h.td , h.d , h.edit_F2.String );
if isempty(h.ft_)
    %         h.ft_ = [nan nan]; h.fs_ = [nan nan];
    [h.ty,h.sy,h.d2hi,h.d2lo,h.td2,h.d2] = signal_peaks(h.td ,h.d  );
else
    [h.ty,h.sy,h.d2hi,h.d2lo,h.td2,h.d2] = signal_peaks(h.ft_,h.fs_);
end

function plot_(h)
xl = h.axes.XLim; yl = h.axes.YLim;
hold off
if isempty(h.ft)
    if isempty(h.ft_)
        plot( h.axes ...
            , h.t0 , h.s0 , '-k' ...    % signal
            , h.t  , h.s  , 'ok'  ...   % sampled
            , h.td , h.d , 'x:b' ...    % first derivative
            , h.td2, h.d2, 'x:r' ...    % second derivative
            , h.xgrid,h.ygrid , ':k' ...
            , h.tx , h.sx   , 'pb' ...  % s_max = peak amplitude of s
            , h.tx , h.dhi  , '^b' ...  % local maxima around s_max
            , h.tx , h.dlo  , 'vb' ...  % local minima around s_max
            , kron(h.tx,[1 1 1]) , kron(h.sx,[0 1 nan]) , '--b' ...                               % see note_1
            , kron(h.tx,[1 1 1]) , kron(h.dlo,[1 0 nan]) + kron(h.dhi,[0 1 nan]) , '-b' ...       % link note_2
            , h.ty , h.sy   , 'pr' ...
            , h.ty , h.d2hi , '^r' ...
            , h.ty , h.d2lo , 'vr' ...
            , kron(h.ty,[1 1 1]) , kron(h.sy,[0 1 nan]) , '--r' ...
            , kron(h.ty,[1 1 1]) , kron(h.d2lo,[1 0 nan]) + kron(h.d2hi,[0 1 nan]) , '-r' ...
            )
        legend({'Signal','Sampled','D1','D2'});
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
            )
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
            )
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
            )
        legend({'Signal','Sampled','FD1','D1','FD2','D2'});                %*
    end
end
h.axes.XLim = xl; h.axes.YLim = yl;

function [tx,sx,dhi,dlo,td,d] = signal_peaks(t,s)
d  =  s(2:end) -  s(1:end-1);               % derivative
td = (  t(2:end) +  t(1:end-1) ) / 2;       % shift timeline
kx = d > 0;                                 % look at maxima
kx = find(kx(1:end-1) & ~kx(2:end));        % find derivative zero crossing: indices i where d(i+1)=0 and d(i)~=0
sx = s(kx+1);                               % asign s_max to next index
tx = td(kx) + (td(kx+1)-td(kx)) .* d(kx)./(d(kx)-d(kx+1));      % tx = td(kx) - d(kx)/D(kx) where D is derivative with respect to td
dhi = d(kx);
dlo = d(kx+1);
for k = 1:length(kx)
    i = kx(k)-1;   while i > 0         && d(i) >= dhi(k); dhi(k) = d(i); i = i-1; end    % search for local maxima at left
    i = kx(k)+2;   while i < length(d) && d(i) <= dlo(k); dlo(k) = d(i); i = i+1; end    % search for local minima at right
end

function detect_points(h) %#ok<DEFNU>
%     ECG_analysis(h); return  %TODO
%
%     if h.ecg
%         d  = h.ecg(2:end) - h.ecg(1:end-1);  td  = ( h.t0(2:end) + h.t0(1:end-1) ) / 2;
%         d2 =     d(2:end) -     d(1:end-1);  td2 = (   td(2:end) +   td(1:end-1) ) / 2;
%         [tx,sx, dhi, dlo] = signal_peaks(h.t0,h.ecg);
%         [ty,sy,d2hi,d2lo] = signal_peaks(  td,    d);
%         xl = h.axes.XLim;
%         yl = h.axes.YLim;
%         hold off
%         plot( h.axes ...
%             , h.t0 , h.ecg, '-k' ...
%             , td , d , 'x:' ...
%             , td2, d2, 'x:' ...
%             , [h.t0(1) h.t0(end)] , [0 0] , ':k' ...
%             , tx , sx   , 'pr' ...
%             , tx , dhi  , '^r' ...
%             , tx , dlo  , 'vr' ...
%             , kron(tx,[1 1 1]) , kron(sx,[0 1 nan]) , '--r' ...
%             , kron(tx,[1 1 1]) , kron( dlo,[1 0 nan]) + kron( dhi,[0 1 nan]) , '-r' ...
%             , ty , sy   , 'pb' ...
%             , ty , d2hi , '^b' ...
%             , ty , d2lo , 'vb' ...
%             , kron(ty,[1 1 1]) , kron(sy,[0 1 nan]) , '--b' ...
%             , kron(ty,[1 1 1]) , kron(d2lo,[1 0 nan]) + kron(d2hi,[0 1 nan]) , '-b' ...
%             )
%     else
d  = h.s(2:end) - h.s(1:end-1);  td  = ( h.t(2:end) + h.t(1:end-1) ) / 2;       % first derivative
d2 =   d(2:end) -   d(1:end-1);  td2 = (  td(2:end) +  td(1:end-1) ) / 2;       % second derivative
[tx,sx, dhi, dlo] = signal_peaks(h.t,h.s);                                      % detect peaks of signal
[ty,sy,d2hi,d2lo] = signal_peaks( td,  d);                                      % detect peaks of first derivative
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
    )
%     end
%     h.axes.XLim = xl;
%     h.axes.YLim = yl;
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
[tx,sx, dhi, dlo] = signal_peaks(thf,ecg_hf);
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
    )
%     hold off
%     plot( h.axes ...
%         , h.t0 , h.ecg , 'x-k' ...
%         , h.t0(k+1:end-k) , ecg_hf, '-r' ...
%         )


function callback_infile(h)  %#ok<DEFNU>
h = update_infile(h);
guidata(h.output, h);


function callback_sampling(h) %#ok<DEFNU>
h = quantisize_input(h);
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



function edit13_Callback(hObject, eventdata, handles)
% hObject    handle to edit_t_int (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_t_int as text
%        str2double(get(hObject,'String')) returns contents of edit_t_int as a double

    
% --- Executes during object creation, after setting all properties.
function edit_t_int_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_t_int (see GCBO)
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

% Hint: get(hObject,'Value') returns toggle state of checkbox5



function edit14_Callback(hObject, eventdata, handles)
% hObject    handle to edit14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit14 as text
%        str2double(get(hObject,'String')) returns contents of edit14 as a double


% --- Executes during object creation, after setting all properties.
function edit14_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox6.
function checkbox6_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox6


function edit_t_int_Callback(hObject, eventdata, handles)
% hObject    handle to edit_t_int (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_t_int as text
%        str2double(get(hObject,'String')) returns contents of edit_t_int as a double

function edit_f_sample_Callback(hObject, eventdata, handles)
% hObject    handle to edit_f_sample (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_f_sample as text
%        str2double(get(hObject,'String')) returns contents of edit_f_sample as a double


% --- Executes on button press in checkbox_dNdS.
function checkbox_dNdS_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_dNdS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_dNdS



function edit_dNdS_Callback(hObject, eventdata, handles)
% hObject    handle to edit_dNdS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_dNdS as text
%        str2double(get(hObject,'String')) returns contents of edit_dNdS as a double


% --- Executes on selection change in popupmenu_files.
function popupmenu_files_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_files (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_files contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_files


% --- Executes on button press in checkbox_ecg.
function checkbox_ecg_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_ecg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_ecg


% --- Executes on button press in checkbox_detect.
function checkbox_detect_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_detect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_detect


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function edit_F1_Callback(hObject, eventdata, handles)
% hObject    handle to edit_F1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_F1 as text
%        str2double(get(hObject,'String')) returns contents of edit_F1 as a double



function edit_F2_Callback(hObject, eventdata, handles)
% hObject    handle to edit_F2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_F2 as text
%        str2double(get(hObject,'String')) returns contents of edit_F2 as a double


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
