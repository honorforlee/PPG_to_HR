Name = '3900497m';
load(strcat(Name, '.mat'));
fid = fopen(strcat(Name, '.info'), 'rt');
fgetl(fid); fgetl(fid); fgetl(fid);
dt0 = sscanf(fgetl(fid), 'Sampling frequency: %f Hz  Sampling interval: %f sec');
dt0 = dt0(2);
fgetl(fid); sig = fgetl(fid);

ppg = 0; ecg = 0; meas = 0;
while sig
    [k,signal,~,~,~] = strread(sig,'%d%s%f%f%s','delimiter','\t'); %#ok<DSTRRD>
    if strcmp(signal,'II');    ecg = k; end
    if strcmp(signal,'PLETH'); ppg = k; end
    if strcmp(signal,'MEAS'); meas = k; end
    sig = fgetl(fid);
end
fclose(fid);

val(isnan(val)) = [];
t0 = (1:length(val)) * dt0;            % timeline
s0 = val(ecg,1:length(val));

s0  = (s0  - mean(s0 ))/sqrt(var(s0));        % rescale s on 0 (standard score of signal)

t = []; s=[]; dt = 1;
[t0_, s0_, ~, ~] = time_div(t0,s0,dt0,t0,s0,dt0,2,1);
 
% %       -- Peaks identification --
% [kx,tx,sx, dhi,dlo, td,d, kx_n,tx_N,sx_N, note_x] = signal_peaks(t0_,s0_);
% 
% %       -- Minimum variance algorithm --
% [kx_major,tx_major,sx_major, T] = min_variance(t0_,s0_, td,d, kx,tx,sx,note_x,0.001);

addpath(fullfile(matlabroot,'/help/matlab/math/examples'));
fftgui(s0_)


%%
figure(2);
plot(t0_,s0_);

%%
% lpf = low-pass filter,  or "sacling"   (for db4: lpf = [1  3  3  1] + sqrt(3)*[ 1  1 -1 -1])
% hpf = high-pass filter, or "wavelet"   (for db4: hpf = [1 -3  3 -1] + sqrt(3)*[-1  1  1 -1])
% max_level = max level of the discrete wavelet transform
% returns  { hpf level 1, hpf level 2, ..., hpf level max_level, lpf level max_level }


% lpf = [1  3  3  1] + sqrt(3)*[ 1  1 -1 -1];
% hpf = [1 -3  3 -1] + sqrt(3)*[-1  1  1 -1];
%     
% lpf = fliplr(lpf);
% hpf = fliplr(hpf);

% w = {s0_};
% for l = 1:5
%     w{l+1} = conv( w{l} ,lpf,'valid');
%     w{l}   = conv( w{l} ,hpf,'valid');
% end