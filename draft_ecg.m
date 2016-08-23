Name = '3916979m (1)';
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

s0  = -(s0  - mean(s0 ))/sqrt(var(s0));        % rescale s on 0 (standard score of signal)

%  - Peaks identification -
[kx,tx,sx, dhi,dlo, td,d, kx_n,tx_N,sx_N, note_x] = signal_peaks(t0,s0,'ecg');

frame_init =20; frame_end = 25;

index_x = find(tx >= frame_init & tx <= frame_end);
sx_N_frame = sx_N(index_x);
dhi_frame = dhi(index_x);
dlo_frame = dlo(index_x);

index0 = find(t0 >= frame_init & t0 <= frame_end);
t0_frame = t0(index0);
s0_frame = s0(index0);

[kx_frame,tx_frame,sx_frame,note_x_frame] = frame_select(kx,tx,sx,note_x, frame_init,frame_end);

kx = kx_frame; tx = tx_frame; sx = sx_frame; note_x = note_x_frame;

[kx_major,tx_major,sx_major, T,warning] = min_variance_ecg(kx_frame,tx_frame,sx_frame, note_x_frame, 0.1);

%%
%   - Plots -
figure(2);
plot( tx_frame , sx_frame   , 'dc','MarkerSize',12);
hold on
plot(t0_frame,s0_frame,'-k');
plot(tx_frame, dhi_frame,'^b','MarkerSize',10);
plot(tx_frame, dlo_frame,'vb','MarkerSize',10);
plot( kron(tx_frame,[1 1 1]) , kron(dlo_frame,[1 0 nan]) + kron(dhi_frame,[0 1 nan]), '-b');       % link note_2
plot( tx_major , sx_major, 'pr','MarkerSize',20);
plot( tx_frame,sx_N_frame, 'dc','MarkerSize',10);
plot(kron(tx_frame,[1 1 1]), kron(sx_N_frame,[1 0 nan]) + kron(sx_frame,[0 1 nan]),'-c');
hold off


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