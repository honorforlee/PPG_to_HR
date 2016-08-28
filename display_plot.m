Name = '3919370m (1)';
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

if meas == 0
    val(isnan(val)) = [];
    t0 = (1:length(val)) * dt0;            % timeline
    s0 = val(ppg,1:length(val));
    
else
    Vout(isnan(Vout)) = [];
    t0 = (1:length(Vout)) * dt0;            % timeline
    s0 = Vout';
end

s0  = (s0  - mean(s0 ))/sqrt(var(s0));        % rescale s on 0 (standard score of signal)

%   - Timeline, noise, integration, quantization -
dt = 1/20;                           % sampling time: dt >> dt0
t_int = dt * (1/3);                  % integration time: dt0 <= t_int < dt
quant = 0.1;                         % LSB: vertical step

[t,s] = integration(t0,s0,dt0,dt,t_int,quant,0);

%  - Peaks identification -
[kx,tx,sx, dhi,dlo, td,d, kx_n,tx_N,sx_N, note_x] = signal_peaks(t,s);

frame_init =5; frame_end = 10;

index_x = find(tx >= frame_init & tx <= frame_end);
sx_N_frame = sx_N(index_x);
dhi_frame = dhi(index_x);
dlo_frame = dlo(index_x);

index0 = find(t0 >= frame_init & t0 <= frame_end);
t0_frame = t0(index0);
s0_frame = s0(index0);

index = find(t >= frame_init & t <= frame_end);
t_frame = t(index);
s_frame = s(index);

[kx_frame,tx_frame,sx_frame,note_x_frame] = frame_select(kx,tx,sx,note_x, frame_init,frame_end);

kx = kx_frame; tx = tx_frame; sx = sx_frame; note_x = note_x_frame;




%%
figure(1);
plot(tx_r,note_r,'pr','MarkerSize',15);
hold on
plot(tx_b,note_b,'pb','MarkerSize',15);
plot(kron(tx_r,[1 1 1]), kron(null_r,[1 0 nan]) + kron(note_r,[0 1 nan]),'-k');
plot(kron(tx_b,[1 1 1]), kron(null_b,[1 0 nan]) + kron(note_b,[0 1 nan]),'-k');

%%
%   - Plots -
figure(2);
plot( tx_frame , sx_frame   , 'dc','MarkerSize',12);
hold on
plot(t_frame,s_frame,'ok','LineWidth',.2);
plot(t0_frame,s0_frame,'-k');
plot(tx_frame, dhi_frame,'^b','MarkerSize',10);
plot(tx_frame, dlo_frame,'vb','MarkerSize',10);
plot( kron(tx_frame,[1 1 1]) , kron(dlo_frame,[1 0 nan]) + kron(dhi_frame,[0 1 nan]), '-b');       % link note_2
plot( tx_major , sx_major, 'pr','MarkerSize',20);
plot( tx_frame,sx_N_frame, 'dc','MarkerSize',10);
plot(kron(tx_frame,[1 1 1]), kron(sx_N_frame,[1 0 nan]) + kron(sx_frame,[0 1 nan]),'-c');
hold off


%%
null = zeros(1,length(tx_frame));

note_r = note_x_frame(note_x_frame>1);
note_b =  note_x_frame(note_x_frame<1);

tx_r = tx_frame(note_x_frame>1);
tx_b = tx_frame(note_x_frame<1);

null_r = zeros(1,length(note_r));
null_b = zeros(1,length(note_b));