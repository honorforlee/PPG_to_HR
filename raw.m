% Ivan NY HANITRA - Master thesis
%       -- Clustering peaks, discrimination, compute PPG frequency --

%   - Load file and data -
Name = '3801060_0007m';   % row 1
%Name = '3900497m';     % row 6 - BPM = 95

load(strcat(Name, '.mat'));
fid = fopen(strcat(Name, '.info'), 'rt');
fgetl(fid);
fgetl(fid);
fgetl(fid);
[interval] = sscanf(fgetl(fid), 'Sampling frequency: %f Hz  Sampling interval: %f sec');
interval = interval(2);              % data acquisition rate (interval = 1/f_spl_u = 0.5903 ms in practice)

fclose(fid);

t0 = (1:length(val)) * interval;            % timeline
s0 = val(1,1:length(val));
s0  = (s0  - mean(s0 ))/sqrt(var(s0 ));        % rescale s on 0 (standard score of signal)

%   - Timeline, noise, integration, quantization -
dt = 1/10;                           % sampling time: dt >> interval
t_int = dt * (1/3);                  % integration time: interval <= t_int < dt
quant = .1;                        % LSB: vertical step

[t,s] = integration(t0,s0,interval,dt,t_int,quant,0);

%   - Derivative - 
d = s(2:end) -  s(1:end-1);
td = (  t(2:end) +  t(1:end-1) ) / 2;

kx = d > 0;
kx = find(kx(1:end-1) & ~kx(2:end));       % k_{x}:index where d > 0; d( k_{x} + 1 ) <= 0

%   - Local maxima sx, maximum slope around sx -
[tx,sx, dhi,dlo, tx_N,sx_N, note_x] = peaks_processing(t,s,kx);

%   - Agglomerative clustering -
[clust, clust_index] = agglo_clustering(note_x,5);

%   - Remove oultiers -
kx = outlier(kx,clust_index, floor (0.05*length(kx)));

[tx,sx, dhi,dlo, tx_N,sx_N, note_x] = peaks_processing(t,s,kx);

plot(t0, s0,'k-','MarkerSize',8,'LineWidth',.5);               % siganl s
hold on
plot(t, s,'ko--','MarkerSize',10,'LineWidth',1);     % sampled signal s_n
plot(td, d,'g--','MarkerSize',10,'LineWidth',1);     % derivative of s_n
plot(tx,sx,'rd','MarkerSize',12,'LineWidth',2);
plot(tx_N,sx_N,'bd','MarkerSize',12,'LineWidth',2);
plot(kron(tx,[1 1 1]), kron(sx_N,[1 0 nan]) + kron(sx,[0 1 nan]),'r-');
hold off


%%
%plot(F);

figure(1);
plot(note_1,'.');
hold on
plot(ones(1,length(kx)) .* mean(note_1),'r-');
hold off

% figure(2);
% %subplot(2,1,2);
% plot(note_2,'.');

figure(3);
subplot(2,1,1);
plot(delta,'.');
hold on
plot(ones(1,length(kx)) .* mean(delta),'r-');
hold off

subplot(2,1,2);
plot(normlist(note_x),'.');
hold on
plot(ones(1,length(kx)) .* mean(note_x),'r-');
hold off


figure(4);
plot(t0, s0,'k-','MarkerSize',8,'LineWidth',.5);               % siganl s
hold on
plot(t, s,'ko--','MarkerSize',10,'LineWidth',1);     % sampled signal s_n
plot(td, d,'g--','MarkerSize',10,'LineWidth',1);     % derivative of s_n
plot(tx,sx,'rd','MarkerSize',12,'LineWidth',2);
plot(tx_N,sx_N,'bd','MarkerSize',12,'LineWidth',2);
plot(kron(tx,[1 1 1]), kron(sx_N,[1 0 nan]) + kron(sx,[0 1 nan]),'r-');
hold off



