% Ivan NY HANITRA - Master thesis
%       -- Plot signal --

Name = '3919370m (1)';
load(strcat(Name, '.mat'));

dt0 = 0.59e-3;

val(isnan(val)) = [];
t0 = (1:length(val)) * dt0;            % timeline
s0 = val(1,1:length(val));
s0  = (s0  - mean(s0 ))/sqrt(var(s0));        % rescale s on 0 (standard score of signal)

%   - Timeline, noise, integration, quantization -
dt = 1/10;                           % sampling time: dt >> dt0
t_int = dt * (1/3);                  % integration time: dt0 <= t_int < dt
quant = 0.1;                         % LSB: vertical step

[t,s] = integration(t0,s0,dt0,dt,t_int,quant,0);

[t0_ s0_ t_ s_] = time_div(t0,s0,dt0, t,s,dt,5,2);

%  - Peaks identification -
[kx,tx,sx, dhi,dlo, td,d, kx_n,tx_N,sx_N, note_x] = signal_peaks(t_,s_);

%  - Local maxima sx, maximum slope around sx -
[tx,sx, dhi,dlo, kx_n,tx_N,sx_N, note_x] = peaks_processing(t_,s_,kx);

%   - Minimum variance algorithm -
%[kx_major,tx_major,sx_major, T, warning] = min_variance(kx,tx,sx, note_x, 0.1);
%   - Plots -

plot( tx , sx   , 'dr','MarkerSize',12);
hold on
plot(t0_,s0_,'-k','LineWidth',.1);
plot(t_,s_,'ok','LineWidth',.2);
plot(tx, dhi,'c^','MarkerSize',12);
plot(tx, dlo,'cv','MarkerSize',12);
plot( kron(tx,[1 1 1]) , kron(dlo,[1 0 nan]) + kron(dhi,[0 1 nan]), '-c');
plot( tx,sx_N, 'dr','MarkerSize',12);
plot(kron(tx,[1 1 1]), kron(sx_N,[1 0 nan]) + kron(sx,[0 1 nan]),'r-');
%plot( tx_major , sx_major, 'pr','MarkerSize',20);
hold off

