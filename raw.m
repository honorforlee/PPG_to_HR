Name = '3900679m';      % row 5

load(strcat(Name, '.mat'));
fid = fopen(strcat(Name, '.info'), 'rt');
fgetl(fid);
fgetl(fid);
fgetl(fid);
[interval] = sscanf(fgetl(fid), 'Sampling frequency: %f Hz  Sampling interval: %f sec');
interval = interval(2);              % data acquisition rate (interval = 1/f_spl_u = 0.5903 ms in practice)

fclose(fid);

t0 = (1:length(val)) * interval;            % timeline
s0 = val(5,1:length(val));
s0  = (s0  - mean(s0 ))/sqrt(var(s0 ));        % rescale s on 0 (standard score of signal)

%   - Timeline, noise, integration, quantization -
dt = 1/10;                           % sampling time: dt >> interval
t_int = dt * (1/3);                  % integration time: interval <= t_int < dt
quant = .1;                        % LSB: vertical step

[t,s] = integration(t0,s0,interval,dt,t_int,quant,0);

%   - Derivative -
d = s(2:end) -  s(1:end-1);
%td = (  t(2:end) +  t(1:end-1) ) / 2;      % timeline of derivative shifted by t_sample/2
td = t(2:end);

kx = d > 0;
kx = find(kx(1:end-1) & ~kx(2:end));       % k_{x}:index where d > 0; d( k_{x} + 1 ) <= 0

%   - Local maxima sx, maximum slope around sx -
[tx,sx, dhi,dlo, tx_N,sx_N, note_x] = peaks_processing(t,s,kx);


delta = delta_tx(tx);
[T,eps,R,plot_reg] = periodicity(tx);

plot_reg;

eps_note_x = 1;
eps_per = 0.05;

med = median(note_x);
med_ = find(note_x == med);

% if var([note_x(1) note_x(2)]) < 1
%     %periodic_note_x(1) = note_x(1); periodic_note_x(2) = note_x(2);
%     for k = 2 : length(kx)-1
%     var( [tx(k+1)-tx(k) tx(k)-tx(k-1)]) < 1
%         clust_per(k-1:k+1) = note_x(k-1:k+1);
%     
%         
%     end    
% else
%     periodic_note_x(1) = note_x_(1);
%     
%     
% end
