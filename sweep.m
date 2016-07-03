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
quant = 0.1;                          % LSB: vertical step

[t,s] = integration(t0,s0,interval,dt,t_int,quant,0);

range = (1 : (10/dt)) * dt;

%   - Divide timeline -
for k = 0 : length(s) / length(range) - 1
    
    t_div (k+1,:) =  t(  k*(length(range)) + 1 : (k+1)*(length(range)) ) ;
    s_div(k+1,:)= s(  k*(length(range)) + 1 : (k+1)*(length(range)) ) ;
    
end

%   - Peaks identification -
[kx,tx,sx, dhi,dlo, td,d, tx_N,sx_N, note_x] = signal_peaks(t_div(1,:),s_div(1,:));

%   - Sweep and minimum varaiance -
eps = 0.1;
% idx(1,1) = kx(1);
% note(1,1) = note_x(1);
% per(1,1) = tx(1);

VAR = nan(length(kx),length(kx));
kx_min = kx;

for k = 2:length(kx)
    i = 1;  
    while k > i
    VAR(i,k) = var([note_x(k) note_x(k-i)],1); 
    i = i+1;
    end
    
    [minimun(k) ind_(k)] = min(VAR(1:end,k));  
    idx(k) = (k - ind_(k));

end

for 

