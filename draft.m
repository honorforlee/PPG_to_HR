Name = '3900679m';       % row 5

load(strcat(Name, '.mat'));
fid = fopen(strcat(Name, '.info'), 'rt');
fgetl(fid);
fgetl(fid);
fgetl(fid);
[interval] = sscanf(fgetl(fid), 'Sampling frequency: %f Hz  Sampling interval: %f sec');
dt0 = interval(2);              % data acquisition rate (interval = 1/f_spl_u = 0.5903 ms in practice)

fclose(fid);

t0 = (1:length(val)) * dt0;            % timeline
s0 = val(5,1:length(val));
s0  = (s0  - mean(s0 ))/sqrt(var(s0 ));        % rescale s on 0 (standard score of signal)

%   - Timeline, noise, integration, quantization -
dt = 1/20;                           % sampling time: dt >> dt0
t_int = dt * (1/3);                  % integration time: dt0 <= t_int < dt
quant = 0.1;                         % LSB: vertical step

[t,s] = integration(t0,s0,dt0,dt,t_int,quant,0);

frame_length = 10;

range0 = (1 : (frame_length/dt0)) * dt0;
range = (1 : (frame_length/dt)) * dt;

%   - Divide timeline -
for k = 0 : (length(s0) / length(range0)) - 1
    
    t0_div (k+1,:) =  t0(  k*(length(range0)) + 1 : (k+1)*(length(range0)) ) ;
    s0_div(k+1,:)= s0(  k*(length(range0)) + 1 : (k+1)*(length(range0)) ) ;
    
end
%t0_ = t0_div(frame,:);    s0_ = s0_div(frame,:);   

for k = 0 : (length(s) / length(range)) - 1
    
    t_div (k+1,:) =  t(  k*(length(range)) + 1 : (k+1)*(length(range)) ) ;
    s_div(k+1,:)= s(  k*(length(range)) + 1 : (k+1)*(length(range)) ) ;
    
end
%t_ = t_div(frame,:);    s_ = s_div(frame,:);   