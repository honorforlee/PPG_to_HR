% Ivan NY HANITRA - Master thesis
%       -- Plot signal --

Name = 'test1-1';           
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

plot(t,s,'ok');
