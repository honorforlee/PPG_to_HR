load(strcat('3899985_0005m.mat'));
dt0 = 8e-3;
val(isnan(val)) = [];

t0 = (1:length(val)) * dt0;            % timeline
s0 = val(1,1:length(val));

s0  = (s0  - mean(s0 ))/sqrt(var(s0));        % rescale s on 0 (standard score of signal)

%   - Timeline, noise, integration, quantization -
dt = 1/10;                           % sampling time: dt >> dt0
t_int = dt * (1/3);                  % integration time: dt0 <= t_int < dt
quant = 0.1;                         % LSB: vertical step

[t,s] = integration(t0,s0,dt0,dt,t_int,quant,0);

[kx,tx,sx, dhi,dlo, td,d, kx_n,tx_N,sx_N, note_x] = signal_peaks(t,s);


periodicity(tx);


%%
plot_reg=plot(j,tx,'r.', j,t_reg,'b-');
hold on 
title('Linear regression of t_{x,k}');
xlabel('k');
ylabel('t_{x,k}, s');
legend('sampled t_{x,k}','linear regression of t_{x,k}','Location','northwest');

xlim=get(gca,'XLim');
ylim=get(gca,'YLim');

polyfit_str = ['\bf \tau = ' num2str(mean(tx) - T*mean(j)) ' + k*' num2str(T)];

equation = text(0.8*xlim(1)+0.15*xlim(2),0.3*ylim(1)+0.95*ylim(2),polyfit_str);
equation.Color = 'blue';
equation.FontSize = 14;
hold off

