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

[t_,s_] = integration(t0,s0,dt0,dt,t_int,quant,0);
 
frame_init = 5; frame_end = 10;

index = find(t_ >= frame_init & t_ <= frame_end);
t = t_(index);
s = s_(index);


%%
%   - Derivative -
d = s(2:end) -  s(1:end-1);
td = t(2:end);

kx = d > 0;
kx = find(kx(1:end-1) & ~kx(2:end));       % k_{x}:index where d > 0; d( k_{x} + 1 ) <= 0

sx = s(kx+1);                          % local maxima
tx = td(kx) + (td(kx+1)-td(kx)) .* d(kx)./(d(kx)-d(kx+1));      % linear interpolation of dhi and dho to get tx (@zero crossing)

kx_n = d < 0;                               % search for local minima
kx_n = find(kx_n(1:end-1) & ~kx_n(2:end));

    for k = 1:length(kx)                    % compute minima
        kx_index(k) = max( find( kx_n < kx(k) ) );
    end
sx_N = s(kx_n( kx_index ) + 1);
    
%   - Peaks notation -
note_1 = sx;
for k = 2:length(kx)-1
    note_1(k) = sx(k) - ( sx(k+1) + sx(k-1) )/2;                                % average peak value
end

note_3 = sx - sx_N;

note_x = 0.2*note_1 + 0.8*note_3;




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

