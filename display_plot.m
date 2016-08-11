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


j = [1:length(tx)];
A_matrix = cov(j,tx,1);                                          % output is normalized by the number of observations length(tx)
T_matrix = A_matrix(3)/A_matrix(1);                              % peaks eriodicity
T_0_matrix = mean(tx) - T_matrix * mean(j); 

R_2_matrix  = ( A_matrix(3) / sqrt(A_matrix(1)*A_matrix(4)))^2;  % coefficient of determination

%   - with loops -
for k = 1:length(tx)
   A(k) = k * tx(k);
end   
sum_A = cumsum(A);

cov_xy = (sum_A(length(tx)) / length(tx)) - mean([1:length(tx)])*mean(tx);
cov_xx = (length(tx)^2 -1)/12;

T = cov_xy/cov_xx;                                              % peaks eriodicity
t_reg = mean(tx) - T*mean(j) + j*T;

for k = 1:length(tx)
   B(k) = ( tx(k) - mean(tx) + ( mean([1:length(tx)]) - k )*T )^2;
   C(k) = (mean(tx)-T*mean(j) + k*T  - mean(tx))^2;
   D(k) = (tx(k) - mean(tx))^2;
end
sum_B = cumsum(B);                                              % Mean Squared Error
sum_C = cumsum(C); sum_D = cumsum(D);                                             
    
eps = sum_B(length(tx)) / (length(tx)*T^2);                     % MSE normalized with respect to #samples and period
R_sq = sum_C(length(tx)) / sum_D(length(tx));                   % coefficient of determination R^2

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

