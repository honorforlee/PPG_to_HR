% Ivan Ny Hanitra - Master thesis
%       -- Linear regression of tx (3 methods) - output period, normalized Mean Squared Error, coefficient of determination - regression plot --

function [T,eps,R_sq,plot_reg] = periodicity(tx)
%   - with linear regression function -
% tbl = table([1:length(tx)]', tx','VariableNames',{'k','tx'});
% mdl = fitlm(tbl,'tx~k');
% T_est = mdl.Coefficients(2);    % period estimate
% RMSE_est = mdl.RMSE;            % Square root of the mean squared error, which estimates the standard deviation of the error distribution.

%   - with covariance matrix - 
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

%   - plot -
plot_reg=plot(j,tx,'r.', j,t_reg,'b-');
hold on 
title('Linear regression of t_{x,k}');
xlabel('k');
ylabel('t_{x,k},s');
legend('sampled t_{x,k}','linear regression of t_{x,k}','Location','northwest');

xlim=get(gca,'XLim');
ylim=get(gca,'YLim');

polyfit_str = ['\bf \tau = ' num2str(mean(tx) - T*mean(j)) ' + k*' num2str(T)];

equation = text(0.8*xlim(1)+0.15*xlim(2),0.3*ylim(1)+0.95*ylim(2),polyfit_str);
equation.Color = 'blue';
equation.FontSize = 14;
hold off

end