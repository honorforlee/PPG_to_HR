tx = [1.286 2.23091 3.00448 4.26824 4.74249 6.08068 6.83541 7.75297 8.96053 10.168];

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

plot_reg=plot(j,tx,'r.', j,t_reg,'b-');
hold on 
title('Linear regression of t_{x,k}');
xlabel('k');
ylabel('t_{x,k},s');
legend('sampled t_{x,k}','linear regression of t_{x,k}');
hold off
