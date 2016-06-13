function [T,eps] = periodicity(tx)
% k = [1:length(tx)];
% A = cov(k,tx,1);                    % output is normalized by the number of observations length(tx)
% T = A(3)/A(1);                      % peaks eriodicity
% R_2  = ( A(3) / sqrt(A(1)*A(4)))^2;  % coefficient of determination

N = length(tx);

for k = 1:N
   A(k) = k * tx(k);
end   
sum_A = cumsum(A);

cov_xy = (sum_A(N) / N) - mean([1:N])*mean(tx);
cov_xx = (N^2 -1)/12;

T = cov_xy/cov_xx;                  % peaks eriodicity

for k = 1:N
   B(k) = ( tx(k) - mean(tx) + ( mean([1:N]) - k )*T )^2;
end
sum_B = cumsum(B);                  % Mean Squared Error

eps = sum_B(N) / (N*T^2);           % MSE normalized with respect to # samples and period

end