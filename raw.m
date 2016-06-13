tx = [1.286 2.23091 3.00448 4.26824 4.74249 6.08068 6.83541 7.75297 8.96053 10.168];

N = length(tx);

for k = 1:N
   A(k) = k * tx(k);
end   
sum_A = cumsum(A);

cov_xy = (sum_A(N) / N) - mean([1:N])*mean(tx);
cov_xx = (N^2 -1)/12;

T = cov_xy/cov_xx;

for k = 1:N
   B(k) = ( tx(k) - mean(tx) + ( mean([1:N]) - k )*T )^2;
end
sum_B = cumsum(B);

eps = sum_B(N) / (N*T^2);