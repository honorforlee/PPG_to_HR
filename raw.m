kx = 100;
tx = ones(1,kx);
t_sum =0;
t_ = tx;

t_(1)= 10;
for k = 2:length(kx)
    t_(k) = k * 7*tx(k);
end
t_sum = t_(length(kx))/ length(kx);

T = (12*t_sum - 6* (length(kx)+1) * mean(tx)) / (length(kx)*length(kx) - 1);    % linear regression of peaks period (T-periodic)