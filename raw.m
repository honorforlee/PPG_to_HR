tx_major = [1 2 3 3.2 4 5 6];
interval = delta_tx(tx_major);
T = median(interval);

for k = 1:length(interval)
    if interval(k) <= 0.5*T
        j = k+1;
        tx_sum = 0;
  
        while ( ~(0.8*T < tx_sum < 1.2*T)  && tx_sum < 2*T ) || tx_sum == 0
            tx_sum = sum(interval(k:j));
            j = j+1;
        end
    end
end