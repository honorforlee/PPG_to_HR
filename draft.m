kx_add = [nan nan 0 nan nan 0 nan nan];

kx_add_temp = kx_add;
kx_add_temp(isnan(kx_add_temp)) = [];   % remove NaN to find zeros

if find(kx_add==0)                   % one imaginary peak to create
    zeros = find(kx_add==0);             % indexes to create peak
    a =1;
else 
    a =0;
end

