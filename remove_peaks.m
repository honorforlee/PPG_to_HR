% Ivan NY HANITRA - Master thesis
%       -- remove peaks to major cluster --

function [kx_major, tx_major, sx_major, T] = remove_peaks(t_,s_,td, tx_neg)
for k = 1:length(tx_neg)
    if tx_neg(k) < T - T*0.5        % remove peaks - another loop because of matrix size inconsistency
        kx_major(k+1) = nan;
    end
end

kx_major(isnan(kx_major)) = [];             % remove NaN values
kx_major = unique(kx_major);                % sort

tx_major = td(kx_major) + (td(kx_major+1)-td(kx_major)) .* d(kx_major)./(d(kx_major)-d(kx_major+1));      % linear interpolation of dhi and dho to get tx (@zero crossing)
sx_major = s_(kx_major+1);          % local maxima
T = mean(delta_tx(tx_major));