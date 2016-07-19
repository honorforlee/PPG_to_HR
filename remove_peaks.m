% Ivan NY HANITRA - Master thesis
%       -- Remorve peaks from major cluster  --

function [kx_major,tx_major,sx_major,T,warning] = remove_peaks(kx_major,tx_major,sx_major, T, kx,note_x)
if length(kx_major) > 2
    loop = 0;
    tx_neg = delta_tx(tx_major);
    loop_ = length(tx_neg);
    i=1;
    
    while loop < loop_ && length(kx_major) > 2
        for k = i:length(tx_neg)
            if tx_neg(k) < T - T/3 || tx_neg(k) < 1/3.17
                
                if note_x(kx==kx_major(k)) > note_x(kx==kx_major(k+1))      % compare which peak is more relevant
                    kx_major(k+1) = [];
                    tx_major(k+1) = [];
                    sx_major(k+1) = [];
                else
                    kx_major(k) = [];
                    tx_major(k) = [];
                    sx_major(k) = [];
                end
                tx_neg = delta_tx(tx_major);        % recompute tx_neg and T
                T = mean(delta_tx(tx_major));
                i=k;                                % start after peak removal
                break                               % exit from for loop
                
            end
            loop = loop +1;
        end
        
        loop = loop +1;
    end
    if length(kx_major) > 2
        warning = 0;
    else
        warning = 1;
        tx_major = nan;
        sx_major = nan;
        T = nan;
    end
else
    warning = 1;
    tx_major = nan;
    sx_major = nan;
    T = nan;
end