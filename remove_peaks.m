% Ivan NY HANITRA - Master thesis
%       -- Remorve peaks from major cluster  --   

function [kx_major,tx_major,sx_major,T] = remove_peaks(kx_major,tx_major,sx_major, T, kx,note_x)
    loop = 0;
    tx_neg = delta_tx(tx_major);
    loop_ = length(tx_neg);
    i=1;
    
    while loop < loop_
        for k = i:length(tx_neg)
            if tx_neg(k) < T - T/3
                if similarity(note_x(k), note_x(k+1), 'relative') > 0.2              % 20% note_x tolerance to suppress peak
                    
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
            end
            loop = loop +1;
        end
        
        loop = loop +1;
    end
