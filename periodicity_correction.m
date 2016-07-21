%% PERIODICITY CORRECTION


    
    %       - Search for missing peaks -
    % Miss peaks inside the frame: add/create 2 peaks max in a hole
    loop = 0;
    
    while loop < 2
        tx_pos = delta_tx(tx_major);
        
        % Search for missing peaks inside the frame
        [kx_add,tx_pos] = missing_peaks(kx,tx, kx_major,tx_major, tx_pos,T, note_x,NOTE_major,eps);
        
        % Add/create peak to major cluster
        [kx_major, tx_major, sx_major, T] = add_peaks(kx,sx,tx, kx_major, kx_add,tx_pos);
       
        clearvars tx_pos kx_add;
        loop = loop+1;
    end
    
%     % Miss first peak
%     insert = @(a, x, n)cat(2,  x(1:n), a, x(n+1:end));      % insert(element inserted,array,position)
%     
%     if similarity(T,abs( tx(1) - tx_major(1) ), 'relative') < 0.2 || abs( tx(1) - tx_major(1) ) > T + T/5             % 20% relative error or more than 20% error above T
%         kx_major = insert(kx(1),kx_major,0);
%         tx_major = td(kx_major) + (td(kx_major+1)-td(kx_major)) .* d(kx_major)./(d(kx_major)-d(kx_major+1));          % linear interpolation of dhi and dho to get tx (@zero crossing)
%         sx_major = insert(sx_major(1),sx_major,0);                                                                    % for visibility
%         T = mean(delta_tx(tx_major));
%     end
%     
%     % Miss last peak
%     if similarity(T,abs( tx(end) - tx_major(end) ), 'relative') < 0.2  || abs( tx(end) - tx_major(end) ) > T + T/5    % 20% relative error or more than 20% error above T
%         kx_major = insert(kx(end),kx_major,length(kx_major));
%         tx_major = td(kx_major) + (td(kx_major+1)-td(kx_major)) .* d(kx_major)./(d(kx_major)-d(kx_major+1));          % linear interpolation of dhi and dho to get tx (@zero crossing)
%         sx_major = insert(sx_major(end),sx_major,length(sx_major));                                                   % for visibility
%         T = mean(delta_tx(tx_major));
%     end
    
    %   - Remove peaks from major cluster -
    [kx_major,tx_major,sx_major,T,warning] = remove_peaks(kx_major,tx_major,sx_major, T, kx, note_x);
    if warning == 1
        display('No peaks detected')
        return
    end