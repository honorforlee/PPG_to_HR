% Ivan NY HANITRA - Master thesis
%       -- Search for missing peaks inside major cluster --

function [kx_add,tx_pos] = missing_peaks(kx,tx, kx_major,tx_major, tx_pos,T, note_x,NOTE_major, eps)
kx_add = nan(1,length(kx_major));       % for horizontal concatenation

for k = 1:length(tx_pos)                % assume ONE missing/skipped peak
    if tx_pos(k) > T + T/3            % need enough large frame length to give weight to T
        left(k) = kx_major(k);
        right(k) = kx_major(k+1);
        kx_add_ = kx( kx(1,:) > left(k) & kx(1,:) < right(k));
        
        if length(kx_add_) == 1         % one peak present in the hole
            tx_pos_temp_left = tx(kx==kx_add_) - tx_major(k);
            tx_pos_temp_right = tx_major(k+1) - tx(kx==kx_add_);
            
            if similarity(NOTE_major, note_x(kx==kx_add_), 'variance') < 5*eps || note_x(kx==kx_add_) > NOTE_major || (similarity(T,tx_pos_temp_left,'relative') < 0.1 && similarity(T,tx_pos_temp_right,'relative') < 0.1)
                
                kx_add(k) = kx_add_;
                
            else
                tx_pos(k) = nan;            % to compute T not affected by missing tx_major
                kx_add(k) = 0;
            end
            clearvars kx_add_ tx_pos_temp;
            
        elseif length(kx_add_) >= 2     % more than one peak present in the hole
            for i = 1:length(kx_add_)
                kx_add_idx(i) = find(kx == kx_add_(i));
                kx_add_note(i) = note_x(kx_add_idx(i));
            end
            
            [value kx_add_max] = max(kx_add_note);
            tx_pos_temp_left = tx(kx_add_idx(kx_add_max)) - tx_major(k);         % delta_tx from major peak to added peak
            tx_pos_temp_right = tx_major(k+1) - tx(kx_add_idx(kx_add_max));
            
            if similarity(NOTE_major, value, 'variance') < 5*eps || value > NOTE_major || (similarity(T,tx_pos_temp_left,'relative') < 0.1 && similarity(T,tx_pos_temp_right,'relative') < 0.1)
                
                kx_add(k) = kx_add_(kx_add_max);        % max note_x index added to major cluster
                
            else
                tx_pos(k) = nan;            % to compute T not affected by missing tx_major
                kx_add(k) = 0;
            end
            clearvars kx_add_ kx_add_idx kx_add_note kx_add_max value tx_pos_temp ;
            
        else                            % no peak present in the hole => create peak
            tx_pos(k) = nan;            % to compute T not affected by missing tx_major
            kx_add(k) = 0;
            
        end
    end
    
end