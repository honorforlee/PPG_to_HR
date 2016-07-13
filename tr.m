%   - Select major cluster + merge sub-major clusters -
Nrows = max(cellfun(@numel,clust_cell));
X = nan(Nrows(1),L(2));
for iCol = 1:L(2)
    X(1:numel(clust_cell{iCol}),iCol) = clust_cell{iCol};       % copy idx values of each cluster into X
end
clust_merge = nan(Nrows(1),L(2));

if L(2) >= 2                                             % more than 1 cluster
    if all(clust_note == 0)
        [NOTE_major major_idx] = max(NOTE);
        kx_major = clust_cell{major_idx,1}';
        
        NOTE_comp = NOTE_major;                    % EMP
        
        if all(NOTE <= NOTE_comp/2)                % EMP
            for k = 1:L(2)
                if var([NOTE_major NOTE(k)],1) < 7*eps                                            % EMPIRICAL: compare NOTE to NOTE of major cluster
                    NOTE_major = NOTE(k) * ( Nrows(1) - sum(isnan(X(:,k))) ) + NOTE_major * ( L(2)*Nrows(1) - sum(sum(isnan(clust_merge))) );
                    clust_merge(:,k) = X(:,k);                                                    % merge cluster to major cluster
                    NOTE_major = NOTE_major / ( L(2)*Nrows(1) - sum(sum(isnan(clust_merge))) );   % recompute NOTE_major
                   
                    NOTE_comp = NOTE_major;        % EMP
                end
            end
        else
            for k = 1:L(2)
                if var([NOTE_major NOTE(k)],1) < 7*eps && NOTE(k) > NOTE_comp                                 % EMPIRICAL: compare NOTE to NOTE of major cluster
                    NOTE_major = NOTE(k) * ( Nrows(1) - sum(isnan(X(:,k))) ) + NOTE_major * ( L(2)*Nrows(1) - sum(sum(isnan(clust_merge))) );
                    clust_merge(:,k) = X(:,k);                                                      % merge cluster to major cluster
                    NOTE_major = NOTE_major / ( L(2)*Nrows(1) - sum(sum(isnan(clust_merge))) );     % recompute NOTE_major
                   
                    NOTE_comp = NOTE_major;       % EMP
                end
            end
        end
    else
        if all(NOTE <= max(NOTE)/3)               % EMP
            [NOTE_major major_idx] = max(NOTE);
            kx_major = clust_cell{major_idx,1}';
            
            NOTE_comp = NOTE_major;               % EMP
            
            for k = 1:L(2)
                if var([NOTE_major NOTE(k)],1) < 7*eps                                               % EMPIRICAL: compare NOTE to NOTE of major cluster
                    NOTE_major = NOTE(k) * ( Nrows(1) - sum(isnan(X(:,k))) ) + NOTE_major * ( L(2)*Nrows(1) - sum(sum(isnan(clust_merge))) );
                    clust_merge(:,k) = X(:,k);                                                       % merge cluster to major cluster
                    NOTE_major = NOTE_major / ( L(2)*Nrows(1) - sum(sum(isnan(clust_merge))) );      % recompute NOTE_major
                    
                    NOTE_comp = NOTE_major;       % EMP
                end
            end
            
        else
            [clust_note_temp idx_temp] = max(clust_note);
            
            if NOTE(idx_temp)  > max(NOTE)/3      % EMP
                major_idx = idx_temp;
                kx_major = clust_cell{major_idx,1}';
                clust_note_major = clust_note_temp;
                NOTE_major = NOTE(major_idx);
                
                NOTE_comp = NOTE_major;           % EMP
                
                for k = 1:L(2)
                    if clust_note(k) > NOTE_comp && var([NOTE_major NOTE(k)],1) < 7*eps && NOTE(k) > NOTE_comp      % EMPIRICAL: compare NOTE to NOTE of major cluster
                        NOTE_major = NOTE(k) * ( Nrows(1) - sum(isnan(X(:,k))) ) + NOTE_major * ( L(2)*Nrows(1) - sum(sum(isnan(clust_merge))) );
                        clust_merge(:,k) = X(:,k);                                                  % merge cluster to major cluster
                        NOTE_major = NOTE_major / ( L(2)*Nrows(1) - sum(sum(isnan(clust_merge))) ); % recompute NOTE_major
                        
                        NOTE_comp = NOTE_major;           % EMP
                    elseif clust_note(k) == 0 && var([NOTE_major NOTE(k)],1) < 7*eps && NOTE(k) > NOTE_comp
                        NOTE_major = NOTE(k) * ( Nrows(1) - sum(isnan(X(:,k))) ) + NOTE_major * ( L(2)*Nrows(1) - sum(sum(isnan(clust_merge))) ); % EMPIRICAL: compare NOTE to NOTE of major cluster (that have max clust_note)
                        clust_merge(:,k) = X(:,k);                                                  % merge cluster to major cluster
                        NOTE_major = NOTE_major / ( L(2)*Nrows(1) - sum(sum(isnan(clust_merge))) ); % recompute NOTE_major
                        
                        NOTE_comp = NOTE_major;           % EMP
                    end
                end
                
            else
                [NOTE_major major_idx] = max(NOTE);
                kx_major = clust_cell{major_idx,1}';
                
                NOTE_comp = NOTE_major;           % EMP
                
                for k = 1:L(2)
                    if var([NOTE_major NOTE(k)],1) < 7*eps  && NOTE(k) > NOTE_comp                              % EMPIRICAL: compare NOTE to NOTE of major cluster
                        NOTE_major = NOTE(k) * ( Nrows(1) - sum(isnan(X(:,k))) ) + NOTE_major * ( L(2)*Nrows(1) - sum(sum(isnan(clust_merge))) );
                        clust_merge(:,k) = X(:,k);                                                      % merge cluster to major cluster
                        NOTE_major = NOTE_major / ( L(2)*Nrows(1) - sum(sum(isnan(clust_merge))) );     % recompute NOTE_major
                        
                        NOTE_comp = NOTE_major;           % EMP
                    end
                end
            end
        end
    end
    
    clust_merge(isnan(clust_merge)) = [];         % remove NaN values
    clust_merge = unique(clust_merge);            % remove repeated elements ans sort array
    kx_major(1,1:length(clust_merge)) = clust_merge;
    
else                                                    % one cluster only
    major_idx = 1;
    kx_major = clust_cell{major_idx,1}';
    clust_note_major = clust_note;
    NOTE_major = NOTE;
end

if length(kx_major) >= 2