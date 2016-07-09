% Ivan NY HANITRA - Master thesis
%       -- Clustering according to minimum variance of note_x  --

function [kx_major,tx_major,sx_major, T] = min_variance(t_,s_, td,d, kx,tx,sx,note_x, eps)
kx_ = kx;
for i = 1:length(kx)
    if kx_(i) ~= 0
        
        idx(1,i) = kx(i);
        per(1,i) = tx(i);
        note(1,i) = note_x(i);
        j = 2;
        
        for k = i + 1 : length(kx)
            
            if  var([note_x(i) note_x(k)],1) < eps;
                
                idx(j,i) = kx(k);
                per(j,i) = tx(k);
                note(j,i) = note_x(k);
                
                kx_(k) = 0;
                j = j+1;
            else
                idx(j,i)=nan;
                per(j,i)=nan;
                note(j,i)=nan;
                
                j = j+1;
            end
            
        end
    end
end

%   - Remove columns -
zero = find(~idx(1,:));

for k = 1:length(zero)
    idx(:,zero(k))=[];
    per(:,zero(k))=[];
    note(:,zero(k))=[];
    zero = bsxfun(@minus ,zero,ones(1,length(zero))) ;
end

%   - Create cells: kx-tx-note_x -
L = size(idx);
for k = 1:L(2)
    NAN_ = ~isnan(idx(:,k));         % extract non NAN value of idx, per, note
    NAN_idx = idx(NAN_,k);
    NAN_per = per(NAN_,k);
    NAN_note = note(NAN_,k);
    
    idx_ = NAN_idx(find(NAN_idx));     % extract non zero value of idx, per, note
    per_ = NAN_per(find(NAN_per));
    note_ = NAN_note(find(NAN_note));
    
    clust_cell{k,1} = idx_;            % kx
    clust_cell{k,2} = per_;            % tx
    clust_cell{k,3} = note_;           % note_x
    
    clear NAN_ NAN_idx NAN_per NAN_note idx_ per_ note_
end

%   - Cluster notation: size, note_x, periodicity -
for k = 1:L(2)
    if length(clust_cell{k,1}) > 2
        
        SIZE(k) = length(clust_cell{k,1});
        [PER_T(k),PER_eps(k), PER_R(k)] = periodicity(clust_cell{k,2});    % note
        
        if PER_eps(k) <= 0.1            % best periodicity note set to 0.1 otherwise increase too much the cluster note 
            PER_eps(k) = 0.1;
        end
        
        NOTE(k) = mean(clust_cell{k,3});
        clust_note(k) = (0.7 * NOTE(k) + 0.2 * SIZE(k)) / (PER_eps(k)/0.1);             % cluster note EMPIRICAL
        
    else
        SIZE(k) = length(clust_cell{k,1});
        PER_T(k) = 0; PER_eps(k) = 0; PER_R(k) = 0;
        NOTE(k) = mean(clust_cell{k,3});
        clust_note(k) = 0;
        
    end
end

tbl_note = table([1:L(2)]', SIZE', PER_T',PER_eps', PER_R', NOTE', clust_note','VariableNames',{'Cluster','Size','T','eps','R','Note_x','Cluster_note'})

%   - Major cluster -
if L(2) >= 2            % more than 1 cluster
for k = 1:L(2)
    if NOTE(k) > 0
        clust_note_pos(k) = clust_note(k);
    end
end
if all(clust_note_pos == 0)      % all cluster have zero note (all size <= 2)
[NOTE_major major_idx] = max(NOTE);
kx_major = clust_cell{major_idx,1}'; 
else 
[clust_note_major major_idx] = max(clust_note_pos);
kx_major = clust_cell{major_idx,1}';
NOTE_major = NOTE(major_idx);
end

else
    major_idx=1;
    clust_note_pos = clust_note;
    clust_note_major = clust_note;
    NOTE_major = NOTE;
    kx_major = clust_cell{1,1}';
end



tx_major = td(kx_major) + (td(kx_major+1)-td(kx_major)) .* d(kx_major)./(d(kx_major)-d(kx_major+1));      % linear interpolation of dhi and dho to get tx (@zero crossing) DEBUG
sx_major = s_(kx_major+1);                                                                                % DEBUG

%   - Merge sub-major clusters -
Nrows = max(cellfun(@numel,clust_cell));
X = nan(Nrows(1),L(2));
for iCol = 1:L(2)   
    X(1:numel(clust_cell{iCol}),iCol) = clust_cell{iCol};       % copy idx values of each cluster into X
end
clust_merge = nan(Nrows(1),L(2));
clust_merge(:,major_idx) = X(:,major_idx);

for k = 1:L(2)
    if NOTE(k) > 0 && ~all(clust_note_pos == 0)
        if var([clust_note_major clust_note(k)],1) < 7*eps && k ~= major_idx       % EMPIRICAL: compare cluter_note to max(cluster_note)
            if var([NOTE_major NOTE(k)],1) < 7*eps && k ~= major_idx               % EMPIRICAL: compare NOTE to NOTE of major cluster
                NOTE_major = NOTE(k) * ( Nrows(1) - sum(isnan(X(:,k))) ) + NOTE_major * ( L(2)*Nrows(1) - sum(sum(isnan(clust_merge))) );
                clust_merge(:,k) = X(:,k);
                NOTE_major = NOTE_major / ( L(2)*Nrows(1) - sum(sum(isnan(clust_merge))) ); % new NOTE_major after merging
            end
        else
            if var([NOTE_major NOTE(k)],1) < 7*eps && k ~= major_idx               % EMPIRICAL: compare NOTE to NOTE of major cluster - case cluster of 1/2 elements containing major peaks
                NOTE_major = NOTE(k) * ( Nrows(1) - sum(isnan(X(:,k))) ) + NOTE_major * ( L(2)*Nrows(1) - sum(sum(isnan(clust_merge))) );
                clust_merge(:,k) = X(:,k);
                NOTE_major = NOTE_major / ( L(2)*Nrows(1) - sum(sum(isnan(clust_merge))) );
            end
            
        end
    elseif  NOTE(k) > 0 && all(clust_note == 0)
        if var([NOTE_major NOTE(k)],1) < 7*eps && k ~= major_idx               % EMPIRICAL: compare NOTE to NOTE of major cluster - case cluster of 1/2 elements containing major peaks
            NOTE_major = NOTE(k) * ( Nrows(1) - sum(isnan(X(:,k))) ) + NOTE_major * ( L(2)*Nrows(1) - sum(sum(isnan(clust_merge))) );
            clust_merge(:,k) = X(:,k);
            NOTE_major = NOTE_major / ( L(2)*Nrows(1) - sum(sum(isnan(clust_merge))) );
        end
        
    end
end

clust_merge(isnan(clust_merge)) = [];         % remove NaN values
clust_merge = unique(clust_merge);            % remove repeated elements ans sort array

%   - Major peaks -
kx_major(1,1:length(clust_merge)) = clust_merge; 

tx_major = td(kx_major) + (td(kx_major+1)-td(kx_major)) .* d(kx_major)./(d(kx_major)-d(kx_major+1));      % linear interpolation of dhi and dho to get tx (@zero crossing)
sx_major = s_(kx_major+1);          % local maxima
T = mean(delta_tx(tx_major));

%   - Modify cluster according to periodicity -
% ADD PEAK TO MAJOR CLUSTER
loop = 0;

while loop < 2
tx_pos = delta_tx(tx_major);
kx_add = nan(1,length(kx_major));       % for horizontal concatenation

    for k = 1:length(tx_pos)                % assume ONE missing/skipped peak
    if tx_pos(k) > T + 0.5*T            % need enough large frame length to give weight to T
        left(k) = kx_major(k);
        right(k) = kx_major(k+1);
        kx_add_ = kx( kx(1,:) > left(k) & kx(1,:) < right(k));
        
        if length(kx_add_) == 1         % one peak present in the hole
            kx_add(k) = kx_add_;
            
        elseif length(kx_add_) >= 2     % more than one peak present in the hole
            for i = 1:length(kx_add_)
                kx_add_idx(i) = find(kx == kx_add_(i));
                kx_add_note(i) = note_x(kx_add_idx(i));
            end
            
            [value kx_add_max] = max(kx_add_note);
            kx_add(k) = kx_add_(kx_add_max);        % max note_x index added to major cluster
            
            clearvars kx_add_idx kx_add_note kx_add_max value ;
            
        else                            % no peak present in the hole => create peak
            kx_add_ = nan;
            tx_pos(k) = nan;            % to compute T not affected by missing tx_major
            kx_add(k) = 0;
            
        end
        clearvars kx_add_;
    end
    
end

[kx_major, tx_major, sx_major, T] = add_peaks(t_,s_,td,d, tx_pos,kx_major,tx_major,kx_add);
loop = loop+1;
end

% REMOVE PEAK FROM MAJOR CLUSTER
tx_neg = delta_tx(tx_major);
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