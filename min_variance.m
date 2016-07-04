function [kx_major,tx_major,sx_major, T] = min_variance(t_,s_, td,d, kx,tx,sx,note_x, eps)
%   - Clustering according to minimum variance of note_x -
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
    NAN_ = ~isnan(per(:,k));         % extract non NAN value of idx, per, note
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
        
        if PER_eps(k) <= 0.1
            PER_eps(k) = 0.1;
        end
                
        NOTE(k) = mean(clust_cell{k,3});       
        clust_note(k) = (0.4 * NOTE(k) + 0.3 * SIZE(k)) / (PER_eps(k)/0.3);
        
    else
        SIZE(k) = length(clust_cell{k,1});
        PER_T(k) = 0; PER_eps(k) = 0; PER_R(k) = 0;
        NOTE(k) = mean(clust_cell{k,3});
        clust_note(k) = 0;
        
    end
end

tbl_note = table([1:L(2)]', SIZE', PER_T',PER_eps', PER_R', NOTE', clust_note','VariableNames',{'Cluster','Size','T','eps','R','Note_x','Cluster_note'})

%   - Major cluster -
for k = 1:L(2)
    if NOTE(k) > 0
        clust_note_pos(k) = clust_note(k);
    end
end
[clust_note_max major_idx] = max(clust_note_pos);
kx_major = clust_cell{major_idx,1};

%   - Merge sub-major clusters -
Nrows = max(cellfun(@numel,clust_cell));
X = nan(Nrows(1),L(2));
for iCol = 1:L(2)
    X(1:numel(clust_cell{iCol}),iCol) = clust_cell{iCol};
end
clust_merge = nan(Nrows(1),L(2));
clust_merge(:,major_idx) = X(:,major_idx);

for k = 1:L(2)
    if NOTE(k) > 0
        if var([max(NOTE) NOTE(k)],1) < eps && k ~= major_idx
            clust_merge(:,k) = X(:,k);
        end
    end
end

clust_merge(isnan(clust_merge)) = [];   % remove NaN values

%   - Major peaks -
kx_major = unique(clust_merge);      % remove repeated elements
tx_major = td(kx_major) + (td(kx_major+1)-td(kx_major)) .* d(kx_major)./(d(kx_major)-d(kx_major+1));      % linear interpolation of dhi and dho to get tx (@zero crossing)
sx_major = s_(kx_major+1);          % local maxima
T = mean(delta_tx(tx_major));

%   - Modify cluster according to periodicity -
% ADD PEAK
tx_pos = delta_tx(tx_major);
kx_add = nan(1,length(kx_major));

for k = 1:length(tx_pos)                % assume ONE missing/skipped peak
    if tx_pos(k) > T + 0.5*T            % add a major peak
        left(k) = kx_major(k);
        right(k) = kx_major(k+1);
        kx_add_ = kx( kx(1,:) > left(k) & kx(1,:) < right(k));
               
        if length(kx_add_) == 1
            kx_add(k) = kx_add_;
        elseif length(kx_add_) >= 2
            for i = 1:length(kx_add_)
                kx_add_idx(i) = find(kx == kx_add_(i));
                kx_add_note(i) = note_x(kx_add_idx(i));
            end
            
            [value kx_add_max] = max(kx_add_note);
            kx_add(k) = kx_add_(kx_add_max);
            
            clearvars kx_add_idx kx_add_note kx_add_max value ;
        end
   
    end
    clearvars kx_add_;
end

kx_major = horzcat(kx_major,kx_add);        % add peak to major cluster
kx_major(isnan(kx_major)) = [];             % remove NaN values
kx_major = unique(kx_major);                % sort

tx_major = td(kx_major) + (td(kx_major+1)-td(kx_major)) .* d(kx_major)./(d(kx_major)-d(kx_major+1));      % linear interpolation of dhi and dho to get tx (@zero crossing)
sx_major = s_(kx_major+1);          % local maxima
T = mean(delta_tx(tx_major));

% REMOVE PEAK
tx_neg = delta_tx(tx_major);
for k = 1:length(tx_neg)
    if tx_neg(k) < T - T*0.5        % remove peaks - another loop because of matrix size issue 
        kx_major(k+1) = nan;
    end
end

kx_major(isnan(kx_major)) = [];             % remove NaN values
kx_major = unique(kx_major);                % sort

tx_major = td(kx_major) + (td(kx_major+1)-td(kx_major)) .* d(kx_major)./(d(kx_major)-d(kx_major+1));      % linear interpolation of dhi and dho to get tx (@zero crossing)
sx_major = s_(kx_major+1);          % local maxima
T = mean(delta_tx(tx_major));