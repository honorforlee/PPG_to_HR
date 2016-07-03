function [kx_major] = min_variance(kx,tx,sx,note_x, eps)
%   - Clustering according to minimum variance of note_x -
kx_ = kx;

for i = 1:length(kx)
    if kx_(i) ~= 0
        
        idx(1,i) = kx(i);
        note(1,i) = note_x(i);
        per(1,i) = tx(i);
        j = 2;
        
        for k = i + 1 : length(kx)
            
            if  var([note_x(i) note_x(k)],1) < eps;
                
                note(j,i) = note_x(k);
                per(j,i) = tx(k);
                idx(j,i) = kx(k);
                
                kx_(k) = 0;
                
                j = j+1;
            else
                note(j,i)=nan;
                per(j,i)=nan;
                idx(j,i)=nan;
                
                j = j+1;
            end
            
        end
    end
end

%   - Remove columns -
zero = find(~idx(1,:));

for k = 1:length(zero)
    note(:,zero(k))=[];
    per(:,zero(k))=[];
    idx(:,zero(k))=[];
    zero = bsxfun(@minus ,zero,ones(1,length(zero))) ;
end

%   - Create cells: kx-tx-note_x -
L = size(idx);
for k = 1:L(2)
    NAN_ = ~isnan(per(:,k));         % extract non NAN value of per, note
    NAN_idx = idx(NAN_,k);
    NAN_per = per(NAN_,k);
    NAN_note = note(NAN_,k);
    
    idx_ = NAN_idx(find(NAN_idx));     % extract non zero value of per, note
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
        
        [PER_T(k),PER_eps(k), PER_R(k)] = periodicity(clust_cell{k,2});    % note
        
        if PER_eps(k) <= 0.01
            PER_eps(k) = 0.01;
        end
        
        NOTE(k) = mean(clust_cell{k,3});
        SIZE(k) = length(clust_cell{k,1});
        
        clust_note(k) = (0.2 * NOTE(k) + 0.4 * SIZE(k)) / (PER_eps(k)/0.4);
        
    else
        
        PER_T(k) = 0; PER_eps(k) = 0; PER_R(k) = 0;
        NOTE(k) = mean(clust_cell{k,2});
        SIZE(k) = length(clust_cell{k,1});
        clust_note(k) = 0;
        
    end
end

tbl_note = table([1:L(2)]', PER_T',PER_eps', PER_R', NOTE', SIZE', clust_note','VariableNames',{'Cluster','T','eps','R','Note_x','Size','Cluster_note'})

%   - Major cluster -
[note_major major_idx] = max(clust_note);   
kx_major = clust_cell{major_idx,1};

Nrows = max(cellfun(@numel,clust_cell));
X = nan(Nrows(1),L(2));
for iCol = 1:L(2)
  X(1:numel(clust_cell{iCol}),iCol) = clust_cell{iCol};
end

%   - Merge sub-major clusters -
for k = 1:L(2)
    
    if var([note_major clust_note(k)],1) < 5*eps && k ~= major_idx

        kx_major = vertcat(kx_major,clust_cell{k,1});

    end
        
end