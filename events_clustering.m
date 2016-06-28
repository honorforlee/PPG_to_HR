%   - Hierarchical clustering (agglomerative) -
function [tx,sx, dhi,dlo, td,d, tx_N,sx_N, note_x, clust_note_x, clust_tx, clust_periodicity, kmax, tx_major,sx_major ] = events_clustering(t,s)
%   - Derivative -
d = s(2:end) -  s(1:end-1);
%td = (  t(2:end) +  t(1:end-1) ) / 2;      % timeline of derivative shifted by t_sample/2
td = t(2:end);

kx = d > 0;
kx = find(kx(1:end-1) & ~kx(2:end));       % k_{x}:index where d > 0; d( k_{x} + 1 ) <= 0

%   - Local maxima sx, maximum slope around sx -
[tx,sx, dhi,dlo, tx_N,sx_N, note_x] = peaks_processing(t,s,kx);

% Initialization
kmax_init = 6;
[clust_index,  ~,~,  ~,~,  kmax, diff] = agglo_clustering(note_x, tx, kmax_init);

% Remove oultiers
kx = outlier(kx,clust_index, floor (0.05*length(kx)));      % remove cluster containing population <= 5% length(kx)

[tx,sx, dhi,dlo, tx_N,sx_N, note_x] = peaks_processing(t,s,kx);

if diff(2,2) >= 1    % EMPIRICAL: no clustering if 2-clustering clusters are too close
    
    % Initialization with outliers removed
    [clust_index,  clust_note_x,mean_clust,  clust_tx,clust_periodicity,  kmax, diff] = agglo_clustering(note_x, tx, kmax_init);
    
    % Search for best number of clusters
    div = 2;            % EMPIRICAL: merge clusters that are too close (    min(mean_clust difference) <= (max(mean_clust) - min(mean_clust)) / div
    while min( diff(2:end,kmax) ) <= ( max(mean_clust(:,kmax)) - min(mean_clust(:,kmax)) )/div && kmax >= 3
        
        kmax = kmax - 1;
        [clust_index,  clust_note_x,mean_clust,  clust_tx,clust_periodicity,  kmax, diff] = agglo_clustering(note_x, tx, kmax);
        
    end
    
    [~,clust_major_index] = max(mean_clust(:,kmax));
    kx_major = clust_index{clust_major_index,kmax};
    tx_major = tx(kx_major);
    sx_major = sx(kx_major);
    
else
    clust_note_x = nan;
    clust_tx = nan;
    [clust_periodicity(1),clust_periodicity(2),clust_periodicity(3)] = periodicity(tx);
    kmax = 1;
    tx_major = tx;
    sx_major = sx;
end
