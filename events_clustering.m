% Ivan Ny Hanitra - Master thesis
%       -- Algorithm for clustering the events through in an agglomerative way --

function [clust_note_x, clust_tx, clust_periodicity, kmax, tx_major,sx_major ] = events_clustering(kx,tx,sx,note_x)
% Initialization
kmax_init = 6;
[clust_index,  ~,~,  ~,~,  kmax, diff] = agglo_clustering(note_x, tx, kmax_init);

% Remove oultiers
[kx,tx,sx,note_x] = outlier(kx,tx,sx,note_x, clust_index, floor (0.05*length(kx)));      % remove cluster containing population <= 5% length(kx)

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