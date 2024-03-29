% Ivan NY HANITRA - Master thesis
%       -- Agglomerative clustering with Ward's criterion as linkage (decrease in variance for the cluster being merged) --

function [clust_index,  clust_note_x,mean_clust,  clust_tx,clust_periodicity,  kmax_value, diff] = agglo_clustering(note_x, tx, kmax)
kmax_value = kmax;
X = [note_x ]';
clust_note_x{1} = note_x;

mean_clust = nan(kmax,kmax);
mean_clust(1,1) = mean(note_x);

for k = 2:kmax
    c = clusterdata(X,'linkage','ward','savememory','on','maxclust',k);
    
    %     uv = unique(c);                                    % list of clutsters [1 .. k]
    %     n  = histc(c,uv);                                  % number of elements in each cluster (vector)

    
    for i = 1 : k     % inter clust
        clust_index{i,k} = find(c == i);
        clust_note_x{i,k} = note_x (clust_index{i,k});   % cluster note_x partition      
        clust_tx{i,k} = tx (clust_index{i,k});           % cluster tx partition
        
        mean_clust(i,k) = mean(clust_note_x{i,k});       % cluster average note_x
        [clust_periodicity{i,k}(1),clust_periodicity{i,k}(2),clust_periodicity{i,k}(3)]=periodicity(clust_tx{i,k}); % cluster peiodicity output T,eps,R_sq
        
        n = cellfun(@length,clust_note_x);
        
        num_F_(i) =( n(i,k) * (distance(mean(clust_note_x{i,k}), mean(note_x), 2))^2 ) / (k - 1);              % distance INTER - clust
        
        for j = 1 : n(i,k)     % intra clust
            den_F_d(j) = distance( clust_note_x{i,k}(j), mean(clust_note_x{i,k}), 2)^2 / (length(note_x) - k);        % distance INTRA - clust j
        end
        
        den_F_(i) = sum(den_F_d);
        clearvars den_F_d;
    end
    
    %         figure(k);
    %         plot(clust{i,k} , '.');
    %         hold on
    %     end
    %         hold off
    
    num_F(k) = sum(num_F_);
    den_F(k) = sum(den_F_);
    F(k) = num_F(k) / den_F(k);             % F-statistics notation
    warning('off','all');
end    

%  Search best number of clusters
mean_clust_ = sort(mean_clust);      

for k = 2 : kmax
    for i = 2:kmax
    diff(i,k) = mean_clust_(i,k)-mean_clust_(i-1,k);
    end
end