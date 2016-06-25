% Ivan NY HANITRA - Master thesis
%       -- Agglomerative clustering with Ward's criterion as linkage (decrease in variance for the cluster being merged) --

function [clust, clust_index] = agglo_clustering(note_x,kmax)
X = [note_x ]';
clust{1} = note_x;

for k = 2:kmax
    c = clusterdata(X,'linkage','ward','savememory','on','maxclust',k); 
    
    %     uv = unique(c);                                    % list of clutsters [1 .. k]
    %     n  = histc(c,uv);                                  % number of elements in each cluster (vector)
    
    for i = 1 : k     % inter clust
        clust_index{i,k} = find(c == i);
        clust{i,k} = note_x (clust_index{i,k});   % clust partition
        
        n = cellfun(@length,clust);
        
        num_F_(i) =( n(i,k) * (distance(mean(clust{i,k}), mean(note_x), 2))^2 ) / (k - 1);              % distance INTER - clust
        
        for j = 1 : n(i,k)     % intra clust
            den_F_d(j) = distance( clust{i,k}(j), mean(clust{i,k}), 2)^2 / (length(note_x) - k);        % distance INTRA - clust j
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
