% Ivan NY HANITRA - Master thesis
%       -- Remove outliers according to cluster population -

% function X = outlier(note_x)
% std_score = abs ( normlist(note_x) );
%
% for k = 1:length(note_x)
%     if abs(note_x(k)) <= 2
%         X(k) = note_x(k);
%     else X(k) = nan;
%     end
% end

function kx = outlier(kx,clust_index)
size_ = size(clust_index);
clust_outlier = zeros(length(kx),size_(1)*size_(2));
k = 1;
for i = 1:size_(1)
    for j = 1:size_(2)
        
%        if length( clust_index{i,j} ) <= floor ( 0.05 * length(kx) ) && length( clust_index{i,j} ) ~ 0
         if length( clust_index{i,j} ) <= floor ( 0.1 * length(kx) ) && length( clust_index{i,j} ) ~ 0
            clust_outlier(1:length(clust_index{i,j}),k) =  clust_index{i,j};
            k = k+1;
        else
            k = k+1;
        end
        
    end
end

out_index = find(clust_outlier);
out = clust_outlier(out_index);
kx(out)=[];