% Ivan NY HANITRA - Master thesis
        % -- Compute distance matrix between tx --

function delta = delta_tx(tx)
for k = 1:length(tx)-1
    delta(k) = tx(k+1) - tx(k);
    
end