% Ivan NY HANITRA - Master thesis
% -- Compute distance between tx --

function delta = delta_tx(tx,step,varargin)
if nargin == 1  % default value of step = 1
    for k = 1:length(tx)-1
        delta(k) = tx(k+1) - tx(k);
    end
else
    if step == 2
    
    for k = 1:length(tx)-2
        delta(k) = tx(k+2) - tx(k);
    end
    end
end

    