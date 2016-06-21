% Ivan Ny Hanitra - Master thesis
%       -- Distance: manathan/taxicab(L1), euclidian(L2)  --

function [dist]=distance(X,Y,type)
if type == 1                        % norm 1 
    dist = sum( abs( X - Y ) );

elseif type == 2                    % norm 2
    dist = sqrt( sum(  (X - Y) .* (X - Y)    ));
end

