% Ivan NY HANITRA - Master thesis
%       -- Similarity functions  --

function eps = similarity(elm1, elm2, type)
if type == 'variance'
    eps = 0.25 * ( abs(elm1 - elm2) )^2;
    
elseif type == 'relative'
    eps = abs( (elm1-elm2)/elm1 );
    
else
    eps = 0;
end


