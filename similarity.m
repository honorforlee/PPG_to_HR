% Ivan NY HANITRA - Master thesis
%       -- Similarity functions  --

function eps = similarity(elm1, elm2, type, order,varargin)
if nargin == 3  % default
    if type == 'variance'
        eps = 0.25 * ( abs(elm1 - elm2) )^2;
        
    elseif type == 'relative'
        eps = abs( (elm1-elm2)/elm1 );
        
    else
        eps = 0;
    end   
    
elseif nargin == 4
    if type == 'variance'
        eps = 0.25 * ( abs(elm1 - elm2) )^2;
        
    elseif type == 'relative'
        if order == 1
            eps = abs( (elm1-elm2)/elm1 );
        elseif order == 2
            eps = abs( (elm1-elm2)/elm2 );
        end
    end
end

