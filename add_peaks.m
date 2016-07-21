% Ivan NY HANITRA - Master thesis
%       -- Add peaks to major cluster --

function [kx_major, tx_major, sx_major, T] = add_peaks(kx,sx,tx, kx_major, kx_add,tx_pos)
if find(kx_add==0)                      % imaginary peak to create 
    zeros = find(kx_add==0);            % indexes where create a peak
            
    insert = @(a, x, n)cat(2,  x(1:n), a, x(n+1:end));      % insert(element inserted,array,position)
    
    for k = 1:length(zeros)
        kx_major = insert(kx_major(zeros(k)), kx_major, zeros(k) );
        kx_add(zeros(k)) = nan;
        kx_add = insert(nan, kx_add, zeros(k));                             % for concatenation with kx_major
        
        zeros = bsxfun(@plus , zeros, ones(1,length(zeros)));               % shift index of zeros when adding one element in tx_major, sx_major
    end
end 

    kx_major = horzcat(kx_major,kx_add);        % add peak to major cluster
    kx_major(isnan(kx_major)) = [];             % remove NaN values
    kx_major = sort(kx_major);                  % sort
    
    idxs = arrayfun(@(x)find(kx==x,1),kx_major);
          
    tx_pos(isnan(tx_pos)) = [];                 % remove nan values
    T_temp = mean(tx_pos);                      % peaks period (not considering missing peak)
    
    tx_major = tx(idxs);
    for k = 1:length(kx_major)-1                
        if kx_major(k)~=kx_major(k+1)           % peak to add
            tx_major(k+1) = tx(idxs(k+1));
        else                                    % peak to create
            tx_major(k+1) = tx(idxs(k)) + T_temp;
        end      
    end
   sx_major = sx(idxs);
    
T = mean(delta_tx(tx_major));