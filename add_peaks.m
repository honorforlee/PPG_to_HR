% Ivan NY HANITRA - Master thesis
%       -- Add peaks to major cluster --

function [kx_major, tx_major, sx_major, T] = add_peaksadd_peaks(t_,s_,td,d, kx_major,tx_major,sx_major, kx_add,tx_pos)
if find(kx_add==0)                      % imaginary peak to create & peak to add
    zeros = find(kx_add==0);            % indexes where create a peak
    tx_pos(isnan(tx_pos)) = [];         % remove nan values
    T_temp = mean(tx_pos);              % peaks period (not considering missing peak)
    
    insert = @(a, x, n)cat(2,  x(1:n), a, x(n+1:end));      % insert(element inserted,array,position)
    
    for k = 1:length(zeros)
        kx_major = insert(kx_major(zeros(k)), kx_major, zeros(k) );
        kx_add(zeros(k)) = nan;
        kx_add = insert(nan, kx_add, zeros(k));                             % for concatenation with kx_major
        
        zeros = bsxfun(@plus , zeros, ones(1,length(zeros)));               % shift index of zeros when adding one element in tx_major, sx_major
    end
    
    kx_major = horzcat(kx_major,kx_add);        % add peak to major cluster
    kx_major(isnan(kx_major)) = [];             % remove NaN values
    kx_major = sort(kx_major);                  % sort
    
    for k = 1:length(kx_major)-1                % peak to add
        if kx_major(k)~=kx_major(k+1)
            tx_major(k+1) = td(kx_major(k+1)) + (td(kx_major(k+1)+1)-td(kx_major(k+1))) .* d(kx_major(k+1))./(d(kx_major(k+1))-d(kx_major(k+1)+1));      % linear interpolation of dhi and dho to get tx (@zero crossing)
        else                                    % peak to create
            tx_major(k+1) = tx_major(k) + T_temp;
        end
        sx_major(k+1) = s_(kx_major(k+1)+1);        % local maxima
    end
else                                            % peak to add only
    kx_major = horzcat(kx_major,kx_add);        % add peak to major cluster
    kx_major(isnan(kx_major)) = [];             % remove NaN values
    kx_major = unique(kx_major);                % sort
    
    tx_major = td(kx_major) + (td(kx_major+1)-td(kx_major)) .* d(kx_major)./(d(kx_major)-d(kx_major+1));      % linear interpolation of dhi and dho to get tx (@zero crossing)
    sx_major = s_(kx_major+1);                  % local maxima
    
end

T = mean(delta_tx(tx_major));