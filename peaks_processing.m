% Ivan NY HANITRA - Master thesis
%       -- Local maxima sx, maximum slope around sx --

function [tx,sx, dhi,dlo, kx_n,tx_N,sx_N, note_x] = peaks_processing(t,s,kx)

d = s(2:end) -  s(1:end-1);
%td = (  t(2:end) +  t(1:end-1) ) / 2;    % shift derivative timeline of t_spl/2
td = t(2:end);

sx = s(kx+1);                          % local maxima
tx = td(kx) + (td(kx+1)-td(kx)) .* d(kx)./(d(kx)-d(kx+1));      % linear interpolation of dhi and dho to get tx (@zero crossing)

dhi = d(kx);
dlo = d(kx+1);

for k = 1:length(kx)
    i = kx(k)-1;   while i > 0             && d(i) >= dhi(k); dhi(k) = d(i); i = i-1; end    % search for maximum positive slope at kx-
    i = kx(k)+2;   while i < length(d) && d(i) <= dlo(k); dlo(k) = d(i); i = i+1; end    % search for maximmum negative slope at kx+
end

kx_n = d < 0;                               % search for local minima
kx_n = find(kx_n(1:end-1) & ~kx_n(2:end));

if kx_n(1) < kx(1)
    for k = 1:length(kx)
        kx_index(k) = max( find( kx_n < kx(k) ) );
    end
    sx_N = s(kx_n( kx_index ) + 1);
    tx_N = td(kx_n( kx_index )) + (td(kx_n( kx_index )+1)-td(kx_n( kx_index ))) .* d(kx_n( kx_index ))./(d(kx_n( kx_index ))-d(kx_n( kx_index )+1));
else
    j = 1;
    while kx_n(1) > kx(j)                   % find first minima preceding maxima
    kx_index(j) = nan;
    sx_N(j) = nan;
    tx_N(j) = nan;
    j = j+1;
    end
        
    for k = j:length(kx)
        kx_index(k) = max( find( kx_n < kx(k) ) );
        sx_N(k) = s(kx_n( kx_index(k) ) + 1);
        tx_N(k) = td(kx_n( kx_index(k) )) + (td(kx_n( kx_index(k) )+1)-td(kx_n( kx_index(k) ))) .* d(kx_n( kx_index(k) ))./(d(kx_n( kx_index(k) ))-d(kx_n( kx_index(k) )+1));
    end
  clearvars j;  
end

% %   - Peaks notation
% sx_norm = sx * (2/max(sx));                   % normalize amplitude
% note_1 = sx_norm;                     
% 
% note_2 = dhi - dlo;                           % maximum slope difference around peak
% 
% for k = 1:length(tx)                          
%     if tx(k) >= tx_N(k)
%         delta(k) = sx_norm(k) - sx_N(k);
%     else                                     % if minimum out of frame, take first min in the frame 
%         j = k;  
%         while isnan(sx_N(j))                
%             j = j+1;
%         end
%         delta(k) = sx_norm(k) - sx_N(j);
%         clearvars j; 
%     end
% end
% note_3 = delta;
% 
% for k = 2:length(kx)-1
%     note_1(k) = sx_norm(k) - ( sx_norm(k+1) - sx_norm(k-1) )/2;            % average peak value
%     note_3(k) = delta(k) - ( delta(k+1) - delta(k-1) )/2;                  % average peak to peak value
% end
% 
% note_x = 0.1*note_1 + 0.1*note_2 + 0.8*delta;

%   - Peaks notation
note_1 = sx;                     

note_2 = dhi - dlo;                           % maximum slope difference around peak

for k = 1:length(tx)                          
    if tx(k) >= tx_N(k)
        delta(k) = sx(k) - sx_N(k);
    else                                     % if minimum out of frame, take first min in the frame 
        j = k;  
        while isnan(sx_N(j))                
            j = j+1;
        end
        delta(k) = sx(k) - sx_N(j);
        clearvars j; 
    end
end
note_3 = delta;

for k = 2:length(kx)-1
    note_1(k) = sx(k) - ( sx(k+1) - sx(k-1) )/2;                           % average peak value
    note_3(k) = delta(k) - ( delta(k+1) - delta(k-1) )/2;                  % average peak to peak value
end

note_x = 0.1*note_1 + 0.1*note_2 + 0.8*delta;
