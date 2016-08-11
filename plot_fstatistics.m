X = [ note_x(:) ];                        % data
K = 5;

for k = 2: K
    [idx(:,k),C] = kmeans(X,k,'Distance','cityblock',...     % k clustering (idx and C varying after loop)
        'Replicates',5,'Options',statset('Display','final'));
    uv = unique(idx(:,k));                                    % list of clutsters [1 .. k]
    n  = histc(idx(:,k),uv);                                  % number of elements in each cluster (vector)
    
    for i = 1 : k     % inter community
        community_index{i} = find(idx(:,k)==i)';
        community{i} = note_x (community_index{i});   % community partition
        
        num_F_(i) =( n(i) * (distance(mean(community{i}), mean(note_x), 2))^2 ) / (k - 1);              % distance INTER - community
        
        for j = 1 : n(i)     % intra community
            den_F_d(j) = distance( community{i}(j), mean(community{i}), 2)^2 / (length(kx) - k);        % distance INTRA - community j
        end
        
        den_F_(i) = sum(den_F_d);
        clearvars den_F_d;
    end
    
    num_F(k) = sum(num_F_);
    den_F(k) = sum(den_F_);
    F5(k) = num_F(k) / den_F(k);             % F-statistics notation
    
    %   - plot centroïds and clusters -
    figure(k);
    for l = 1:k   
     plot(community{l},'.','MarkerSize',12);  
     hold on
     plot (C(l,1),'xk','MarkerSize',15);
    end
    hold off
    clearvars idx C n uv community_index community num_F_ den_F_ ;
end

clearvars -except F5
load(strcat('3801060_0007m.mat'));
dt0 = 8e-3;
val(isnan(val)) = [];

t0 = (1:length(val)) * dt0;            % timeline
s0 = val(1,1:length(val));

s0  = (s0  - mean(s0 ))/sqrt(var(s0));        % rescale s on 0 (standard score of signal)

%   - Timeline, noise, integration, quantization -
dt = 1/10;                           % sampling time: dt >> dt0
t_int = dt * (1/3);                  % integration time: dt0 <= t_int < dt
quant = 0.1;                         % LSB: vertical step

[t,s] = integration(t0,s0,dt0,dt,t_int,quant,0);

%   - Derivative -
d = s(2:end) -  s(1:end-1);
%td = (  t(2:end) +  t(1:end-1) ) / 2;      % timeline of derivative shifted by t_sample/2
td = t(2:end);

kx = d > 0;
kx = find(kx(1:end-1) & ~kx(2:end));       % k_{x}:index where d > 0; d( k_{x} + 1 ) <= 0

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

%   - Peaks notation -
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
    note_1(k) = sx(k) - ( sx(k+1) + sx(k-1) )/2;                                % average peak value
    note_3(k) = delta(k) - ( delta(k+1) + delta(k-1) )/2;                       % average peak to peak value
end

note_x = 0.1*note_1 + 0.1*note_2 + 0.8*delta;

X = [ note_x(:) ];                        % data
K = 5;

for k = 2: K
    [idx(:,k),C] = kmeans(X,k,'Distance','cityblock',...     % k clustering (idx and C varying after loop)
        'Replicates',5,'Options',statset('Display','final'));
    uv = unique(idx(:,k));                                    % list of clutsters [1 .. k]
    n  = histc(idx(:,k),uv);                                  % number of elements in each cluster (vector)
    
    for i = 1 : k     % inter community
        community_index{i} = find(idx(:,k)==i)';
        community{i} = note_x (community_index{i});   % community partition
        
        num_F_(i) =( n(i) * (distance(mean(community{i}), mean(note_x), 2))^2 ) / (k - 1);              % distance INTER - community
        
        for j = 1 : n(i)     % intra community
            den_F_d(j) = distance( community{i}(j), mean(community{i}), 2)^2 / (length(kx) - k);        % distance INTRA - community j
        end
        
        den_F_(i) = sum(den_F_d);
        clearvars den_F_d;
    end
    
    num_F(k) = sum(num_F_);
    den_F(k) = sum(den_F_);
    F7(k) = num_F(k) / den_F(k);             % F-statistics notation
    
    %   - plot centroïds and clusters -
    figure(k);
    for l = 1:k   
     plot(community{l},'.','MarkerSize',12);  
     hold on
     plot (C(l,1),'xk','MarkerSize',15);
    end
    hold off
    clearvars idx C n uv community_index community num_F_ den_F_ ;
end