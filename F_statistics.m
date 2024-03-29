% Ivan NY HANITRA - Master thesis
%       -- Compute the optimal number of clusters using F-statistics --

%   - Load file and data -
%Name = '3801060_0007m';   % row 1
Name = '3900497m';     % row 6 - BPM = 95

load(strcat(Name, '.mat'));
fid = fopen(strcat(Name, '.info'), 'rt');
fgetl(fid);
fgetl(fid);
fgetl(fid);
[interval] = sscanf(fgetl(fid), 'Sampling frequency: %f Hz  Sampling interval: %f sec');
interval = interval(2);              % data acquisition rate (interval = 1/f_spl_u = 0.5903 ms in practice)

fclose(fid);

t0 = (1:length(val)) * interval;            % timeline
s0 = val(6,1:length(val));
s0  = (s0  - mean(s0 ))/sqrt(var(s0 ));        % rescale s on 0 (standard score of signal)

%   - Timeline, noise, integration, quantization -
dt = 1/10;                           % sampling time: dt >> interval
t_int = dt * (1/3);                  % integration time: interval <= t_int < dt
quant = 0.1;                        % LSB: vertical step

[t,s] = integration(t0,s0,interval,dt,t_int,quant,0);

%   - Derivative, local maxima sx, maximum slope around sx -
d = s(2:end) -  s(1:end-1);
td = (  t(2:end) +  t(1:end-1) ) / 2;

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
    kx_index(1) = nan;
    sx_N(1) = nan;
    tx_N(1) = nan;
    
    for k = 2:length(kx)
        kx_index(k) = max( find( kx_n < kx(k) ) );
        sx_N(k) = s(kx_n( kx_index(k)) + 1);
        tx_N(k) = td(kx_n( kx_index(k) )) + (td(kx_n( kx_index(k) )+1)-td(kx_n( kx_index(k) ))) .* d(kx_n( kx_index(k) ))./(d(kx_n( kx_index(k) ))-d(kx_n( kx_index(k) )+1));
    end
    
end

delta = sx - sx_N;

%   - Peaks notation
note_1 = sx;
for k = 2:length(kx)-1
    note_1(k) = 2*sx(k) - sx(k+1) - sx(k-1);  % average peak value (doubled)
end

note_2 = dhi - dlo;                           % maximum slope difference around peak

note_x = (note_1 + note_2 + 1*delta)/3;

%   - Plot signals -
% figure(1);
% plot(t0, s0,'k-','MarkerSize',8,'LineWidth',.5);               % siganl s
% hold on
% plot(t, s,'ko--','MarkerSize',10,'LineWidth',1);     % sampled signal s_n
% plot(td, d,'g--','MarkerSize',10,'LineWidth',1);     % derivative of s_n 
% plot(tx,sx,'rp','MarkerSize',12,'LineWidth',2);
% plot(tx_N,sx_N,'mp','MarkerSize',12,'LineWidth',2);
% plot(kron(tx,[1 1 1]), kron(sx_N,[1 0 nan]) + kron(sx,[0 1 nan]),'r-');
% hold off

figure(7);
plot(note_x,'r.');

%%
%   - k-means ++ for centro�d research -
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
    F(k) = num_F(k) / den_F(k);             % F-statistics notation
    
    %   - plot centro�ds and clusters -
    figure(k);
    for l = 1:k   
     plot(community{l},'.','MarkerSize',12);  
     hold on
     plot (C(l,1),'xk','MarkerSize',15);
    end
    hold off
    clearvars idx C n uv community_index community num_F_ den_F_ ;
end


figure(1);
plot(F);



