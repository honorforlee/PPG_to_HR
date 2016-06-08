function ppg_freq = clustering_function(val,interval) 

t = (1:length(val)) * interval;              % timeline
s = val(1,1:length(val));
s  = (s  - mean(s ))/sqrt(var(s ));          % rescale s on 0 (standard score of signal)

%   - Timeline, noise, integration, quantisization -
dt = 0.1;                            % sampling time: dt >> interval
t_int = dt * (1/3);                  % integration time: interval <= t_int < dt
quant = 1e-4;                        % LSB: vertical step

subels = (1:round(dt/interval):length(t));
t_spl = t(subels);                           % sample timeline
%s_spl = s(subels);

% Noise
frameNoise = (0:round(dt/interval))';
frameNoise = bsxfun(@minus, subels, frameNoise);
frameNoise_zero = find (frameNoise <= 0);
frameNoise(frameNoise_zero) = 1;

noise_therm = random('Normal',mean(s(frameNoise)),std(s(frameNoise)),1,length(subels));                     % Gaussian distribution (model thermal noise of finite BW)

for k = 1 : length(frameNoise)
 
noise_shot (k) = sum(  poisspdf( fliplr(frameNoise(:,k)') , mean( abs( s( fliplr(frameNoise(:,k)')) )) ) ); % Poisson statistics (model shot noise of Photodiode: independant random events)

end

% Integration
frameInteg = (0:round(t_int/interval))';
frameInteg = bsxfun(@minus, subels, frameInteg);
frameInteg_zero = find (frameInteg <= 0);
frameInteg(frameInteg_zero) = 1;                       % t_int < dt

s_spl = mean( vertcat(s(frameInteg),bsxfun(@plus, noise_therm , noise_shot) )   );          % sampled signal = average of Nint last values + noise during dt

s_spl = quant*floor(s_spl/quant);                      % quantisization

%   - Derivative, local maxima sx, maximum slope around sx -
d_spl = s_spl(2:end) -  s_spl(1:end-1);
td_spl = (  t_spl(2:end) +  t_spl(1:end-1) ) / 2;

kx = d_spl > 0;
kx = find(kx(1:end-1) & ~kx(2:end));       % k_{x}:index where d_spl > 0; d_spl( k_{x} + 1 ) <= 0

sx = s_spl(kx+1);                          % local maxima
tx = td_spl(kx) + (td_spl(kx+1)-td_spl(kx)) .* d_spl(kx)./(d_spl(kx)-d_spl(kx+1));      % linear interpolation of dhi and dho to get tx (@zero crossing)

dhi = d_spl(kx);
dlo = d_spl(kx+1);

% for k = 1:length(kx)                % search for maximum slope 0.5 s around s_max
%
%     for i = 1:floor(0.25/dt)
%         if (kx(k)+1 - i) > 0
%             if d_spl(kx(k)+1 - i) >= dhi(k)
%                 dhi(k) = d_spl(kx(k)+1 - i) ;
%             end
%         end
%         if (kx(k)+1 + i) < length(d_spl)
%             if d_spl(kx(k)+1 + i) <= dlo(k)
%                 dlo(k) = d_spl(kx(k)+1 + i);
%             end
%         end
%         i = i+1;
%     end
% end

for k = 1:length(kx)
    i = kx(k)-1;   while i > 0             && d_spl(i) >= dhi(k); dhi(k) = d_spl(i); i = i-1; end    % search for maximum positive slope at kx-
    i = kx(k)+2;   while i < length(d_spl) && d_spl(i) <= dlo(k); dlo(k) = d_spl(i); i = i+1; end    % search for maximmum negative slope at kx+
end

delta_note2 = dhi - dlo;

%   - k-means clustering of peaks according to sx and delta_note2 -
X = [ sx(:),delta_note2(:) ];                        % data

[idx,C] = kmeans(X,2,'Distance','cityblock',...     % 2 clusters created: minor/major peaks
    'Replicates',5,'Start','plus','Options',statset('Display','final'));  % initialize the replicates 5 times, separately using k-means++ algorithm, choose best arrangement and display final output

cluster1 = find(idx==1)';
cluster2 = find(idx==2)';

if sx(cluster1(1)) > sx(cluster2(1))                % assign major peak cluster
    major_index = cluster1;
else
    major_index = cluster2; 
end

tx_major_c = tx(major_index);                       
sx_major_c = sx(major_index);

%   - Note for major peaks frequency -
for k = 1:length(major_index) - 1
    
   freq_(k) = major_index(k+1) - major_index(k);    % major peak every freq_ peak
 
end

for k = 1 : length(tx) - 1
    
    dtx(k)= tx(k+1) - tx(k);                        % time interval between peaks
    
end

ppg_freq = 1./ (mean(freq_) * mean(dtx)) ;               % PPG frequency
