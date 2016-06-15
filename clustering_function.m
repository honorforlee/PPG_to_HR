function [freq_ppg, BPM, kx, major_index, tx_major, sx_major,note_P] = clustering_function(val, interval, dt, t_int, quant) 

t = (1:length(val)) * interval;              % timeline
s = val(1,1:length(val));
s  = (s  - mean(s ))/sqrt(var(s ));          % rescale s on 0 (standard score of signal)

%   - Timeline, noise, integration, quantization -
dt = 1/10;                           % sampling time: dt >> interval
t_int = dt * (1/3);                  % integration time: interval <= t_int < dt
quant = 1e-4;                        % LSB: vertical step
t_spl = t(1):dt:t(end); 

% Noise
for k = 1:length(t_spl)-1
frameNoise (:,k) = [ floor(t_spl(k)/interval) :  floor(t_spl(k)/interval) + floor(dt/interval) ]; 
end
noise= random('Normal',mean(s(frameNoise)),std(s(frameNoise)),1,length(frameNoise));                     % Gaussian distribution (model thermal noise of finite BW)

% Integration
for k = 2:length(t_spl)
index = min (length(    [floor( (t_spl(k-1)+ dt-t_int)/interval ): floor( t_spl(k)/ interval ) ]    ));
end
index = index-1;
for k = 2:length(t_spl)
frameInteg(:,k-1) = [ floor( t_spl(k)/ interval ) - index : floor( t_spl(k)/ interval ) ];    
frameInteg_(:,k-1)= s(frameInteg(:,k-1)) ;
end
frameInteg_ = vertcat(frameInteg_,noise);

s_spl(1) = s(1);
for k = 2:length(t_spl)
    s_spl(k) = mean(frameInteg_(:,k-1));               % sampled signal = average of Nint last values + noise during dt
end    

s_spl = quant*floor(s_spl/quant);                      % quantization

%   - with time subdivided in indexes -
% subels = (1:round(dt/interval):length(t));
% t_spl = t(subels);                           % sample timeline
% 
% % Noise
% frameNoise = (0:round(dt/interval))';
% frameNoise = bsxfun(@minus, subels, frameNoise);
% frameNoise_zero = find (frameNoise <= 0);
% frameNoise(frameNoise_zero) = 1;
% 
% noise = random('Normal',mean(s(frameNoise)),std(s(frameNoise)),1,length(subels));      % Gaussian distribution (model thermal noise of finite BW)
% 
% % Integration
% frameInteg = (0:round(t_int/interval))';
% frameInteg = bsxfun(@minus, subels, frameInteg);
% frameInteg_zero = find (frameInteg <= 0);
% frameInteg(frameInteg_zero) = 1;                       % t_int < dt
% 
% s_spl = mean( vertcat(s(frameInteg), noise) );         % sampled signal = average of Nint last values + noise during dt
% 
% s_spl = quant*floor(s_spl/quant);                      % quantization


%   - Derivative, local maxima sx, maximum slope around sx -
d_spl = s_spl(2:end) -  s_spl(1:end-1);
td_spl = (  t_spl(2:end) +  t_spl(1:end-1) ) / 2;

kx = d_spl > 0;
kx = find(kx(1:end-1) & ~kx(2:end));       % k_{x}:index where d_spl > 0; d_spl( k_{x} + 1 ) <= 0

sx = s_spl(kx+1);                          % local maxima
tx = td_spl(kx) + (td_spl(kx+1)-td_spl(kx)) .* d_spl(kx)./(d_spl(kx)-d_spl(kx+1));      % linear interpolation of dhi and dho to get tx (@zero crossing)

dhi = d_spl(kx);
dlo = d_spl(kx+1);

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

%   - Compute PPG frequency -
tx_major = tx(major_index);                       
sx_major = sx(major_index);

for k = 1 : length(tx_major) - 1
    
    dtx_major(k)= tx_major(k+1) - tx_major(k);        % time interval between major peaks
    
end

freq_ppg = 1 ./ (mean(dtx_major));
BPM = round(60 * freq_ppg);
note_P = var(1./dtx_major);
