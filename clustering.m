% Ivan NY HANITRA - Master thesis
%       -- Clustering peaks, discrimination, compute PPG frequency --

%   - Load file and data -
%Name = '3987834m';     % BPM = 78
Name = '3801060_0007m'; % BPM = 95

load(strcat(Name, '.mat'));
fid = fopen(strcat(Name, '.info'), 'rt');
fgetl(fid);
fgetl(fid);
fgetl(fid);
[interval] = sscanf(fgetl(fid), 'Sampling frequency: %f Hz  Sampling interval: %f sec');
interval = interval(2);              % data acquisition rate (interval = 1/f_spl_u = 0.5903 ms in practice)

fclose(fid);

t = (1:length(val)) * interval;              % timeline
s = val(1,1:length(val));
s  = (s  - mean(s ))/sqrt(var(s ));          % rescale s on 0 (standard score of signal)

%   - Timeline, noise, integration, quantization -
dt = 1/10;                           % sampling time: dt >> interval
t_int = dt * (1/3);                  % integration time: interval <= t_int < dt
quant = 1e-4;                        % LSB: vertical step

subels = (1:round(dt/interval):length(t));
t_spl = t(subels);                           % sample timeline

% Noise
frameNoise = (0:round(dt/interval))';
frameNoise = bsxfun(@minus, subels, frameNoise);
frameNoise_zero = find (frameNoise <= 0);
frameNoise(frameNoise_zero) = 1;

noise= random('Normal',mean(s(frameNoise)),std(s(frameNoise)),1,length(subels));                     % Gaussian distribution (model thermal noise of finite BW)

% for k = 1 : length(frameNoise)
%  
% noise_shot (k) = sum(  poisspdf( fliplr(frameNoise(:,k)') , mean( abs( s( fliplr(frameNoise(:,k)')) )) ) ); % Poisson statistics (model shot noise of Photodiode: independant random events)
% 
% end

% Integration
frameInteg = (0:round(t_int/interval))';
frameInteg = bsxfun(@minus, subels, frameInteg);
frameInteg_zero = find (frameInteg <= 0);
frameInteg(frameInteg_zero) = 1;                       % t_int < dt

%s_spl = mean( vertcat(s(frameInteg),noise) );          % sampled signal = average of Nint last values + noise during dt
s_spl = mean( s(frameInteg) );                          % signal with no noise

s_spl = quant*floor(s_spl/quant);                      % quantization

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

%   - Peaks notation
note_1 = sx;
for k = 2:length(kx)-1
    note_1(k) = 2*sx(k) - sx(k+1) - sx(k-1);  % average peak value (doubled)
end

note_2 = dhi - dlo;                           % maximum slope difference around peak

[T,eps,R_sq,plot_reg] = periodicity(tx);      % peaks periodicity                          
plot_reg;

%   - k-means clustering of peaks according to sx and note_2 -
X = [ note_1(:),note_2(:) ];                        % data

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
BPM = round(60 * freq_ppg)
note_P = var(1./dtx_major);

%   - Plots -
figure (1);
plot(X(idx==1,1),X(idx==1,2),'r.','MarkerSize',12)  % cluster 1 correspoding sx,delta_note2 (minor)
hold on
plot(X(idx==2,1),X(idx==2,2),'b.','MarkerSize',12)  % cluster 2 correspoding sx,delta_note2 (major)
plot(C(:,1),C(:,2),'kx',...                         % plot centroids
     'MarkerSize',15,'LineWidth',3)
 
title 'Cluster Assignments and Centroids'
legend('Cluster1','Cluster2','Centroids',...
       'Location','NW')
xlabel ('Peak amplitude, a.u')
ylabel ('delta_{note2}, a.u');
hold off

figure(2);
plot(t, s,'k-','MarkerSize',12,'LineWidth',1);               % siganl s
hold on
plot(t_spl, s_spl,'ko--','MarkerSize',10,'LineWidth',1);     % sampled signal s_n
plot(td_spl, d_spl,'g--','MarkerSize',12,'LineWidth',1);     % derivative of s_n 
plot(tx_major,sx_major,'rp','MarkerSize',15,'LineWidth',3);  % major peaks

title('Peaks discrimination for heart rate monitoring');
xlabel('Time, s');
ylabel('Arbitrary units');
legend('s: original signal','s_n: sampled signal','d_n: derivative of s_n','Major peaks','Location','northeastoutside');
hold off
