Name = '3900679m';      % row 5

load(strcat(Name, '.mat'));
fid = fopen(strcat(Name, '.info'), 'rt');
fgetl(fid);
fgetl(fid);
fgetl(fid);
[interval] = sscanf(fgetl(fid), 'Sampling frequency: %f Hz  Sampling interval: %f sec');
interval = interval(2);              % data acquisition rate (interval = 1/f_spl_u = 0.5903 ms in practice)

fclose(fid);

t0 = (1:length(val)) * interval;            % timeline
s0 = val(5,1:length(val));
s0  = (s0  - mean(s0 ))/sqrt(var(s0 ));        % rescale s on 0 (standard score of signal)

%   - Timeline, noise, integration, quantization -
dt = 1/10;                           % sampling time: dt >> interval
t_int = dt * (1/3);                  % integration time: interval <= t_int < dt
quant = .1;                        % LSB: vertical step

[t,s] = integration(t0,s0,interval,dt,t_int,quant,0);

%   - Derivative -
d = s(2:end) -  s(1:end-1);
%td = (  t(2:end) +  t(1:end-1) ) / 2;      % timeline of derivative shifted by t_sample/2
td = t(2:end);

kx = d > 0;
kx = find(kx(1:end-1) & ~kx(2:end));       % k_{x}:index where d > 0; d( k_{x} + 1 ) <= 0

%   - Local maxima sx, maximum slope around sx -
[tx,sx, dhi,dlo, tx_N,sx_N, note_x] = peaks_processing(t,s,kx);

kmax_init = 4;
[clust_index,  ~,~,  ~,~,  kmax, diff] = agglo_clustering(note_x, tx, kmax_init);

if diff(2,2) >= 1    % EMPIRICAL: no clustering if 2-clustering clusters are too close
    
    % Initialization with outliers removed
    [clust_index,  clust_note_x,mean_clust,  clust_tx,clust_periodicity,  kmax, diff] = agglo_clustering(note_x, tx, kmax_init);
    
    % Search for best number of clusters
    div = 2;            % EMPIRICAL: merge clusters that are too close (    min(mean_clust difference) <= (max(mean_clust) - min(mean_clust)) / div
    while min( diff(2:end,kmax) ) <= ( max(mean_clust(:,kmax)) - min(mean_clust(:,kmax)) )/div && kmax >= 3
        
        kmax = kmax - 1;
        [clust_index,  clust_note_x,mean_clust,  clust_tx,clust_periodicity,  kmax, diff] = agglo_clustering(note_x, tx, kmax);
        
    end
    
    [~,clust_major_index] = max(mean_clust(:,kmax));
    kx_major = clust_index{clust_major_index,kmax};
    tx_major = tx(kx_major);
    sx_major = sx(kx_major);
    
else
    clust_note_x = nan;
    clust_tx = nan;
    [clust_periodicity(1),clust_periodicity(2),clust_periodicity(3)] = periodicity(tx);
    kmax = 1;
    tx_major = tx;
    sx_major = sx;
end

plot( kron(tx,[1 1 1]) , kron(dlo,[1 0 nan]) + kron(dhi,[0 1 nan]), '-c');       % link note_2
hold on
plot( tx_major , sx_major   , 'pk','MarkerSize',15);
plot( tx , sx   , 'dr','MarkerSize',12);
plot(tx, dhi,'c^','MarkerSize',12);
plot(tx, dlo,'cv','MarkerSize',12);
plot( tx,sx_N, 'dr','MarkerSize',12);
plot(kron(tx,[1 1 1]), kron(sx_N,[1 0 nan]) + kron(sx,[0 1 nan]),'r-');
hold off

% plot clust_note_x
for i = 1 : kmax
    figure(2);
    subplot(2,1,1);
    plot(clust_note_x{i,kmax} , '.');
    hold on
end
Legend=cell(kmax,1);
for iter=1:kmax
    Legend{iter}=strcat('cluster ', num2str(iter));
end
legend(Legend);

hold off
subplot(2,1,2);
plot( note_x, '.');

base_array = cellfun(@length,clust_tx);
base_max = max (base_array(:,kmax));
base = [1:base_max];

% plot clust_periodicity
data = nan(base_max,kmax);

for i = 1 : kmax

    data(1:base_array(i,kmax),i) = clust_tx{i,kmax};
    figure(3);
    tx_disp(i) = plot(base,data(:,i),'.');
    hold on
end

Legend=cell(kmax,1);
for iter=1:kmax
    Legend{iter}=strcat('cluster ', num2str(iter),': T = ', num2str(clust_periodicity{iter,kmax}(1)), '; eps = ', num2str(clust_periodicity{iter,kmax}(2)), '; R = ', num2str(clust_periodicity{iter,kmax}(3)));
end
legend(Legend);

title('Linear regression of t_{x,k}');
xlabel('k');
ylabel('t_{x,k}, s');
hold off


%%
delta = delta_tx(tx);
[T,eps,R,plot_reg] = periodicity(tx);

plot_reg;

eps_note_x = 1;
eps_per = 0.05;

med = median(note_x);
med_ = find(note_x == med);

% if var([note_x(1) note_x(2)]) < 1
%     %periodic_note_x(1) = note_x(1); periodic_note_x(2) = note_x(2);
%     for k = 2 : length(kx)-1
%     var( [tx(k+1)-tx(k) tx(k)-tx(k-1)]) < 1
%         clust_per(k-1:k+1) = note_x(k-1:k+1);
%     
%         
%     end    
% else
%     periodic_note_x(1) = note_x_(1);
%     
%     
% end
