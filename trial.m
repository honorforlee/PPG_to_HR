% Ivan NY HANITRA - Master thesis
%       -- Clustering peaks, discrimination, compute PPG frequency --

%   - Load file and data -
Name = '3801060_0007m';   % row 1
%Name = '3900497m';     % row 6 - BPM = 95

load(strcat(Name, '.mat'));
fid = fopen(strcat(Name, '.info'), 'rt');
fgetl(fid);
fgetl(fid);
fgetl(fid);
[interval] = sscanf(fgetl(fid), 'Sampling frequency: %f Hz  Sampling interval: %f sec');
interval = interval(2);              % data acquisition rate (interval = 1/f_spl_u = 0.5903 ms in practice)

fclose(fid);

t = (1:length(val)) * interval;            % timeline
s = val(1,1:length(val));
s  = (s  - mean(s ))/sqrt(var(s ));        % rescale s on 0 (standard score of signal)

%   - Timeline, noise, integration, quantization -
dt = 1/10;                           % sampling time: dt >> interval
t_int = dt * (1/3);                  % integration time: interval <= t_int < dt
quant = 1e-4;                        % LSB: vertical step

[t_spl,s_spl] = integration(t,s,interval,dt,t_int,quant,0);

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

kx_n = d_spl < 0;
kx_n = find(kx_n(1:end-1) & ~kx_n(2:end));

tx_n = td_spl(kx_n) + (td_spl(kx_n+1)-td_spl(kx_n)) .* d_spl(kx_n)./(d_spl(kx_n)-d_spl(kx_n+1)); 
sx_n = s_spl(kx_n+1);

delta = sx;
if kx(1) < kx_n(1)
    delta(1) = sx(1) - s_spl(1);
    delta_plot = kron(sx_n(2:length(kx)),[1 0 nan]) + kron(sx,[0 1 nan]);
    for k = 1:length(kx)-1
        delta(k+1) = sx(k+1) - sx_n(k);
    end
else
    delta_plot = kron(sx_n(1:length(kx)),[1 0 nan]) + kron(sx,[0 1 nan]);
    for k = 1:length(kx)
        delta(k) = sx(k) - sx_n(k);
    end
end

plot(t, s,'k-','MarkerSize',8,'LineWidth',.5);               % siganl s
hold on
plot(t_spl, s_spl,'ko--','MarkerSize',10,'LineWidth',1);     % sampled signal s_n
plot(td_spl, d_spl,'g--','MarkerSize',10,'LineWidth',1);     % derivative of s_n 
plot(tx,sx,'rd','MarkerSize',12,'LineWidth',2);
plot(tx_n,sx_n,'bd','MarkerSize',12,'LineWidth',2);
plot(kron(tx,[1 1 1]),delta_plot , 'r-','LineWidth',3);                     % link delta
hold off



%%
%   - Peaks notation
note_1 = sx;
for k = 2:length(kx)-1
    note_1(k) = 2*sx(k) - sx(k+1) - sx(k-1);  % average peak value (doubled)
end

note_2 = dhi - dlo;                           % maximum slope difference around peak

note_3 = zeros(1,length(kx));                 % peak prominence notation
for k = 1:length(kx)
    i = kx(k) - 1;
    if (i>0 && d_spl(i) > 0)
        note_3(k) = sx(k) - s_spl(i); 
        i = i-1;
    elseif (i>0 && d_spl(i) <= 0)
        note_3(k) = sx(k) - s_spl(i); 
    end
end

note_weighted = note_1 ./ note_2;

% major = kx;
% for k = 1 : length(kx)
%     if abs(note_weighted(k)) <= 1
%         major(k)=kx(k);
%     else
%         major(k) = 0;
%     end
% end
% 
% major_index = find(major);

[T,eps,R_sq,plot_reg] = periodicity(tx);
%[T_2,eps_2,R_sq_2,plot_reg_2] = periodicity(tx(1:2:end));

kx_ = kx;

for k = 1:length(tx) - 1            % discard minor peaks with a frequency f > 3.5 Hz (BPM_max = 210)
    if (tx(k+1)-tx(k)) <= (1/3.5)
        if note_1(k+1) < note_1(k)
            kx_(k+1) = 0;
        else
            kx_(k) = 0;
        end
    end
end

% for k = 1:length(kx) - 1            % discard minor peaks with a frequency f > 3.5 Hz (BPM_max = 210)
%     if (kx(k+1)-kx(k)) <= floor((1/3.5)/dt)
%         kx_(k) = 0;
%     end
% end

kx_major_index= find(kx_(1:end));
kx_major = kx_(kx_major_index);

sx_major = s_spl(kx_major + 1);
tx_major = tx(kx_major_index);

for k = 1 : length(tx_major)-1
diff(k)=tx_major(k+1)-tx_major(k);
end

ax(1) = subplot(2,1,1);
plot(diff,'rp');

ax(2) = subplot(2,1,2);
plot(t, s,'k-','MarkerSize',12,'LineWidth',1);               % siganl s
hold on
plot(t_spl, s_spl,'ko--','MarkerSize',10,'LineWidth',1);     % sampled signal s_n
plot(td_spl, d_spl,'g--','MarkerSize',12,'LineWidth',1);     % derivative of s_n 
plot(tx_major,sx_major,'rp','MarkerSize',15,'LineWidth',3);  % major peaks
hold off

% linkaxes(ax(1:2),'xy');
% axis(ax,[0 60 -5 5]);

%%
%  - k-means clustering of peaks according to sx and note_2 -

X = [ note_weighted(:) ];                        % data

[idx,C] = kmeans(X,2,'Distance','cityblock',...     % 2 clusters created: minor/major peaks
    'Replicates',5,'Options',statset('Display','final'));  % initialize the replicates 5 times, separately using k-means++ algorithm, choose best arrangement and display final output

cluster1 = find(idx==1)';
cluster2 = find(idx==2)';

if note_weighted(cluster1(1)) > note_weighted(cluster2(1))                % assign major peak cluster
    major_index = cluster1;
else
    major_index = cluster2; 
end

tx_major = tx(major_index);                       
sx_major = sx(major_index);

subplot(2,1,1)
plot(t, s,'k-','MarkerSize',12,'LineWidth',1);               % siganl s
hold on
plot(t_spl, s_spl,'ko--','MarkerSize',10,'LineWidth',1);     % sampled signal s_n
plot(td_spl, d_spl,'g--','MarkerSize',12,'LineWidth',1);     % derivative of s_n 
plot(tx,sx,'cp','MarkerSize',12,'LineWidth',1);              % note_1 
plot(tx,note_weighted,'rd','MarkerSize',13,'LineWidth',1);     % note_weighted
plot(tx_major,sx_major,'rp','MarkerSize',15,'LineWidth',3);  % major peaks

plot(tx,dhi,'c^' ,'MarkerSize',12,'LineWidth',1);            % d_max_l
plot(tx,dlo,'cv','MarkerSize',12,'LineWidth',1);             % d_max_r
plot(kron(tx,[1 1 1]) , kron(dlo,[1 0 nan]) + kron(dhi,[0 1 nan]) , 'c-','MarkerSize',12,'LineWidth',1);    % note_2

title('Peaks discrimination for heart rate monitoring');
xlabel('Time, s');
ylabel('arbitrary units');
%legend('s: original signal','s_n: sampled signal','d_n: derivative of s_n','Peaks','Major peaks','Location','southoutside');
hold off

subplot(2,1,2)
plot(X(idx==1),'r.','MarkerSize',12);      % cluster 1 correspoding sx,delta_note2 (minor)
hold on
plot(X(idx==2),'b.','MarkerSize',12);      % cluster 2 correspoding sx,delta_note2 (major)
plot(C,'kx','MarkerSize',15,'LineWidth',3);  % plot centroids
title ('Cluster Assignments and Centroids');
legend('Cluster1','Cluster2','Centroids','Location','NW');
xlabel ('note_{1/2}, arbitrary units');
hold off


%%
%   - plot(2 notes) -
figure(1);
plot(X(idx==1,1),X(idx==1,2),'r.','MarkerSize',12);      % cluster 1 correspoding sx,delta_note2 (minor)
hold on
plot(X(idx==2,1),X(idx==2,2),'b.','MarkerSize',12);      % cluster 2 correspoding sx,delta_note2 (major)
plot(C(:,1),C(:,2),'kx','MarkerSize',15,'LineWidth',3);  % plot centroids
 
title ('Cluster Assignments and Centroids');
legend('Cluster1','Cluster2','Centroids','Location','NW');
xlabel ('note_1, a.u');
ylabel ('{note_2}, a.u');
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

%%

%   - {kx} periodicity -
%[T,eps,R_sq,plot_reg] = periodicity(tx);

tbl = table([1:length(tx)]', tx','VariableNames',{'k','tx'});
mdl = fitlm(tbl,'tx~k');
F_stat = anova(mdl);                        % analyse of variance
F = F_stat.F(1);                            % F = MeanSq(xi)/MeanSq(Error) with MeanSq = SumSq/DF) (DF(xi)=1 , DF(error)=length(kx)-2)


% %[T,eps,R_sq,plot_reg] = periodicity(tx);      % peaks periodicity                          
% %plot_reg;
% 
% random_kx = randsample(  kx, length(kx)/2   );
% random_kx = sort(random_kx);
% random_tx = td_spl(random_kx) + (td_spl(random_kx+1)-td_spl(random_kx)) .* d_spl(random_kx)./(d_spl(random_kx)-d_spl(random_kx+1)); 
% [T,eps,R_sq,plot_reg] = periodicity(random_tx); 
% plot_reg;

%   - Compute PPG frequency -
tx_major = tx(major_index);                       
sx_major = sx(major_index);

for k = 1 : length(tx_major) - 1
    
    dtx_major(k)= tx_major(k+1) - tx_major(k);        % time interval between major peaks
    
end

freq_ppg = 1 ./ (mean(dtx_major));
BPM = round(60 * freq_ppg)
note_P = var(1./dtx_major);


%%
%X = [ note_1(:),note_2(:),note_P(:) ];  
rng default;  % For reproducibility
X = [gallery('uniformdata',[10 3],12);...
    gallery('uniformdata',[10 3],13)+1.2;...
    gallery('uniformdata',[10 3],14)+2.5];

T = clusterdata(X,'distance','cityblock','maxclust',3);
find(T==5)

scatter3(X(:,1),X(:,2),X(:,3),100,T,'filled')


