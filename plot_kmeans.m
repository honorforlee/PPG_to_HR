Name = '3899985_0005m';
load(strcat(Name, '.mat'));
fid = fopen(strcat(Name, '.info'), 'rt');
fgetl(fid); fgetl(fid); fgetl(fid);
dt0 = sscanf(fgetl(fid), 'Sampling frequency: %f Hz  Sampling interval: %f sec');
dt0 = dt0(2);
fgetl(fid); sig = fgetl(fid);

ppg = 0; ecg = 0; meas = 0;
while sig
    [k,signal,~,~,~] = strread(sig,'%d%s%f%f%s','delimiter','\t'); %#ok<DSTRRD>
    if strcmp(signal,'II');    ecg = k; end
    if strcmp(signal,'PLETH'); ppg = k; end
    if strcmp(signal,'MEAS'); meas = k; end
    sig = fgetl(fid);
end
fclose(fid);

if meas == 0
    val(isnan(val)) = [];
    t0 = (1:length(val)) * dt0;            % timeline
    s0 = val(ppg,1:length(val));
    
else
    Vout(isnan(Vout)) = [];
    t0 = (1:length(Vout)) * dt0;            % timeline
    s0 = Vout';
end

s0  = (s0  - mean(s0 ))/sqrt(var(s0));        % rescale s on 0 (standard score of signal)

%   - Timeline, noise, integration, quantization -
dt = 1/20;                           % sampling time: dt >> dt0
t_int = dt * (1/3);                  % integration time: dt0 <= t_int < dt
quant = 0.1;                         % LSB: vertical step

[t,s] = integration(t0,s0,dt0,dt,t_int,quant,0);

%  - Peaks identification -
[kx,tx,sx, dhi,dlo, td,d, kx_n,tx_N,sx_N, note_x] = signal_peaks(t,s);

frame_init =45; frame_end = 50;

index_x = find(tx >= frame_init & tx <= frame_end);
sx_N_frame = sx_N(index_x);
dhi_frame = dhi(index_x);
dlo_frame = dlo(index_x);

index0 = find(t0 >= frame_init & t0 <= frame_end);
t0_frame = t0(index0);
s0_frame = s0(index0);

index = find(t >= frame_init & t <= frame_end);
t_frame = t(index);
s_frame = s(index);

[kx_frame,tx_frame,sx_frame,note_x_frame] = frame_select(kx,tx,sx,note_x, frame_init,frame_end);

kx = kx_frame; tx = tx_frame; sx = sx_frame; note_x = note_x_frame;

%%
% X = [ note_1(:),note_2(:),delta(:) ];                        % data
% 
% [idx,C] = kmeans(X,2,'Distance','cityblock',...     % 2 clusters created: minor/major peaks
%     'Replicates',5,'Start','plus','Options',statset('Display','final'));  % initialize the replicates 5 times, separately using k-means++ algorithm, choose best arrangement and display final output
% 
% plot3(X(idx==1,1),X(idx==1,2),X(idx==1,3),'r.','MarkerSize',15);  % cluster 1 correspoding sx,delta_note2 (minor)
% hold on
% plot3(X(idx==2,1),X(idx==2,2),X(idx==2,3),'b.','MarkerSize',15);  % cluster 2 correspoding sx,delta_note2 (major)
% plot3(C(:,1),C(:,2),C(:,3),'kx','MarkerSize',20,'LineWidth',3);
% 
% title ('Cluster Assignments and Centroids','FontSize',30,'FontWeight','Bold');
% legend('Cluster1','Cluster2','Centroids','Location','NW');
% xlabel('note_1','FontSize',25,'FontWeight','Bold');
% ylabel('note_2','FontSize',25,'FontWeight','Bold');
% zlabel('note_3','FontSize',25,'FontWeight','Bold');
% grid on
% hold off

X = [note_x(:)];
[idx,C] = kmeans(X,2,'Distance','cityblock',...     % 2 clusters created: minor/major peaks
     'Replicates',5,'Start','plus','Options',statset('Display','final'));  % initialize the replicates 5 times, separately using k-means++ algorithm, choose best arrangement and display final output 
%%
hold on
xlim([-1,5]);ylim([-.2,2])

%arrow
[arrowX,arrowY]=dsxy2figxy([-1,5],[0,0]);
annotation('arrow',arrowX,arrowY)

%crosses

plot(X(idx==1),0,'or','markersize',10);
hold on
plot(X(idx==2),0,'ob','markersize',10);
plot(C,0,'kx','MarkerSize',20,'LineWidth',3);
xlabel ('note_x','FontSize',25,'FontWeight','Bold');

%pipes
p=[0.5,0.65];
text(p,[0,0],'$$\vert$$','interpreter','latex')

grid on
axis on
print('-depsc','arrowFigure')
%%
cluster1 = find(idx==1)';
cluster2 = find(idx==2)';

if sx(cluster1(1)) > sx(cluster2(1))                % assign major peak cluster
    major_index = cluster1;
else
    major_index = cluster2;
end

tx_major = tx(major_index);
sx_major = sx(major_index);
plot(t_frame,s_frame,'og');
hold on
plot(t0_frame,s0_frame,'--k');
plot(tx_major,sx_major,'pr','MarkerSize',25);

 
%%

plot(X(idx==1,1),X(idx==1,2),X(idx==1,3),'r.','MarkerSize',15)  % cluster 1 correspoding sx,delta_note2 (minor)
hold on
plot(X(idx==2,1),X(idx==2,2),X(idx==2,3),'b.','MarkerSize',15)  % cluster 2 correspoding sx,delta_note2 (major)
plot(C(:,1),C(:,2),'kx',...                         % plot centroids
    'MarkerSize',20,'LineWidth',3)

title ('Cluster Assignments and Centroids','FontSize',30,'FontWeight','Bold');
legend('Cluster1','Cluster2','Centroids','Location','NW');
xlabel ('note_1','FontSize',25,'FontWeight','Bold');
ylabel ('note_3','FontSize',25,'FontWeight','Bold');
hold off

cluster1 = find(idx==1)';
cluster2 = find(idx==2)';

if sx(cluster1(1)) > sx(cluster2(1))                % assign major peak cluster
    major_index = cluster1;
else
    major_index = cluster2;
end

for k = 1:length(major_index)-1
    
    freq = mean( major_index(k+1) - major_index);
    
end
%%
X = [ note_1(:),note_2(:),delta(:) ];                        % data

c = clusterdata(X,'linkage','ward','savememory','on','maxclust',2);
scatter3(X(:,1),X(:,2),X(:,3),36,c);


%%
plot(t0_frame,s0_frame,'--k');
hold on
plot(t_frame,s_frame,'og');
plot(tx_frame,sx_frame,'pb');