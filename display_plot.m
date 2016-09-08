Name = '3801060_0007m';
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
dt = 1/10;                           % sampling time: dt >> dt0
t_int = dt * (1/3);                  % integration time: dt0 <= t_int < dt
quant = 0.1;                         % LSB: vertical step

[t,s] = integration(t0,s0,dt0,dt,t_int,quant,0);

%  - Peaks identification -
[kx,tx,sx, dhi,dlo, td,d, kx_n,tx_N,sx_N, note_x] = signal_peaks(t,s);

frame_init =25; frame_end = 30;

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

[kx_major,tx_major,sx_major, T, warning] = min_variance(kx_frame,tx_frame,sx_frame, note_x_frame, 0.1);


%%
base= zeros(1,length(tx_frame));
grid = zeros(1,length(t0_frame));
note=note_x_frame;
%%
%   - Plots -
figure(2);
plot(tx_major,note,'p','Color',[1,0.5,0],'MarkerSize',25,'LineWidth',2);
hold on
plot(t0_frame,s0_frame,'--','Color',[0,0,0],'LineWidth',2);
plot(t_frame,s_frame,'o','Color',[0,0.5,0.5],'MarkerSize',10,'LineWidth',2);
plot(kron(tx_frame,[1 1 1]), kron(base,[1 0 nan]) + kron(note,[0 1 nan]),'-k','LineWidth',2);




%plot(t_frame,s_frame,'o','Color',[0,0.5,0.5],'MarkerSize',10,'LineWidth',2);


%plot(t0_frame,grid,'--k');

set(gca,'xtick',[]);
set(gca,'ytick',[]);

%legend({'Detected events','Discriminated events'},'FontSize',15,'Location','NE'); 


% plot( tx_frame , sx_frame , 'dc','MarkerSize',12);
% plot(tx_frame, dhi_frame,'^b','MarkerSize',10);
% plot(tx_frame, dlo_frame,'vb','MarkerSize',10);
% plot( kron(tx_frame,[1 1 1]) , kron(dlo_frame,[1 0 nan]) + kron(dhi_frame,[0 1 nan]), '-b');       % link note_2

% plot(tx_major , sx_major, 'pr','MarkerSize',20);
% plot( tx_frame,sx_N_frame, 'dc','MarkerSize',10);
% plot(kron(tx_frame,[1 1 1]), kron(sx_N_frame,[1 0 nan]) + kron(sx_frame,[0 1 nan]),'-c');

%%
plot(tx_frame,base,'ok','MarkerSize',15,'LineWidth',2);
set(gca,'xtick',[]);
set(gca,'ytick',[]);
%%
start = 25.5; 
finish = 26.75;

tx_ = tx_frame(tx_frame>start & tx_frame<finish); sx_ = sx_frame(tx_frame>start & tx_frame<finish); 

sx_avg = sx_;
for k = 2:length(tx_)-1
    sx_avg (k) = sx_(k) - 0.5*(sx_(k+1)+sx(k-1));
end

dlo_ = dlo_frame(tx_frame>start & tx_frame<finish); dhi_=dhi_frame(tx_frame>start & tx_frame<finish);
sx_N_ = sx_N_frame(tx_frame>start & tx_frame<finish); note_ = note_x_frame(tx_frame>start & tx_frame<finish);
tx_N_ = tx_N(tx_N>start & tx_N<finish);

t_ = t_frame(t_frame>start & t_frame<finish);    s_ = s_frame(t_frame>start & t_frame<finish);
t0_= t0_frame(t0_frame>start & t0_frame<finish); s0_= s0_frame(t0_frame>start & t0_frame<finish);
td_ = t_(2:end); d_ = s_(2:end) -  s_(1:end-1);

tx_maj = tx_major(tx_major>start & tx_major<finish); sx_maj = sx_major(tx_major>start & tx_major<finish);


null = zeros(1,length(tx_));
grid = zeros(1,length(t0_));

tick_p = 0.1*ones(1,length(tx_));
tick_m = -0.1*ones(1,length(tx_));

%% Plot events detection
%   - Plots -
% figure(2);
% plot(t0_,s0_,'--','Color',[0,0,0],'LineWidth',2);
% hold on
% plot(t_,s_,'o','Color',[0,0.5,0.5],'MarkerSize',15,'LineWidth',2);
 
% plot(kron(tx_,[1 1 1]), kron(sx_N_,[1 0 nan]) + kron(sx_,[0 1 nan]),'-r','LineWidth',2);
% plot( tx_ , sx_ , '^r','MarkerSize',20,'LineWidth',2);
% plot( tx_ , sx_N_ , 'vr','MarkerSize',20,'LineWidth',2);

% plot(tx_,note_,'p','Color',[1,0.5,0],'MarkerSize',25,'LineWidth',2);
% plot(kron(tx_,[1 1 1]), kron(null,[1 0 nan]) + kron(note_,[0 1 nan]),'-k','LineWidth',2);

plot(td_,d_,'--s','Color',[0.2,0.1,0.9],'LineWidth',2,'MarkerSize',10);
hold on
plot( kron(tx_,[1 1 1]) , kron(dlo_,[1 0 nan]) + kron(dhi_,[0 1 nan]), '-b','LineWidth',2);       % link note_2
plot(tx_, dhi_,'^b','MarkerSize',15,'LineWidth',2);
plot(tx_, dlo_,'vb','MarkerSize',15,'LineWidth',2);

% plot(tx_,sx_avg,'*r','MarkerSize',20,'LineWidth',2);
% plot(kron(tx_,[1 1 1]), kron(null,[1 0 nan]) + kron(sx_avg,[0 1 nan]),'--r','LineWidth',2);

%plot( tx_ , sx_ , '*r','MarkerSize',20,'LineWidth',2);
% plot( tx_N_ , sx_N_ , '*','Color',[0.5,0.25,0.25],'MarkerSize',20,'LineWidth',2);
% plot(kron(tx_N_,[1 1 1]), kron(null,[1 0 nan]) + kron(sx_N_,[0 1 nan]),'--','Color',[0.5,0.25,0.25],'LineWidth',2);
 
plot(t0_,grid,'--k');

legend({'Derivative','note_2: slope variation'},'FontSize',15,'Location','NW');  
%legend({'Actual signal','Sampled signal','note_x = 0.1 Note_1 + 0.1 Note_2 + 0.8 Note_3'},'FontSize',15,'Location','NW');  
set(gca,'xtick',[])


%% Plot clustering objective
null = zeros(1,length(tx_frame));

note_r = note_x_frame(note_x_frame>1);
note_b =  note_x_frame(note_x_frame<1);

tx_r = tx_frame(note_x_frame>1);
tx_b = tx_frame(note_x_frame<1);

null_r = zeros(1,length(note_r));
null_b = zeros(1,length(note_b));

%% NOTES CLUSTERED
figure(1);
plot(tx_r,note_r,'pr','MarkerSize',25,'LineWidth',2);
hold on
plot(tx_b,note_b,'pb','MarkerSize',25,'LineWidth',2);
plot(kron(tx_r,[1 1 1]), kron(null_r,[1 0 nan]) + kron(note_r,[0 1 nan]),'-k','LineWidth',2);
plot(kron(tx_b,[1 1 1]), kron(null_b,[1 0 nan]) + kron(note_b,[0 1 nan]),'-k','LineWidth',2);
xlabel ('time, s','FontSize',30,'FontWeight','Bold');
ylabel ('note_x','FontSize',30,'FontWeight','Bold');
set(gca,'FontSize',15);
%% NOTES NOT CLUSTERED
figure(2);
plot(tx_r,note_r,'pK','MarkerSize',25,'LineWidth',2);
hold on
plot(tx_b,note_b,'pk','MarkerSize',25,'LineWidth',2);
plot(kron(tx_r,[1 1 1]), kron(null_r,[1 0 nan]) + kron(note_r,[0 1 nan]),'-k','LineWidth',2);
plot(kron(tx_b,[1 1 1]), kron(null_b,[1 0 nan]) + kron(note_b,[0 1 nan]),'-k','LineWidth',2);
xlabel ('time, s','FontSize',30,'FontWeight','Bold');
ylabel ('note_x','FontSize',30,'FontWeight','Bold');
set(gca,'FontSize',15);



