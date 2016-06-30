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

range = (1 : (10/dt)) * dt;

for k = 0 : length(s) / length(range) - 1
    
    t_div (k+1,:) =  t(  k*(length(range)) + 1 : (k+1)*(length(range)) ) ;
    s_div(k+1,:)= s(  k*(length(range)) + 1 : (k+1)*(length(range)) ) ;
    
end

[kx,tx,sx, dhi,dlo, td,d, tx_N,sx_N, note_x] = signal_peaks(t_div(1,:),s_div(1,:));

eps = 0.1;

for k = 1 : length(kx)-1
    Var_1(k) = var( [note_x(k+1) note_x(k)] );
end

for k = 1 : length(kx)-2
    Var_2(k) = var( [note_x(k+2) note_x(k)] );
end

% for k = 2:length(kx)  
%     if var([ref(ref_) note_x(k)],1) < eps
%         clust(k - ref_,ref_) = [note_x(k)];
%         idx(k - ref_,ref_) = [kx(k)];
%         ref(k) = ref(ref_);
%      elseif var([ref(ref_) note_x(k+1)],1) < eps
%         clust(k - ref_,ref_) = [note_x(k+1)];
%         idx(k - ref_,ref_) = [kx(k+1)];
%         ref(k) = ref(ref_);
%     else
%          ref(k) = note_x(k);
%          ref_ = k;        
%        end
% end   

kx_ = kx;

for i = 1:length(kx)
    if kx_(i) ~= 0

        idx(1,i) = kx(i);
        clust(1,i) = note_x(i);
        per(1,i) = tx(i);
        j = 2;
        
        for k = i + 1 : length(kx)
            
            if  var([note_x(i) note_x(k)],1) < eps;
                
                clust(j,i) = note_x(k);
                per(j,i) = tx(k);             
                idx(j,i) = kx(k);
                kx_(k) = 0;
                
                j = j+1;
            else
                clust(j,i)=nan;
                per(j,i)=nan; 
                idx(j,i)=nan;
                
                j = j+1;
            end
            
        end
    end
end

zero = find(~idx(1,:));

for k = 1:length(zero)
   clust(:,zero(k))=[];
   per(:,zero(k))=[];
   idx(:,zero(k))=[];
   zero = bsxfun(@minus ,zero,ones(1,length(zero))) ;
end

L = size(idx);
for k = 1:L(2)
   clust_size(k) = nnz(idx(:,k)) - sum(isnan(idx(:,k)));
%    clust_note(k) = mean(clust(i))*  
end

%%
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
%     periodic_note_x(1) = note_x(1);
%     
%     
% end
