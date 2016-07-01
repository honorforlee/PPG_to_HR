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

%   - Divide timeline -
for k = 0 : length(s) / length(range) - 1
    
    t_div (k+1,:) =  t(  k*(length(range)) + 1 : (k+1)*(length(range)) ) ;
    s_div(k+1,:)= s(  k*(length(range)) + 1 : (k+1)*(length(range)) ) ;
    
end

%   - Peaks identification - 
[kx,tx,sx, dhi,dlo, td,d, tx_N,sx_N, note_x] = signal_peaks(t_div(1,:),s_div(1,:));

%   - noteering according to minimum variance of note_x -
eps = 0.1;      
kx_ = kx;

for i = 1:length(kx)
    if kx_(i) ~= 0

        idx(1,i) = kx(i);           
        note(1,i) = note_x(i);     
        per(1,i) = tx(i);
        j = 2;
        
        for k = i + 1 : length(kx)
            
            if  var([note_x(i) note_x(k)],1) < eps;
                
                note(j,i) = note_x(k);
                per(j,i) = tx(k);             
                idx(j,i) = kx(k);
                kx_(k) = 0;
                
                j = j+1;
            else
                note(j,i)=nan;
                per(j,i)=nan; 
                idx(j,i)=nan;
                
                j = j+1;
            end
            
        end
    end
end

%   - Remove columns -
zero = find(~idx(1,:));

for k = 1:length(zero)
   note(:,zero(k))=[];
   per(:,zero(k))=[];
   idx(:,zero(k))=[];
   zero = bsxfun(@minus ,zero,ones(1,length(zero))) ;
end

%   - Cluster notation: size, mean note, periodicity -


L = size(idx);
for k = 1:L(2)
   clust_size(k) = nnz(idx(:,k)) - sum(isnan(idx(:,k)));
       if clust_size(k) >= 2
            NAN_idx = ~isnan(per(:,k));         % extract non NAN value of per, note
            NAN_per = per(NAN_idx,k);
            NAN_note = note(NAN_idx,k);
            
            per_ = NAN_per(find(NAN_per));      % extract non zero value of per, note
            note_ = NAN_note(find(NAN_note));
                        
        [PER_T(k),PER_eps(k), PER_R(k)] = periodicity(per_);    % note
        NOTE(k) = mean(note_);    
        SIZE(k) = clust_size(k);
        
        clust_note(k) = (0.5 * NOTE(k) + 0.5 * SIZE(k)) / PER_eps(k);
        
        clear NAN_idx NAN_per NAN_note per_ note_
    end
end

tbl_note = table([1:L(2)]', PER_T',PER_eps', PER_R', NOTE', SIZE', clust_note','VariableNames',{'Cluster','T','eps','R','Note_x','Size','Cluster_note'})



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

% plot note_note_x
for i = 1 : kmax
    figure(2);
    subplot(2,1,1);
    plot(note_note_x{i,kmax} , '.');
    hold on
end
Legend=cell(kmax,1);
for iter=1:kmax
    Legend{iter}=strcat('noteer ', num2str(iter));
end
legend(Legend);

hold off
subplot(2,1,2);
plot( note_x, '.');

base_array = cellfun(@length,note_tx);
base_max = max (base_array(:,kmax));
base = [1:base_max];

% plot note_periodicity
data = nan(base_max,kmax);

for i = 1 : kmax

    data(1:base_array(i,kmax),i) = note_tx{i,kmax};
    figure(3);
    tx_disp(i) = plot(base,data(:,i),'.');
    hold on
end

Legend=cell(kmax,1);
for iter=1:kmax
    Legend{iter}=strcat('noteer ', num2str(iter),': T = ', num2str(note_periodicity{iter,kmax}(1)), '; eps = ', num2str(note_periodicity{iter,kmax}(2)), '; R = ', num2str(note_periodicity{iter,kmax}(3)));
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
%         note_per(k-1:k+1) = note_x(k-1:k+1);
%     
%         
%     end    
% else
%     periodic_note_x(1) = note_x(1);
%     
%     
% end
