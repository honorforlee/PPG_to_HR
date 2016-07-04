Name = '3900679m';      % row 5
%Name = '3801060_0007m';  % row 1
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
dt = 1/20;                           % sampling time: dt >> interval
t_int = dt * (1/3);                  % integration time: interval <= t_int < dt
quant = 0.1;                          % LSB: vertical step

[t,s] = integration(t0,s0,interval,dt,t_int,quant,0);

[t_ s_] = time_div(t,s,dt,10,1);

%   - Peaks identification -
[kx,tx,sx, dhi,dlo, td,d, tx_N,sx_N, note_x] = signal_peaks(t_,s_);

%   - Minimum variance algorithm -
% [kx_major,tx_major,sx_major, T] = min_variance(t_,s_, td,d, kx,tx,sx,note_x, 0.1);
eps = 0.1;
kx_ = kx;

for i = 1:length(kx)
    if kx_(i) ~= 0
        
        idx(1,i) = kx(i);
        per(1,i) = tx(i);
        note(1,i) = note_x(i);
        j = 2;
        
        for k = i + 1 : length(kx)
            
            if  var([note_x(i) note_x(k)],1) < eps;
                
                idx(j,i) = kx(k);
                per(j,i) = tx(k);
                note(j,i) = note_x(k);
                
                kx_(k) = 0;
                j = j+1;
            else
                idx(j,i)=nan;
                per(j,i)=nan;
                note(j,i)=nan;
                
                j = j+1;
            end
            
        end
    end
end

%   - Remove columns -
zero = find(~idx(1,:));

for k = 1:length(zero)
    idx(:,zero(k))=[];
    per(:,zero(k))=[];
    note(:,zero(k))=[];
    zero = bsxfun(@minus ,zero,ones(1,length(zero))) ;
end

%   - Create cells: kx-tx-note_x -
L = size(idx);
for k = 1:L(2)
    NAN_ = ~isnan(per(:,k));         % extract non NAN value of idx, per, note
    NAN_idx = idx(NAN_,k);
    NAN_per = per(NAN_,k);
    NAN_note = note(NAN_,k);
    
    idx_ = NAN_idx(find(NAN_idx));     % extract non zero value of idx, per, note
    per_ = NAN_per(find(NAN_per));
    note_ = NAN_note(find(NAN_note));
    
    clust_cell{k,1} = idx_;            % kx
    clust_cell{k,2} = per_;            % tx
    clust_cell{k,3} = note_;           % note_x
    
    clear NAN_ NAN_idx NAN_per NAN_note idx_ per_ note_
end

%   - Cluster notation: size, note_x, periodicity -
for k = 1:L(2)
    if length(clust_cell{k,1}) > 2
        
        SIZE(k) = length(clust_cell{k,1});
        [PER_T(k),PER_eps(k), PER_R(k)] = periodicity(clust_cell{k,2});    % note
        
        if PER_eps(k) <= 0.1
            PER_eps(k) = 0.1;
        end
                
        NOTE(k) = mean(clust_cell{k,3});       
        clust_note(k) = (0.4 * NOTE(k) + 0.3 * SIZE(k)) / (PER_eps(k)/0.3);
        
    else
        SIZE(k) = length(clust_cell{k,1});
        PER_T(k) = 0; PER_eps(k) = 0; PER_R(k) = 0;
        NOTE(k) = mean(clust_cell{k,3});
        clust_note(k) = 0;
        
    end
end

tbl_note = table([1:L(2)]', SIZE', PER_T',PER_eps', PER_R', NOTE', clust_note','VariableNames',{'Cluster','Size','T','eps','R','Note_x','Cluster_note'})

%   - Major cluster -
for k = 1:L(2)
    if NOTE(k) > 0
        clust_note_pos(k) = clust_note(k);
    end
end
[clust_note_max major_idx] = max(clust_note_pos);
kx_major = clust_cell{major_idx,1};

%   - Merge sub-major clusters -
Nrows = max(cellfun(@numel,clust_cell));
X = nan(Nrows(1),L(2));
for iCol = 1:L(2)
    X(1:numel(clust_cell{iCol}),iCol) = clust_cell{iCol};
end
clust_merge = nan(Nrows(1),L(2));
clust_merge(:,major_idx) = X(:,major_idx);

for k = 1:L(2)
    if NOTE(k) > 0
        if var([max(NOTE) NOTE(k)],1) < eps && k ~= major_idx
            clust_merge(:,k) = X(:,k);
        end
    end
end

clust_merge(isnan(clust_merge)) = [];   % remove NaN values

%   - Major peaks -
kx_major = unique(clust_merge);      % remove repeated elements
tx_major = td(kx_major) + (td(kx_major+1)-td(kx_major)) .* d(kx_major)./(d(kx_major)-d(kx_major+1));      % linear interpolation of dhi and dho to get tx (@zero crossing)
sx_major = s_(kx_major+1);          % local maxima
T = mean(delta_tx(tx_major));

%   - Modify cluster according to periodicity -
% ADD PEAK
tx_pos = delta_tx(tx_major);
kx_add = nan(1,length(kx_major));

for k = 1:length(tx_pos)
    if tx_pos(k) > T + 0.5*T            % add a major peak 
        left(k) = kx_major(k);
        right(k) = kx_major(k+1);
        
        kx_add(k) = kx( kx(1,:) > left(k) & kx(1,:) < right(k)); 
    end   
end

kx_major = horzcat(kx_major,kx_add);        % add peak to major cluster
kx_major(isnan(kx_major)) = [];             % remove NaN values
kx_major = unique(kx_major);                % sort

tx_major = td(kx_major) + (td(kx_major+1)-td(kx_major)) .* d(kx_major)./(d(kx_major)-d(kx_major+1));      % linear interpolation of dhi and dho to get tx (@zero crossing)
sx_major = s_(kx_major+1);          % local maxima
T = mean(delta_tx(tx_major));

% REMOVE PEAK
tx_neg = delta_tx(tx_major);
for k = 1:length(tx_neg)
    if tx_neg(k) < T - T*0.5        % remove peaks - another loop because of matrix size issue 
        kx_major(k+1) = nan;
    end
end

kx_major(isnan(kx_major)) = [];             % remove NaN values
kx_major = unique(kx_major);                % sort

tx_major = td(kx_major) + (td(kx_major+1)-td(kx_major)) .* d(kx_major)./(d(kx_major)-d(kx_major+1));      % linear interpolation of dhi and dho to get tx (@zero crossing)
sx_major = s_(kx_major+1);          % local maxima
T = mean(delta_tx(tx_major));



%   - Plots -
figure(2);
plot( tx , sx   , 'dr','MarkerSize',12);
hold on
plot(t_,s_,'-k','LineWidth',.2);
plot(tx, dhi,'c^','MarkerSize',12);
plot(tx, dlo,'cv','MarkerSize',12);
plot( kron(tx,[1 1 1]) , kron(dlo,[1 0 nan]) + kron(dhi,[0 1 nan]), '-c');       % link note_2
plot( tx_major , sx_major, 'pk','MarkerSize',15);
plot( tx,sx_N, 'dr','MarkerSize',12);
plot(kron(tx,[1 1 1]), kron(sx_N,[1 0 nan]) + kron(sx,[0 1 nan]),'r-');
hold off



%%
% %   - Cluster notation: size, mean note, periodicity -
% L = size(idx);
% for k = 1:L(2)
%     clust_size(k) = nnz(idx(:,k)) - sum(isnan(idx(:,k)));
%     if clust_size(k) > 2
%         NAN_idx = ~isnan(per(:,k));         % extract non NAN value of per, note
%         NAN_per = per(NAN_idx,k);
%         NAN_note = note(NAN_idx,k);
%         
%         per_ = NAN_per(find(NAN_per));      % extract non zero value of per, note
%         note_ = NAN_note(find(NAN_note));
%         
%         [PER_T(k),PER_eps(k), PER_R(k)] = periodicity(per_);    % note
%         
%         if PER_eps(k) <= 0.01
%             PER_eps(k) = 0.01;
%         end
%         NOTE(k) = mean(note_);
%         SIZE(k) = clust_size(k);
%         
%         clust_note(k) = (0.2 * NOTE(k) + 0.4 * SIZE(k)) / (PER_eps(k)/0.4);
%         
%         clear NAN_idx NAN_per NAN_note per_ note_
%     else
%         NAN_idx = ~isnan(per(:,k));
%         NAN_note = note(NAN_idx,k);
%         note_ = NAN_note(find(NAN_note));
%         
%         PER_T(k) = 0; PER_eps(k) = 0; PER_R(k) = 0; 
%         NOTE(k) = mean(note_);
%         SIZE(k) = clust_size(k);
%         clust_note(k) = 0;
%         
%         clear NAN_idx NAN_note note_
%     end
% end
% 
% tbl_note = table([1:L(2)]', PER_T',PER_eps', PER_R', NOTE', SIZE', clust_note','VariableNames',{'Cluster','T','eps','R','Note_x','Size','Cluster_note'})

note_major(1) = max(clust_note);

j=1;
for k = 1:L(2)
    
    if var([note_major(1) clust_note(k)],1) < 5*eps
        clust_major(j) = k;
        j = j+1;
    end
    
end

kx_major = kx(clust_major);


%%
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