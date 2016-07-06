%Name = '3900679m';       % row 5
%Name = '3801060_0007m';  % row 1
%Name = '3914288m';       % row 5
%Name = '3919370m';       % row 5
%Name = '3916979m (1)';     % row 6  
Name = '3916979m (3)';     % row 6  
load(strcat(Name, '.mat'));
fid = fopen(strcat(Name, '.info'), 'rt');
fgetl(fid);
fgetl(fid);
fgetl(fid);
[interval] = sscanf(fgetl(fid), 'Sampling frequency: %f Hz  Sampling interval: %f sec');
dt0 = interval(2);              % data acquisition rate (interval = 1/f_spl_u = 0.5903 ms in practice)

fclose(fid);

t0 = (1:length(val)) * dt0;            % timeline
s0 = val(6,1:length(val));
s0  = (s0  - mean(s0 ))/sqrt(var(s0 ));        % rescale s on 0 (standard score of signal)

%   - Timeline, noise, integration, quantization -
dt = 1/20;                           % sampling time: dt >> dt0
t_int = dt * (1/3);                  % integration time: dt0 <= t_int < dt
quant = 0.1;                         % LSB: vertical step

[t,s] = integration(t0,s0,dt0,dt,t_int,quant,0);

[t0_ s0_ t_ s_] = time_div(t0,s0,dt0, t,s,dt,10,3);


%  - Peaks identification -
[kx,tx,sx, dhi,dlo, td,d, tx_N,sx_N, note_x] = signal_peaks(t_,s_); 

%  - Local maxima sx, maximum slope around sx -
[tx,sx, dhi,dlo, tx_N,sx_N, note_x] = peaks_processing(t_,s_,kx);


% d = s_(2:end) -  s_(1:end-1);
% td = t_(2:end);
% kx = d > 0;
% kx = find(kx(1:end-1) & ~kx(2:end));       % k_{x}:index where d > 0; d( k_{x} + 1 ) <= 0
% 
% sx = s_(kx+1);                          % local maxima
% tx = td(kx) + (td(kx+1)-td(kx)) .* d(kx)./(d(kx)-d(kx+1));      % linear interpolation of dhi and dho to get tx (@zero crossing)
% 
% dhi = d(kx);
% dlo = d(kx+1);
% 
% for k = 1:length(kx)
%     i = kx(k)-1;   while i > 0             && d(i) >= dhi(k); dhi(k) = d(i); i = i-1; end    % search for maximum positive slope at kx-
%     i = kx(k)+2;   while i < length(d) && d(i) <= dlo(k); dlo(k) = d(i); i = i+1; end    % search for maximmum negative slope at kx+
% end
% 
% kx_n = d < 0;                               % search for local minima
% kx_n = find(kx_n(1:end-1) & ~kx_n(2:end));
% 
% if kx_n(1) < kx(1)
%     for k = 1:length(kx)
%         kx_index(k) = max( find( kx_n < kx(k) ) );
%     end
%     sx_N = s_(kx_n( kx_index ) + 1);
%     tx_N = td(kx_n( kx_index )) + (td(kx_n( kx_index )+1)-td(kx_n( kx_index ))) .* d(kx_n( kx_index ))./(d(kx_n( kx_index ))-d(kx_n( kx_index )+1));
% else
%     kx_index(1) = nan;
%     sx_N(1) = nan;
%     tx_N(1) = nan;
%     
%     for k = 2:length(kx)
%         kx_index(k) = max( find( kx_n < kx(k) ) );
%         sx_N(k) = s_(kx_n( kx_index(k) ) + 1);
%         tx_N(k) = td(kx_n( kx_index(k) )) + (td(kx_n( kx_index(k) )+1)-td(kx_n( kx_index(k) ))) .* d(kx_n( kx_index(k) ))./(d(kx_n( kx_index(k) ))-d(kx_n( kx_index(k) )+1));
%     end
%     
% end
% 
% %   - Peaks notation
% note_1 = sx;
% for k = 2:length(kx)-1
%     note_1(k) = 2*sx(k) - sx(k+1) - sx(k-1);  % average peak value (doubled)
% end
% 
% note_2 = dhi - dlo;                           % maximum slope difference around peak
% 
% for k = 1:length(tx)
%     if tx(k) >= tx_N(k)
%         delta(k) = sx(k) - sx_N(k);
%     else
%         delta(k)=sx(k) - sx_N(k+1);
%     end
% end
% 
% note_x = 0.2*note_1 + 0.2*note_2 + 0.6*delta;

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
    NAN_ = ~isnan(idx(:,k));         % extract non NAN value of idx, per, note
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
        
        if PER_eps(k) <= 0.1            % best periodicity note set to 0.1 otherwise increase too much the cluster note 
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

if all(clust_note == 0)      % all cluster have zero note (all size <= 2)
[NOTE_major major_idx] = max(NOTE);
kx_major = clust_cell{major_idx,1}'; 
else
[clust_note_max major_idx] = max(clust_note_pos);
kx_major = clust_cell{major_idx,1}';
NOTE_major = NOTE(major_idx);
end

%   - Merge sub-major clusters -
Nrows = max(cellfun(@numel,clust_cell));
X = nan(Nrows(1),L(2));
for iCol = 1:L(2)   
    X(1:numel(clust_cell{iCol}),iCol) = clust_cell{iCol};       % copy idx values of each cluster into X
end
clust_merge = nan(Nrows(1),L(2));
clust_merge(:,major_idx) = X(:,major_idx);

for k = 1:L(2)
    if NOTE(k) > 0
        if var([NOTE_major NOTE(k)],1) < 3*eps && k ~= major_idx   % EMPIRICAL
            clust_merge(:,k) = X(:,k);
            NOTE_major = mean( [NOTE_major NOTE(k)] );              % adjust NOTE_major ti merged cluster
        end
    end     
end

clust_merge(isnan(clust_merge)) = [];         % remove NaN values
clust_merge = unique(clust_merge);            % remove repeated elements ans sort array

%   - Major peaks -
kx_major(1,1:length(clust_merge)) = clust_merge; 

tx_major = td(kx_major) + (td(kx_major+1)-td(kx_major)) .* d(kx_major)./(d(kx_major)-d(kx_major+1));      % linear interpolation of dhi and dho to get tx (@zero crossing)
sx_major = s_(kx_major+1);          % local maxima
T = mean(delta_tx(tx_major));

%   - Modify cluster according to periodicity -
% ADD PEAK TO MAJOR CLUSTER
tx_pos = delta_tx(tx_major);
kx_add = nan(1,length(kx_major));

for k = 1:length(tx_pos)                % assume ONE missing/skipped peak
    if tx_pos(k) > T + 0.5*T            % need enough large frame length to give weight to T
        left(k) = kx_major(k);
        right(k) = kx_major(k+1);
        kx_add_ = kx( kx(1,:) > left(k) & kx(1,:) < right(k));
        
        if length(kx_add_) == 1
            kx_add(k) = kx_add_;
            
        elseif length(kx_add_) >= 2
            for i = 1:length(kx_add_)
                kx_add_idx(i) = find(kx == kx_add_(i));
                kx_add_note(i) = note_x(kx_add_idx(i));
            end
            
            [value kx_add_max] = max(kx_add_note);
            kx_add(k) = kx_add_(kx_add_max);        % max note_x index added to major cluster
            
            clearvars kx_add_idx kx_add_note kx_add_max value ;
            
        else
            kx_add_ = nan;
            tx_pos(k) = nan;            % to compute T not affected by missing tx_major
            kx_add(k) = 0;
            
        end
        
    end
    clearvars kx_add_;
end
kx_add_temp = kx_add;
kx_add_temp(isnan(kx_add_temp)) = [];   % remove NaN to find zeros

if find(~kx_add_temp)                   % one imaginary peak to create
    zeros = find(~kx_add);              % indexes to create peak
    T_temp = mean(tx_pos,'omitnan');    % peaks period omitting peak missing
    
    insert = @(a, x, n)cat(2,  x(1:n), a, x(n+1:end));      % insert: function to insert value in array
    
    for k = 1:length(zeros)
        tx_major = insert(tx_major(zeros(k)) + T_temp, tx_major, zeros(k)); % insert tx_major
        sx_major = insert(sx_major(zeros(k)), sx_major, zeros(k) );         % inset sx_major at added time 
        zeros = bsxfun(@plus ,zeros,ones(1,length(zeros)));                 % shift index of zeros when adding one element in tx_major, sx_major
    end
   
else
    kx_major = horzcat(kx_major,kx_add);        % add peak to major cluster
    kx_major(isnan(kx_major)) = [];             % remove NaN values
    kx_major = unique(kx_major);                % sort
    
    tx_major = td(kx_major) + (td(kx_major+1)-td(kx_major)) .* d(kx_major)./(d(kx_major)-d(kx_major+1));      % linear interpolation of dhi and dho to get tx (@zero crossing)
    sx_major = s_(kx_major+1);                  % local maxima
end

T = mean(delta_tx(tx_major));

% REMOVE PEAK FROM MAJOR CLUSTER
tx_neg = delta_tx(tx_major);
for k = 1:length(tx_neg)
    if tx_neg(k) < T - T*0.5        % remove peaks - another loop because of matrix size inconsistency
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
