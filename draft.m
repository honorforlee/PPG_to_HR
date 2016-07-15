%Name = '3900497mB';        % row 6
%Name = '3900679m';         % row 5
%Name = '3914288m';         % row 5
%Name = '3916979m (4)';     % row 6  (1 : 5)
%Name = '3916979m';         % row 6
Name = '3919370m (1)';     % row 5
%Name = '3919370m';         % row 5
%Name = '3801060_0007m';    % row 1
%Name = '3899985_0005m';    % row 1
%Name = 'a02m';             % row 1

load(strcat(Name, '.mat'));
fid = fopen(strcat(Name, '.info'), 'rt');
fgetl(fid);
fgetl(fid);
fgetl(fid);
[interval] = sscanf(fgetl(fid), 'Sampling frequency: %f Hz  Sampling interval: %f sec');
dt0 = interval(2);              % data acquisition rate (interval = 1/f_spl_u = 0.5903 ms in practice)

fclose(fid);

val(isnan(val)) = [];
t0 = (1:length(val)) * dt0;            % timeline
s0 = val(5,1:length(val));
s0  = (s0  - mean(s0 ))/sqrt(var(s0));        % rescale s on 0 (standard score of signal)

%   - Timeline, noise, integration, quantization -
dt = 1/20;                           % sampling time: dt >> dt0
t_int = dt * (1/3);                  % integration time: dt0 <= t_int < dt
quant = 0.1;                         % LSB: vertical step

[t,s] = integration(t0,s0,dt0,dt,t_int,quant,0);

[t0_, s0_, t_, s_] = time_div(t0,s0,dt0, t,s,dt,5,7);

%  - Peaks identification -
%[kx,tx,sx, dhi,dlo, td,d, kx_n,tx_N,sx_N, note_x] = signal_peaks(t_,s_);

d = s_(2:end) -  s_(1:end-1);
td = t_(2:end);
kx = d > 0;
kx = find(kx(1:end-1) & ~kx(2:end));       % k_{x}:index where d(k(x) > 0; d( k(x) + 1 ) <= 0

sx = s_(kx+1);                          % local maxima
tx = td(kx) + (td(kx+1)-td(kx)) .* d(kx)./(d(kx)-d(kx+1));      % linear interpolation of dhi and dho to get tx (@zero crossing)

dhi = d(kx);
dlo = d(kx+1);

for k = 1:length(kx)
    i = kx(k)-1;   while i > 0             && d(i) >= dhi(k); dhi(k) = d(i); i = i-1; end    % search for maximum positive slope at kx-
    i = kx(k)+2;   while i < length(d) && d(i) <= dlo(k); dlo(k) = d(i); i = i+1; end    % search for maximmum negative slope at kx+
end

kx_n = d < 0;                               % search for local minima
kx_n = find(kx_n(1:end-1) & ~kx_n(2:end));

if kx_n(1) < kx(1)
    for k = 1:length(kx)
        kx_index(k) = max( find( kx_n < kx(k) ) );
    end
    sx_N = s_(kx_n( kx_index ) + 1);
    tx_N = td(kx_n( kx_index )) + (td(kx_n( kx_index )+1)-td(kx_n( kx_index ))) .* d(kx_n( kx_index ))./(d(kx_n( kx_index ))-d(kx_n( kx_index )+1));
else
    j = 1;
    while kx_n(1) > kx(j)                   % find first minima preceding maxima
        kx_index(j) = nan;
        sx_N(j) = nan;
        tx_N(j) = nan;
        j = j+1;
    end
    for k = j:length(kx)
        kx_index(k) = max( find( kx_n < kx(k) ) );
        sx_N(k) = s_(kx_n( kx_index(k) ) + 1);
        tx_N(k) = td(kx_n( kx_index(k) )) + (td(kx_n( kx_index(k) )+1)-td(kx_n( kx_index(k) ))) .* d(kx_n( kx_index(k) ))./(d(kx_n( kx_index(k) ))-d(kx_n( kx_index(k) )+1));
    end
    
end

%   - Peaks notation -
note_1 = sx;
note_2 = dhi - dlo;                           % maximum slope difference around peak

for k = 1:length(tx)
    if tx(k) >= tx_N(k)
       delta(k) = sx(k) - sx_N(k);       
    else                                     % if minimum out of frame, take first min in the frame
        j = k;
        while isnan(sx_N(j))
               j = j+1;
        end
        
        delta(k) = sx(k) - sx_N(j);
        clearvars j;
    end
end
note_3 = delta;

for k = 2:length(kx)-1
    note_1(k) = sx(k) - ( sx(k+1) + sx(k-1) )/2;                                % average peak value
    note_3(k) = delta(k) - ( delta(k+1) + delta(k-1) )/2;                       % average peak to peak value
end

note_x = 0.1*note_1 + 0.1*note_2 + 0.8*delta;

%       -- Minimum variance algorithm --
% [kx_major,tx_major,sx_major, T] = min_variance(t_,s_, td,d, kx,tx,sx,note_x, 0.1);

eps = 0.1;
% Ivan NY HANITRA - Master thesis
%       -- Clustering according to minimum variance of note_x  --

eps = 0.1;
kx_ = kx;

%   - Clustering according to minimum variance of note_x -
for i = 1:length(kx)
    if kx_(i) ~= 0
        
        idx(1,i) = kx(i);
        per(1,i) = tx(i);
        note(1,i) = note_x(i);
        j = 2;
        
        for k = i + 1 : length(kx)
            
            if  similarity(note_x(i),note_x(k),'variance') < eps;
                
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

%   - Cluster notation: size, tx periodicity, note_x -
for k = 1:L(2)
    if length(clust_cell{k,1}) > 2
        
        SIZE(k) = length(clust_cell{k,1});                                 % size
        [PER_T(k),PER_eps(k), PER_R(k)] = periodicity(clust_cell{k,2});    % tx periodicity
        
        if PER_eps(k) <= 0.1            % best periodicity note set to 0.1 otherwise increase too much the cluster note
            PER_eps(k) = 0.1;
        end
        
        NOTE(k) = mean(clust_cell{k,3});                                             % average note_x
        clust_note(k) = (0.7 * NOTE(k) + 0.2 * SIZE(k)) / (PER_eps(k)/0.1);          % cluster note EMPIRICAL
        
    else
        SIZE(k) = length(clust_cell{k,1});
        PER_T(k) = 0; PER_eps(k) = 0; PER_R(k) = 0;
        NOTE(k) = mean(clust_cell{k,3});
        clust_note(k) = 0;
        
    end
end

tbl_note = table([1:L(2)]', SIZE', PER_T',PER_eps', PER_R', NOTE', clust_note','VariableNames',{'Cluster','Size','T','eps','R','Note_x','Cluster_note'});

%   - Select major cluster + merge sub-major clusters -
Nrows = max(cellfun(@numel,clust_cell));
X = nan(Nrows(1),L(2));
for iCol = 1:L(2)
    X(1:numel(clust_cell{iCol}),iCol) = clust_cell{iCol};       % copy idx values of each cluster into X
end
clust_merge = nan(Nrows(1),L(2));

if L(2) >= 2                                             % more than 1 cluster
    if all(clust_note == 0)
        [NOTE_major major_idx] = max(NOTE);
        kx_major = clust_cell{major_idx,1}';
                
        if all(NOTE <= 1)                                                                         % EMP
            for k = 1:L(2)
                if similarity(NOTE_major, NOTE(k),'variance') < 7*eps                             % EMPIRICAL: compare NOTE to NOTE of major cluster
                    NOTE_major = NOTE(k) * ( Nrows(1) - sum(isnan(X(:,k))) ) + NOTE_major * ( L(2)*Nrows(1) - sum(sum(isnan(clust_merge))) );
                    clust_merge(:,k) = X(:,k);                                                    % merge cluster to major cluster
                    NOTE_major = NOTE_major / ( L(2)*Nrows(1) - sum(sum(isnan(clust_merge))) );   % recompute NOTE_major
                end
            end
        else
            for k = 1:L(2)
                if similarity(NOTE_major, NOTE(k),'variance') < 7*eps && NOTE(k) > 1                % EMPIRICAL: compare NOTE to NOTE of major cluster
                    NOTE_major = NOTE(k) * ( Nrows(1) - sum(isnan(X(:,k))) ) + NOTE_major * ( L(2)*Nrows(1) - sum(sum(isnan(clust_merge))) );
                    clust_merge(:,k) = X(:,k);                                                      % merge cluster to major cluster
                    NOTE_major = NOTE_major / ( L(2)*Nrows(1) - sum(sum(isnan(clust_merge))) );     % recompute NOTE_major                  
                end
            end
        end
        
    else
        if all(NOTE <= 1)                                                                            % EMP
            [NOTE_major major_idx] = max(NOTE);
            kx_major = clust_cell{major_idx,1}';
            
            for k = 1:L(2)
                if similarity(NOTE_major, NOTE(k),'variance') < 7*eps                                % EMPIRICAL: compare NOTE to NOTE of major cluster
                    NOTE_major = NOTE(k) * ( Nrows(1) - sum(isnan(X(:,k))) ) + NOTE_major * ( L(2)*Nrows(1) - sum(sum(isnan(clust_merge))) );
                    clust_merge(:,k) = X(:,k);                                                       % merge cluster to major cluster
                    NOTE_major = NOTE_major / ( L(2)*Nrows(1) - sum(sum(isnan(clust_merge))) );      % recompute NOTE_major
                end
            end
            
        else
            [clust_note_temp idx_temp] = max(clust_note);
            
            if NOTE(idx_temp) > 1                                                                  % EMP
                major_idx = idx_temp;
                kx_major = clust_cell{major_idx,1}';
                clust_note_major = clust_note_temp;
                NOTE_major = NOTE(major_idx);
                
                for k = 1:L(2)
                    if similarity(NOTE_major, NOTE(k),'variance') < 7*eps && clust_note(k) > 1 && NOTE(k) > 1           % EMPIRICAL: compare NOTE to NOTE of major cluster
                        NOTE_major = NOTE(k) * ( Nrows(1) - sum(isnan(X(:,k))) ) + NOTE_major * ( L(2)*Nrows(1) - sum(sum(isnan(clust_merge))) );
                        clust_merge(:,k) = X(:,k);                                                                      % merge cluster to major cluster
                        NOTE_major = NOTE_major / ( L(2)*Nrows(1) - sum(sum(isnan(clust_merge))) );                     % recompute NOTE_major

                    elseif similarity(NOTE_major, NOTE(k),'variance') < 7*eps && clust_note(k) == 0 && NOTE(k) > 1      % EMPIRICAL
                        NOTE_major = NOTE(k) * ( Nrows(1) - sum(isnan(X(:,k))) ) + NOTE_major * ( L(2)*Nrows(1) - sum(sum(isnan(clust_merge))) );
                        clust_merge(:,k) = X(:,k);                                                                      % merge cluster to major cluster
                        NOTE_major = NOTE_major / ( L(2)*Nrows(1) - sum(sum(isnan(clust_merge))) );                     % recompute NOTE_major
                    end
                end
                
            else
                [NOTE_major major_idx] = max(NOTE);
                kx_major = clust_cell{major_idx,1}';
                
                for k = 1:L(2)
                    if similarity(NOTE_major, NOTE(k),'variance') < 7*eps && NOTE(k) > 1                             % EMPIRICAL: compare NOTE to NOTE of major cluster
                        NOTE_major = NOTE(k) * ( Nrows(1) - sum(isnan(X(:,k))) ) + NOTE_major * ( L(2)*Nrows(1) - sum(sum(isnan(clust_merge))) );
                        clust_merge(:,k) = X(:,k);                                                      % merge cluster to major cluster
                        NOTE_major = NOTE_major / ( L(2)*Nrows(1) - sum(sum(isnan(clust_merge))) );     % recompute NOTE_major
                    end
                end
            end
        end
    end
    
    clust_merge(isnan(clust_merge)) = [];         % remove NaN values
    clust_merge = unique(clust_merge);            % remove repeated elements ans sort array
    kx_major(1,1:length(clust_merge)) = clust_merge;
    
else                                              % one cluster only
    major_idx = 1;
    kx_major = clust_cell{major_idx,1}';
    clust_note_major = clust_note;
    NOTE_major = NOTE;
end

if length(kx_major) >= 2
    
    %   - Major peaks -
    kx_major = unique(kx_major);
    tx_major = td(kx_major) + (td(kx_major+1)-td(kx_major)) .* d(kx_major)./(d(kx_major)-d(kx_major+1));      % linear interpolation of dhi and dho to get tx (@zero crossing)
    sx_major = s_(kx_major+1);          % local maxima
    T = mean(delta_tx(tx_major));
    
    % - Remove some merged peaks -
    [kx_major,tx_major,sx_major,T] = remove_peaks(kx_major,tx_major,sx_major, T, kx,note_x);
    
    %   - Add peaks to major cluster considering peaks periodicity -
    % Periodic peaks in row
    tx_rect = delta_tx(tx);
    
    kx_major_ = nan(1,length(kx)+length(kx_major));
    kx_major_(1:length(kx_major)) = kx_major;        % length(kx_major) < length(kx)
    
    for k = 1:length(tx_rect)
        if similarity(T,tx_rect(k), 'relative') < 0.5             % less than 50% relative error from T
            if similarity(note_x(k), note_x(k+1), 'relative') < 0.5 && (  ~any(kx_major == kx(k)) || ~any(kx_major == kx(k+1)) )           % similar note_x and kx not present in kx_major
                kx_major_(length(kx_major)+k) = kx(k);
                kx_major_(length(kx_major)+k+1) = kx(k+1);
                
                %             elseif abs( tx_rect(k) - T )/T < 0.2  && (  ~any(kx_major == kx(k)) || ~any(kx_major == kx(k+1)) )
                %                 kx_major_(length(kx_major)+k) = kx(k);
                %                 kx_major_(length(kx_major)+k+1) = kx(k+1);
            end
        end
    end
    
    kx_major_(isnan(kx_major_))=[];
    kx_major = unique(kx_major_);
    
    tx_major = td(kx_major) + (td(kx_major+1)-td(kx_major)) .* d(kx_major)./(d(kx_major)-d(kx_major+1));      % linear interpolation of dhi and dho to get tx (@zero crossing)
    sx_major = s_(kx_major+1);          % local maxima
    T = mean(delta_tx(tx_major));
    
    % Periodic peaks separated by minor peak
    tx_rect2 = delta_tx(tx,2);
    
    kx_major_ = nan(1,length(kx)+length(kx_major));
    kx_major_(1:length(kx_major)) = kx_major;        % length(kx_major) < length(kx)
    
    for k = 1:length(tx_rect2)
        if similarity(T, tx_rect2(k), 'relative') < 0.5            % less than 50% relative error from T and avoid periodic minor peaks
            if similarity(note_x(k), note_x(k+2), 'relative') < 0.5 && similarity(NOTE_major, note_x(k), 'relative') < 0.2 && ( ~any(kx_major == kx(k)) || ~any(kx_major == kx(k+2)) )
                kx_major_(length(kx_major)+k) = kx(k);
                kx_major_(length(kx_major)+k+1) = kx(k+2);
            end
        end
    end
    
    kx_major_(isnan(kx_major_))=[];
    kx_major = unique(kx_major_);
    
    tx_major = td(kx_major) + (td(kx_major+1)-td(kx_major)) .* d(kx_major)./(d(kx_major)-d(kx_major+1));      % linear interpolation of dhi and dho to get tx (@zero crossing)
    sx_major = s_(kx_major+1);          % local maxima
    T = mean(delta_tx(tx_major));
    
    %       - Search for missing peaks -
    % Miss peaks inside the frame: add/create 2 peaks max in a hole
    loop = 0;
    
    while loop < 2
        tx_pos = delta_tx(tx_major);
        
        % Search for missing peaks inside the frame
        [kx_add,tx_pos] = missing_peaks(kx,tx, kx_major,tx_major, tx_pos,T, note_x,NOTE_major);
        
        % Add/create peak to major cluster
        [kx_major, tx_major, sx_major, T] = add_peaks(t_,s_,td,d, kx_major,tx_major,sx_major, kx_add,tx_pos);
        
        clearvars tx_pos kx_add;
        loop = loop+1;
    end
    
    % Miss first peak
    insert = @(a, x, n)cat(2,  x(1:n), a, x(n+1:end));      % insert(element inserted,array,position)
    
    if similarity(T,abs( tx(1) - tx_major(1) ), 'relative') < 0.2 || abs( tx(1) - tx_major(1) ) > T + T/5             % 20% relative error or more than 20% error above T
        kx_major = insert(kx(1),kx_major,0);
        tx_major = td(kx_major) + (td(kx_major+1)-td(kx_major)) .* d(kx_major)./(d(kx_major)-d(kx_major+1));          % linear interpolation of dhi and dho to get tx (@zero crossing)
        sx_major = insert(sx_major(1),sx_major,0);                                                                    % for visibility
        T = mean(delta_tx(tx_major));
    end
    
    % Miss last peak
    if similarity(T,abs( tx(end) - tx_major(end) ), 'relative') < 0.2  || abs( tx(end) - tx_major(end) ) > T + T/5    % 20% relative error or more than 20% error above T
        kx_major = insert(kx(end),kx_major,length(kx_major));
        tx_major = td(kx_major) + (td(kx_major+1)-td(kx_major)) .* d(kx_major)./(d(kx_major)-d(kx_major+1));          % linear interpolation of dhi and dho to get tx (@zero crossing)
        sx_major = insert(sx_major(end),sx_major,length(sx_major));                                                   % for visibility
        T = mean(delta_tx(tx_major));
    end
    
    %   - Remove peaks from major cluster -
    [kx_major,tx_major,sx_major,T] = remove_peaks(kx_major,tx_major,sx_major, T, kx, note_x);
    
else
    display('Not enough points');
    tx_major = nan;
    sx_major = nan;
    T = nan;
    
end
%%
%   - Plots -
figure(2);
plot( tx , sx   , 'dc','MarkerSize',12);
hold on
plot(t_,s_,'ok','LineWidth',.2);
plot(t0_,s0_,'-k');
plot(tx, dhi,'^b','MarkerSize',10);
plot(tx, dlo,'vb','MarkerSize',10);
plot( kron(tx,[1 1 1]) , kron(dlo,[1 0 nan]) + kron(dhi,[0 1 nan]), '-b');       % link note_2
plot( tx_major , sx_major, 'pr','MarkerSize',20);
plot( tx,sx_N, 'dc','MarkerSize',10);
plot(kron(tx,[1 1 1]), kron(sx_N,[1 0 nan]) + kron(sx,[0 1 nan]),'-c');
hold off
