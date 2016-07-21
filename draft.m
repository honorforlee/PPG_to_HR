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

frame_init = 50; frame_end = 55;

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
    
    clust_cell_temp{k,1} = idx_;            % kx
    clust_cell_temp{k,2} = per_;            % tx
    clust_cell_temp{k,3} = note_;           % note_x
    NOTE_mean(k) = mean(clust_cell_temp{k,3});
    
    clear NAN_ NAN_idx NAN_per NAN_note idx_ per_ note_
end

% Sort clusters by NOTE
NOTE = sort(NOTE_mean,'descend');                               % average of cluster note_x sorted

for k = 1:L(2)
    
    clust_cell{k,1} = clust_cell_temp{NOTE_mean == NOTE(k),1};      % kx sorted
    clust_cell{k,2} = clust_cell_temp{NOTE_mean == NOTE(k),2};      % tx sorted
    clust_cell{k,3} = clust_cell_temp{NOTE_mean == NOTE(k),3};      % note_x sorted
    
end

%   - Cluster notation: size, tx periodicity, note_x -
for k = 1:L(2)
    if length(clust_cell{k,1}) > 2
        
        SIZE(k) = length(clust_cell{k,1});                                 % size
        [PER_T(k),PER_eps(k), PER_R(k)] = periodicity(clust_cell{k,2});    % tx periodicity
        
        if PER_eps(k) <= 0.1            % best periodicity note set to 0.1 otherwise increase too much the cluster note
            PER_eps(k) = 0.1;
        end
        clust_note(k) = (0.7 * NOTE(k) + 0.2 * SIZE(k)) / (PER_eps(k)/0.1);          % cluster note EMPIRICAL
        
    else
        SIZE(k) = length(clust_cell{k,1});
        PER_T(k) = 0; PER_eps(k) = 0; PER_R(k) = 0;
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
        major_idx = 1;
        NOTE_major = NOTE(major_idx);
        kx_major = clust_cell{major_idx,1}';
        
        if all(NOTE <= 10*eps)                                                                         % EMP
            for k = 1:L(2)
                if similarity(NOTE_major, NOTE(k),'variance') < 7*eps                             % EMPIRICAL: compare NOTE to NOTE of major cluster
                    NOTE_major = NOTE(k) * ( Nrows(1) - sum(isnan(X(:,k))) ) + NOTE_major * ( L(2)*Nrows(1) - sum(sum(isnan(clust_merge))) );
                    clust_merge(:,k) = X(:,k);                                                    % merge cluster to major cluster
                    NOTE_major = NOTE_major / ( L(2)*Nrows(1) - sum(sum(isnan(clust_merge))) );   % recompute NOTE_major
                end
            end
        else
            for k = 1:L(2)
                if similarity(NOTE_major, NOTE(k),'variance') < 7*eps && NOTE(k) > 10*eps                % EMPIRICAL: compare NOTE to NOTE of major cluster
                    NOTE_major = NOTE(k) * ( Nrows(1) - sum(isnan(X(:,k))) ) + NOTE_major * ( L(2)*Nrows(1) - sum(sum(isnan(clust_merge))) );
                    clust_merge(:,k) = X(:,k);                                                      % merge cluster to major cluster
                    NOTE_major = NOTE_major / ( L(2)*Nrows(1) - sum(sum(isnan(clust_merge))) );     % recompute NOTE_major
                end
            end
        end
        
    else
        if all(NOTE <= 10*eps)                                                                            % EMP
            major_idx = 1;
            NOTE_major = NOTE(major_idx);
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
            
            if NOTE(idx_temp) > 10*eps                                                                  % EMP
                major_idx = idx_temp;
                NOTE_major = NOTE(major_idx);
                kx_major = clust_cell{major_idx,1}';
                
                for k = 1:L(2)
                    if similarity(NOTE_major, NOTE(k),'variance') < 7*eps && NOTE(k) > 10*eps           % EMPIRICAL: compare NOTE to NOTE of major cluster
                        NOTE_major = NOTE(k) * ( Nrows(1) - sum(isnan(X(:,k))) ) + NOTE_major * ( L(2)*Nrows(1) - sum(sum(isnan(clust_merge))) );
                        clust_merge(:,k) = X(:,k);                                                                      % merge cluster to major cluster
                        NOTE_major = NOTE_major / ( L(2)*Nrows(1) - sum(sum(isnan(clust_merge))) );                     % recompute NOTE_major
                    end
                end
            else
                clearvars clust_note_temp idx_temp major_idx;
                clust_note_temp_ = max(   clust_note(clust_note < max(clust_note))    );          % second max clust_note
                idx_temp_ = find(clust_note == clust_note_temp_);
                NOTE_temp_ = NOTE(idx_temp_);
                
                NOTE_temp = max(NOTE_temp_);
                idx_temp = find(NOTE == NOTE_temp);
                
                if NOTE(idx_temp) > 10*eps                                                                  % EMP
                    major_idx = idx_temp;
                    NOTE_major = NOTE(major_idx);
                    kx_major = clust_cell{major_idx,1}';
                    
                    for k = 1:L(2)
                        if similarity(NOTE_major, NOTE(k),'variance') < 7*eps  && NOTE(k) > 10*eps           % EMPIRICAL: compare NOTE to NOTE of major cluster
                            NOTE_major = NOTE(k) * ( Nrows(1) - sum(isnan(X(:,k))) ) + NOTE_major * ( L(2)*Nrows(1) - sum(sum(isnan(clust_merge))) );
                            clust_merge(:,k) = X(:,k);                                                                      % merge cluster to major cluster
                            NOTE_major = NOTE_major / ( L(2)*Nrows(1) - sum(sum(isnan(clust_merge))) );                     % recompute NOTE_major
                        end
                    end
                    
                else
                    major_idx = 1;
                    NOTE_major = NOTE(major_idx);
                    kx_major = clust_cell{major_idx,1}';
                    
                    for k = 1:L(2)
                        if similarity(NOTE_major, NOTE(k),'variance') < 7*eps && NOTE(k) > 10*eps                             % EMPIRICAL: compare NOTE to NOTE of major cluster
                            NOTE_major = NOTE(k) * ( Nrows(1) - sum(isnan(X(:,k))) ) + NOTE_major * ( L(2)*Nrows(1) - sum(sum(isnan(clust_merge))) );
                            clust_merge(:,k) = X(:,k);                                                      % merge cluster to major cluster
                            NOTE_major = NOTE_major / ( L(2)*Nrows(1) - sum(sum(isnan(clust_merge))) );     % recompute NOTE_major
                        end
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
    NOTE_major = NOTE;
end

if length(kx_major) >= 2
    warning = 0;
    
    %   - Major peaks -
    kx_major = unique(kx_major);
    idxs = arrayfun(@(x)find(kx==x,1),kx_major);
    tx_major = tx(idxs);
    sx_major = sx(idxs);
    
    T = delta_tx(tx_major);
    clearvars idxs;

%    % - Remove some merged peaks -
    %[kx_major,tx_major,sx_major,T,warning] = remove_peaks(kx_major,tx_major,sx_major, T, kx,note_x);
%     if warning == 1
%         display('No peaks detected')
%         return
%     end
    
    %   - Add peaks to major cluster considering peaks periodicity -
    % Periodic peaks in row
    tx_rect = delta_tx(tx);
    T_rect = mean(tx_rect);
    
    kx_major_ = nan(1,length(kx)+length(kx_major));
    kx_major_(1:length(kx_major)) = kx_major;        % length(kx_major) < length(kx)
    
    for k = 1:length(tx_rect)
        if similarity(T_rect,tx_rect(k), 'relative') < 0.2
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
    
    idxs = arrayfun(@(x)find(kx==x,1),kx_major);
    tx_major = tx(idxs);
    sx_major = sx(idxs);

    T = mean(delta_tx(tx_major));
    clearvars idxs;
    
    % Periodic peaks separated by minor peak
    tx_rect2 = delta_tx(tx,2);
    T_rect2 = mean(tx_rect2);
    
    kx_major_ = nan(1,length(kx)+length(kx_major));
    kx_major_(1:length(kx_major)) = kx_major;        % length(kx_major) < length(kx)
    
    for k = 1:length(tx_rect2)
        if similarity(T_rect2, tx_rect2(k), 'relative') < 0.2            % less than 50% relative error from T and avoid periodic minor peaks
            if similarity(note_x(k), note_x(k+2), 'relative') < 0.5  && sx(k) > sx(k+1) && sx(k+2) > sx(k+1) && ( ~any(kx_major == kx(k)) || ~any(kx_major == kx(k+2)) )
                kx_major_(length(kx_major)+k) = kx(k);
                kx_major_(length(kx_major)+k+1) = kx(k+2);
            end
        end
    end
    
    kx_major_(isnan(kx_major_))=[];
    kx_major = unique(kx_major_);
    
    idxs = arrayfun(@(x)find(kx==x,1),kx_major);
    tx_major = tx(idxs);
    sx_major = sx(idxs);

    T = mean(delta_tx(tx_major));
    clearvars idxs;
    
    [kx_major,tx_major,sx_major,T,warning] = remove_peaks(kx_major,tx_major,sx_major, T, kx, note_x);
    if warning == 1
        display('No peaks detected')
        return
    end
    
    
        
else
    warning = 1;
    tx_major = nan;
    sx_major = nan;
    T = nan;
    display('No peaks detected')
end



%%
%   - Plots -
figure(2);
plot( tx_frame , sx_frame   , 'dc','MarkerSize',12);
hold on
plot(t_frame,s_frame,'ok','LineWidth',.2);
plot(t0_frame,s0_frame,'-k');
plot(tx_frame, dhi_frame,'^b','MarkerSize',10);
plot(tx_frame, dlo_frame,'vb','MarkerSize',10);
plot( kron(tx_frame,[1 1 1]) , kron(dlo_frame,[1 0 nan]) + kron(dhi_frame,[0 1 nan]), '-b');       % link note_2
plot( tx_major , sx_major, 'pr','MarkerSize',20);
plot( tx_frame,sx_N_frame, 'dc','MarkerSize',10);
plot(kron(tx_frame,[1 1 1]), kron(sx_N_frame,[1 0 nan]) + kron(sx_frame,[0 1 nan]),'-c');
hold off