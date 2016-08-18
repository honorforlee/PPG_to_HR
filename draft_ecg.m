Name = '3900497m';
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

val(isnan(val)) = [];
t0 = (1:length(val)) * dt0;            % timeline
s0 = val(ecg,1:length(val));

s0  = (s0  - mean(s0 ))/sqrt(var(s0));        % rescale s on 0 (standard score of signal)

%  - Peaks identification -
[kx,tx,sx, dhi,dlo, td,d, kx_n,tx_N,sx_N, note_x] = signal_peaks(t0,s0,'ecg');

frame_init =14; frame_end = 21;

index_x = find(tx >= frame_init & tx <= frame_end);
sx_N_frame = sx_N(index_x);
dhi_frame = dhi(index_x);
dlo_frame = dlo(index_x);

index0 = find(t0 >= frame_init & t0 <= frame_end);
t0_frame = t0(index0);
s0_frame = s0(index0);

[kx_frame,tx_frame,sx_frame,note_x_frame] = frame_select(kx,tx,sx,note_x, frame_init,frame_end);

kx = kx_frame; tx = tx_frame; sx = sx_frame; note_x = note_x_frame;

%[kx_major,tx_major,sx_major, T,warning] = min_variance(kx_frame,tx_frame,sx_frame, note_x_frame, 0.1);
kx_ = kx;
eps = 1;

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
    
    if clust_cell_temp{k,3} ~= 0            % case note_x = 0
    NOTE_mean(k) = mean(clust_cell_temp{k,3});
    else 
        NOTE_mean(k)=0;
    end    
    
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
    
        SIZE(k) = length(clust_cell{k,1});                                 % size
        [PER_T(k),PER_eps(k), PER_R(k)] = periodicity(clust_cell{k,2});    % tx periodicity
        
        if PER_eps(k) <= 0.1            % best periodicity note set to 0.1 otherwise increase too much the cluster note
            PER_eps(k) = 0.1;
        end
        
        clust_note(k) = (NOTE(k));          % cluster note EMPIRICAL
                
end

tbl_note = table([1:L(2)]', SIZE', PER_T',PER_eps', PER_R', NOTE', clust_note','VariableNames',{'Cluster','Size','T','eps','R','Note_x','Cluster_note'});

[major_clust major_idx] = max(clust_note);
kx_major = clust_cell{major_idx,1}';
idxs = arrayfun(@(x)find(kx==x,1),kx_major);
tx_major = tx(idxs);
sx_major = sx(idxs);


%%
%   - Plots -
figure(2);
plot( tx_frame , sx_frame   , 'dc','MarkerSize',12);
hold on
plot(t0_frame,s0_frame,'-k');
plot(tx_frame, dhi_frame,'^b','MarkerSize',10);
plot(tx_frame, dlo_frame,'vb','MarkerSize',10);
plot( kron(tx_frame,[1 1 1]) , kron(dlo_frame,[1 0 nan]) + kron(dhi_frame,[0 1 nan]), '-b');       % link note_2
plot( tx_major , sx_major, 'pr','MarkerSize',20);
plot( tx_frame,sx_N_frame, 'dc','MarkerSize',10);
plot(kron(tx_frame,[1 1 1]), kron(sx_N_frame,[1 0 nan]) + kron(sx_frame,[0 1 nan]),'-c');
hold off


%%
% lpf = low-pass filter,  or "sacling"   (for db4: lpf = [1  3  3  1] + sqrt(3)*[ 1  1 -1 -1])
% hpf = high-pass filter, or "wavelet"   (for db4: hpf = [1 -3  3 -1] + sqrt(3)*[-1  1  1 -1])
% max_level = max level of the discrete wavelet transform
% returns  { hpf level 1, hpf level 2, ..., hpf level max_level, lpf level max_level }


% lpf = [1  3  3  1] + sqrt(3)*[ 1  1 -1 -1];
% hpf = [1 -3  3 -1] + sqrt(3)*[-1  1  1 -1];
%     
% lpf = fliplr(lpf);
% hpf = fliplr(hpf);

% w = {s0_};
% for l = 1:5
%     w{l+1} = conv( w{l} ,lpf,'valid');
%     w{l}   = conv( w{l} ,hpf,'valid');
% end