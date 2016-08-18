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
dt = 1/20;                           % sampling time: dt >> dt0
t_int = dt * (1/3);                  % integration time: dt0 <= t_int < dt
quant = 0.1;                         % LSB: vertical step

[t,s] = integration(t0,s0,dt0,dt,t_int,quant,0);

%  - Peaks identification -
[kx,tx,sx, dhi,dlo, td,d, kx_n,tx_N,sx_N, note_x] = signal_peaks(t,s);

frame_init =14; frame_end = 21;

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


add = 0;
mult = 0;
div = 0;
comp = 0;
isnan_comp = 0;
abs = 0;

sort_length = 0;
sort_count = 0;

periodicity_count = 0;

%   - Clustering according to minimum variance of note_x -
for i = 1:length(kx)
    comp = comp + 1;
    
    if kx_(i) ~= 0
        comp = comp + 1;
        
        idx(1,i) = kx(i);
        per(1,i) = tx(i);
        note(1,i) = note_x(i);
        j = 2;
        
        for k = i + 1 : length(kx)
            comp = comp + 1;
            
            if  similarity(note_x(i),note_x(k),'variance') < eps;
                comp = comp +1;
                
                idx(j,i) = kx(k);
                per(j,i) = tx(k);
                note(j,i) = note_x(k);
                
                kx_(k) = 0;
                j = j+1;
                add = add +1;
            else
                comp = comp + 1;
                
                idx(j,i)=nan;
                per(j,i)=nan;
                note(j,i)=nan;
                
                j = j+1;
                add = add +1;
            end
        end
        comp = comp + 1;
    end
end
comp = comp +1;

%   - Remove columns -
zero = find(~idx(1,:));
comp = comp + length(idx);

for k = 1:length(zero)
    comp = comp + 1;
    
    idx(:,zero(k))=[];
    per(:,zero(k))=[];
    note(:,zero(k))=[];
    zero = bsxfun(@minus ,zero,ones(1,length(zero))) ;
    add = add + 1;
end
comp = comp + 1;

%   - Create cells: kx-tx-note_x -
L = size(idx);
for k = 1:L(2)
    comp = comp + 1;
    
    NAN_ = ~isnan(idx(:,k));         % extract non NAN value of idx, per, note
    NAN_idx = idx(NAN_,k);
    NAN_per = per(NAN_,k);
    NAN_note = note(NAN_,k);
    isnan_comp = isnan_comp + 1;
    
    
    idx_ = NAN_idx(find(NAN_idx));     % extract non zero value of idx, per, note
    per_ = NAN_per(find(NAN_per));
    note_ = NAN_note(find(NAN_note));
    
    clust_cell_temp{k,1} = idx_;            % kx
    clust_cell_temp{k,2} = per_;            % tx
    clust_cell_temp{k,3} = note_;           % note_x
    
    if clust_cell_temp{k,3} ~= 0            % case note_x = 0
        comp = comp + 1;
        
    NOTE_mean(k) = mean(clust_cell_temp{k,3});
    add = add + length(clust_cell_temp{k,3});
    div = div + 1;
    
    else 
        comp = comp + 1;
        
        NOTE_mean(k)=0;
    end    
    
    clear NAN_ NAN_idx NAN_per NAN_note idx_ per_ note_
end
comp = comp + 1;

% Sort clusters by NOTE
NOTE = sort(NOTE_mean,'descend');                               % average of cluster note_x sorted
sort_length = sort_length + length(NOTE_mean);
sort_count = sort_count + 1;

for k = 1:L(2)
    comp = comp + 1;
    
    clust_cell{k,1} = clust_cell_temp{NOTE_mean == NOTE(k),1};      % kx sorted
    clust_cell{k,2} = clust_cell_temp{NOTE_mean == NOTE(k),2};      % tx sorted
    clust_cell{k,3} = clust_cell_temp{NOTE_mean == NOTE(k),3};      % note_x sorted
    
end
comp = comp + 1;

%   - Cluster notation: size, tx periodicity, note_x -
for k = 1:L(2)
    comp = comp + 1;
    if length(clust_cell{k,1}) > 2
        comp = comp +1;
        
        SIZE(k) = length(clust_cell{k,1});                                 % size
        [PER_T(k),PER_eps(k), PER_R(k)] = periodicity(clust_cell{k,2});    % tx periodicity
        
        periodicity_count = periodicity_count + 1 ;
        comp = comp + 2*(length(clust_cell{k,2})+1);
        mult = mult + length(clust_cell{k,2}) + 8;
        add = add + 3*length(clust_cell{k,2}) + 4;
        div = div +5;
        
        if PER_eps(k) <= 0.1            % best periodicity note set to 0.1 otherwise increase too much the cluster note
            comp = comp + 1;
            
            PER_eps(k) = 0.1;
        end
        comp = comp + 1;
        
        clust_note(k) = (0.7 * NOTE(k) + 0.2 * SIZE(k)) / (PER_eps(k)/0.1);          % cluster note EMPIRICAL
        mult = mult + 3; 
        div = div + 3;
        add = add + 1;
        
    else
        comp = comp + 1;
        
        SIZE(k) = length(clust_cell{k,1});
        PER_T(k) = 0; PER_eps(k) = 0; PER_R(k) = 0;
        clust_note(k) = 0;
        
    end
end
comp = comp + 1;

tbl_note = table([1:L(2)]', SIZE', PER_T',PER_eps', PER_R', NOTE', clust_note','VariableNames',{'Cluster','Size','T','eps','R','Note_x','Cluster_note'});

tbl_complexity = table(comp',mult',add',div',abs',sort_count',sort_length',isnan_comp',periodicity_count','VariableNames',{'comp','mult','add','div','abs','sort_c','sort_l','isnan_comp','periodicity_count'});


