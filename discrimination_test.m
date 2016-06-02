%% Discrimination
Name = 'a07m';

load(strcat(Name, '.mat'));
fid = fopen(strcat(Name, '.info'), 'rt');
fgetl(fid);
fgetl(fid);
fgetl(fid);
[interval] = sscanf(fgetl(fid), 'Sampling frequency: %f Hz  Sampling interval: %f sec');
interval = interval(2);

fclose(fid);

t = (1:500) * interval;          %timeline
s = val(1:500);
s  = (s  - mean(s ))/sqrt(var(s ));

% Quantisize
dt = 8e-2;                            %t_sample
quant = 0.1;                          %vertical step
subels = 1:round(dt/interval):length(t);
subt = t(subels); subval = s(subels); %sample timeline
subval = quant*floor(subval/quant);  %sampled value

% Derivative, local maxima s_max, maximum slope around s_max
d = s(2:end) -  s(1:end-1);
td = (  t(2:end) +  t(1:end-1) ) / 2;

kx = d > 0;
kx = find(kx(1:end-1) & ~kx(2:end));       % k_x:index where d>0 and k_x+1<=0

sx = s(kx+1);
tx = td(kx) + (td(kx+1)-td(kx)) .* d(kx)./(d(kx)-d(kx+1));

dhi = d(kx);
dlo = d(kx+1);

for k = 1:length(kx)
    i = kx(k)-1;   while i > 0         && d(i) >= dhi(k); dhi(k) = d(i); i = i-1; end    % search for local maxima at left
    i = kx(k)+2;   while i < length(d) && d(i) <= dlo(k); dlo(k) = d(i); i = i+1; end    % search for local minima at right
end

% Filter
f = [-0.5 0.5];
ft = conv( t , ones(size(f)) , 'valid' ) / length(f) ;
fs = conv( s , fliplr(f)     , 'valid' ) ;

% Discrimination
dhi_max = max(dhi);
dlo_min = min(dlo);
kd = zeros(1,length(kx));

for j = 1:length(kx)
    if (dhi(j) >= dhi_max/2) && (dlo(j) <= dlo_min/2)
        kd(j)=kx(j);
    end
end

kd_index = find(kd(1:end));
kd = kd(kd_index);

t_d=tx(kd_index);

for l = 1:length(kd_index)
   s_d(l) = s(kd(l)); 
end


% Plot
plot(t, s,'k-'...               % siganl s
    ,subt, subval,'bo--'...     % sampled signal s_n
    ,td,d,'gx--'...             % d(s)
    ,tx,sx,'cp' ...             % s_max
    ,tx,dhi,'c^' ...            % d_max_l
    ,tx,dlo,'cv' ...            % d_max_r
    ,kron(tx,[1 1 1]) , kron(sx,[0 1 nan]) , 'c-' ...                              % link note_1
    ,kron(tx,[1 1 1]) , kron(dlo,[1 0 nan]) + kron(dhi,[0 1 nan]) , 'c-' ...       % link note_2
    ,ft,fs,'m+'...              % filter
    ,t_d,s_d,'rd' ...
    );
xlabel('Time (sec)');


