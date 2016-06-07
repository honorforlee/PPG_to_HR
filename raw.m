%Discrimination with sampled signal

%Name = '3987834m';
Name = '3801060_0007m';

load(strcat(Name, '.mat'));
fid = fopen(strcat(Name, '.info'), 'rt');
fgetl(fid);
fgetl(fid);
fgetl(fid);
[interval] = sscanf(fgetl(fid), 'Sampling frequency: %f Hz  Sampling interval: %f sec');
interval = interval(2);

fclose(fid);

t = (1:7500) * interval;             % timeline
s = val(1,1:7500);
s  = (s  - mean(s ))/sqrt(var(s ));  % rescale s on 0 (standard score of signal)

% Quantisize and integration
dt = 8e-3;                            % t_sample
quant = 1e-4;                         % vertical step

subels = 1:round(dt/interval):length(t);

t_spl = t(subels);      % timeline
s_spl = s(subels); 

t_int = 1e-3;

for k = 1:length(subels)-1
    
    subels_int(k) = subels(k) + (subels(k+1)-subels(k))/(dt/t_int);
    
end

subels_tot = horzcat(subels,subels_int);
subels_tot = sort(subels_tot);



for k = 1:length(subels)-1
    
    s_int;
    
end


noise = poisspdf(subels,5);
s_spl = quant*floor((s_spl + noise)/quant);  % sampled value

%%
% Derivative, local maxima s_max, maximum slope around s_max, major maxima
d = s(2:end) -  s(1:end-1);
td = (  t(2:end) +  t(1:end-1) ) / 2;

d_spl = s_spl(2:end) -  s_spl(1:end-1);
td_spl = (  t_spl(2:end) +  t_spl(1:end-1) ) / 2;

kx = d_spl > 0;
kx = find(kx(1:end-1) & ~kx(2:end));       % k_{x}:index where d>0 and k_{x} +1<=0

sx = s_spl(kx+1);           % local maxima
sx_d = sx;                  % major local maxima

tx = td_spl(kx) + (td_spl(kx+1)-td_spl(kx)) .* d_spl(kx)./(d_spl(kx)-d_spl(kx+1));

dhi = d_spl(kx);
dlo = d_spl(kx+1);

for k = 1:length(kx)
    
    for i = 1:round(0.25/dt)         % search for maximum slope 0.5 s around s_max
        if (kx(k)+1 - i) > 0
            if d_spl(kx(k)+1 - i) >= d_spl(kx(k)+1)
                dhi(k) = d_spl(kx(k)+1 - i) ;
            end
        end
        if (kx(k)+1 + i) < length(d_spl)
            if d_spl(kx(k)+1 + i) <= d_spl(kx(k)+1)
                dlo(k) = d_spl(kx(k)+1 + i);
            end
        end
        i = i+1;
    end
end


% Filter
f = [-0.5 0.5];
ft = conv( t_spl , ones(size(f)) , 'valid' ) / length(f) ;
fs = conv( s_spl , fliplr(f)     , 'valid' ) ;

% Discrimination of peaks
% - note 1

kx_ = kx;       % auxiliary array

if sx(1) < 0.5 * sx(2)
    kx_(1) = 0;
end

for k = 1:length(kx) - 1
    if (kx(k+1)-kx(k)) <= round((1/3.5)/dt)
        if sx(k+1) < sx(k)                  % discard minor peaks with a frequency f > 3.5 Hz (BPM_max = 210)
            kx_(k+1) = 0;
        end
    end
    
    if (kx(k+1)-kx(k)) <= round((1/dt))     % discard minor amplitude peaks
        if sx(k+1) < 0.5 * sx(k)
            kx_(k+1) = 0;
        end
    end
end

kx_major_index= find(kx_(1:end));
kx_major = kx_(kx_major_index);

sx_major = s_spl(kx_major + 1);
tx_major = tx(kx_major_index);

% - note 2

% dhi_= d_spl(kx_major);
% dlo_ =d_spl(kx_major + 1);
%
% for k = 1:length(kx_major)
%     i = kx_major(k)-1;   while i > 0         && d_spl(i) >= dhi_(k); dhi_(k) = d_spl(i); i = i-1; end    % search for local maxima at left
%     i = kx_major(k)+2;   while i < length(d) && d_spl(i) <= dlo_(k); dlo_(k) = d_spl(i); i = i+1; end    % search for local minima at right
% end

delta_note2 = dhi - dlo;
normalisation = normlist(delta_note2);      % standard score of delta_note2

kd = zeros(1,length(kx_major));

for k = 1:length(kx_major)           
    if (normalisation(k) >= -1)   % discard peaks wich delta_note2 is 1 standard deviation below the mean delta_note2      
        kd(k)=kx_major(k);
    end
end

kd_index = find(kd(1:end));
kd = kd(kd_index);               % indices of discriminated peaks

t_d=tx_major(kd_index);

for k = 1:length(kd_index)
    s_d(k) = s_spl(kd(k)+1);     % value of discriminated peaks
end

% Measure HR
HR = (60*length(kd))/t(end);
result = sprintf('Average heart Rate: %d bpm', HR);
disp(result);

% Plot
plot(t, s,'k-'...               % siganl s
    ,t_spl, s_spl,'bo--'...     % sampled signal s_n
    ,td_spl,d_spl,'gx--'...     % derivative of s_n
    ,ft,fs,'m+'...              % filter s
    ,tx_major,sx_major,'bp' ... % Major peaks
    ,tx,dhi,'c^' ...            % d_max_l
    ,tx,dlo,'cv' ...            % d_max_r
    ,t_d,s_d,'rp' ...           % detected peaks
    ,kron(tx,[1 1 1]) , kron(dlo,[1 0 nan]) + kron(dhi,[0 1 nan]) , 'c-' ...       % link note_2
    );
%,kron(tx,[1 1 1]) , kron(sx,[0 1 nan]) , 'c-' ...                                 % link note_1

title('Peaks discrimination for heart rate monitoring');
xlabel('Time, s');
ylabel('Arbitrary units');
legend('s: original signal'...
    ,'s_n: sampled signal'...
    ,'derivative of s_n'...
    ,'filtering of s'...
    ,'Major peaks'...
    ,'maximum positive slope around s_{max}'...
    ,'maximum negative slope around s_{max}'...
    ,'detected peaks'...
    ,'Location','northeastoutside');
