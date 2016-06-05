%Discrimination with sampled signal

%Name = 'a07m';
Name = '3987834m';

load(strcat(Name, '.mat'));
fid = fopen(strcat(Name, '.info'), 'rt');
fgetl(fid);
fgetl(fid);
fgetl(fid);
[interval] = sscanf(fgetl(fid), 'Sampling frequency: %f Hz  Sampling interval: %f sec');
interval = interval(2);

fclose(fid);

t = (1:7500) * interval;          %timeline
s = val(5,1:7500);
s  = (s  - mean(s ))/sqrt(var(s ));

% Quantisize
dt = 8e-3;                            %t_sample
quant = 1e-4;                          %vertical step
subels = 1:round(dt/interval):length(t);
t_spl = t(subels); s_spl = s(subels); %sample timeline
s_spl = quant*floor(s_spl/quant);  %sampled value

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
    
    for i = 1:floor(0.25/dt)         % search for maximum slope 0.5 s around s_max
        if d_spl(kx(k)+1 - i) >= d_spl(kx(k)+1)
            dhi(k) = d_spl(kx(k)+1 - i) ;
        end
        if d_spl(kx(k)+1 + i) <= d_spl(kx(k)+1)
            dlo(k) = d_spl(kx(k)+1 + i);
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

for k = 1:length(kx) - 1
    if (kx(k+1)-kx(k)) <= floor((1/3.5)/dt)
        if sx(k) > sx(k+1)                  % discard minor peaks with a frequency f > 3.5 Hz (BPM_max = 210)
            kx_(k+1) = 0;
        end
    end
    if (kx(k+1)-kx(k)) <= floor((1/dt))     % discard minor amplitude peaks
        if sx(k+1) < 0.8 * sx(k)
            kx_(k+1) = 0;
        end
    end
end

kx_major_index= find(kx_(1:end));
kx_major = kx_(kx_major_index);

sx_major = s_spl(kx_major + 1);
tx_major = tx(kx_major_index);

% - note 2

dhi_max = max(dhi);
dlo_min = min(dlo);

delta_note2 = dhi - dlo;
normalisation = normlist(delta_note2);      % standard score of delta_note2

kd = zeros(1,length(kx_major));

for k = 1:length(kx_major)
    if (abs(normalisation(k)) <= 1)
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
    ,tx_major,sx_major,'rp' ... % Major peaks
    ,tx,dhi,'c^' ...            % d_max_l
    ,tx,dlo,'cv' ...            % d_max_r
    ,kron(tx,[1 1 1]) , kron(sx,[0 1 nan]) , 'c-' ...                              % link note_1
    ,kron(tx,[1 1 1]) , kron(dlo,[1 0 nan]) + kron(dhi,[0 1 nan]) , 'c-' ...       % link note_2
    );

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
    ,'Location','northeastoutside');


% ,t_d,s_d,'rp' ...           % detected peaks
%  ,'detected peaks'...