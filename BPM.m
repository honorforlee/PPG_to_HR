% Ivan NY HANITRA - Master thesis
%       -- Compute signal BPM with dt=interval  --


Name = '3987834m';           % BPM = 78
%Name = '3801060_0007m';     % BPM = 95

load(strcat(Name, '.mat'));
fid = fopen(strcat(Name, '.info'), 'rt');
fgetl(fid);
fgetl(fid);
fgetl(fid);
[interval] = sscanf(fgetl(fid), 'Sampling frequency: %f Hz  Sampling interval: %f sec');
interval = interval(2);              % data acquisition rate (interval = 1/f_spl_u = 0.5903 ms in practice)

fclose(fid);

t = (1:length(val)) * interval;              % timeline
s = val(5,1:length(val));
s  = (s  - mean(s ))/sqrt(var(s ));          % rescale s on 0 (standard score of signal)

% Quantisize
dt = interval;                       % sampling time: dt >> interval
%t_int = dt * (1/3);                 % integration time: interval <= t_int < dt
quant = 1e-4;                        % LSB: vertical step

subels = (1:round(dt/interval):length(t));
t_spl = t(subels);                    % sample timeline
s_spl = s(subels);
s_spl = quant*floor(s_spl/quant);       % quantisation

% Derivative, local maxima sx, maximum slope around sx
d = s(2:end) -  s(1:end-1);
td = (  t(2:end) +  t(1:end-1) ) / 2;

d_spl = s_spl(2:end) -  s_spl(1:end-1);
td_spl = (  t_spl(2:end) +  t_spl(1:end-1) ) / 2;

kx = d_spl > 0;
kx = find(kx(1:end-1) & ~kx(2:end));       % k_{x}:index where d>0 and k_{x} +1<=0

sx = s_spl(kx+1);                          % local maxima
tx = td_spl(kx) + (td_spl(kx+1)-td_spl(kx)) .* d_spl(kx)./(d_spl(kx)-d_spl(kx+1));      % linear interpolation of dhi and dho to get tx (@zero crossing)

% dhi = d_spl(kx);
% dlo = d_spl(kx+1);
% 
% for k = 1:length(kx)                % search for maximum slope 0.5 s around s_max
%     
%     for i = 1:floor(0.25/dt)         
%         if (kx(k)+1 - i) > 0
%             if d_spl(kx(k)+1 - i) >= dhi(k)
%                 dhi(k) = d_spl(kx(k)+1 - i) ;
%             end
%         end
%         if (kx(k)+1 + i) < length(d_spl)
%             if d_spl(kx(k)+1 + i) <= dlo(k)
%                 dlo(k) = d_spl(kx(k)+1 + i);
%             end
%         end
%         i = i+1;
%     end
% end

dhi_ = d_spl(kx);
dlo_ = d_spl(kx+1);

for k = 1:length(kx)
    i = kx(k)-1;   while i > 0         && d(i) >= dhi_(k); dhi_(k) = d(i); i = i-1; end    % search for local maxima at left
    i = kx(k)+2;   while i < length(d) && d(i) <= dlo_(k); dlo_(k) = d(i); i = i+1; end    % search for local minima at right
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
    if (kx(k+1)-kx(k)) <= floor((1/3.5)/dt)
        if sx(k+1) < sx(k)                  % discard minor peaks with a frequency f > 3.5 Hz (BPM_max = 210)
            kx_(k+1) = 0;
        end
    end
    
    if (kx(k+1)-kx(k)) <= floor((1/dt))     % discard minor amplitude peaks
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

delta_note2 = dhi_ - dlo_;
normalisation = normlist(delta_note2);      % standard score of delta_note2

kd = zeros(1,length(kx_major));

for k = 1:length(kx_major)
    if (normalisation(k) >= -2)   % discard peaks wich delta_note2 is 1 standard deviation below the mean delta_note2
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

