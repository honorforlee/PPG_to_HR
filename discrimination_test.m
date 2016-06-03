%% Discrimination with siganl

%Name = 'a07m';
Name = '3916979m';

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
dt = 8e-2;                            %t_sample
quant = 0.1;                          %vertical step
subels = 1:round(dt/interval):length(t);
t_spl = t(subels); s_spl = s(subels); %sample timeline
s_spl = quant*floor(s_spl/quant);  %sampled value

% Derivative, local maxima s_max, maximum slope around s_max
d = s(2:end) -  s(1:end-1);
d_spl = s_spl(2:end) -  s_spl(1:end-1);
td = (  t(2:end) +  t(1:end-1) ) / 2;

kx = d > 0;
kx = find(kx(1:end-1) & ~kx(2:end));       % k_x: index where d>0 and k_{x} +1<=0
kx_d = zeros(1,length(kx));                % kx_d: index of major peak

sx = s(kx+1);
tx = td(kx) + (td(kx+1)-td(kx)) .* d(kx)./(d(kx)-d(kx+1));
sx_d = sx;                                 % major peak

dhi = d(kx);
dlo = d(kx+1);

% for k = 1:length(kx)
%     i = kx(k)-1;   while i > 0         && d(i) >= dhi(k); dhi(k) = d(i); i = i-1; end    % maximum slope before peak
%     i = kx(k)+2;   while i < length(d) && d(i) <= dlo(k); dlo(k) = d(i); i = i+1; end    % maximum slope after peak
% end

for k = 1:length(kx)
 
    for i = 1:floor(0.25/dt)         % search for maximum slope 0.5 s around s_max
           if d(kx(k)-i) >= d(kx(k))
               dhi(k) = d(kx(k)- i) ;
           end
           if d(kx(k) +i) <= d(kx(k))
               dlo(k) = d(kx(k)+i);
           end
         i = i+1;
    end
       
       for i_ = 1:floor(0.5/dt)       % search for maxima 1 s around s_max
           if s(kx(k)-i_) >= s(kx(k))
               sx_d(k) = s(kx(k)- i_) ;
           else   kx_d(k) = 1;
                           
           end
           if s(kx(k) +i_) >= s(kx(k))
               sx_d(k) = s(kx(k)+ i_);
           else kx_d(k) = 1;
           end
         i_ = i_+1;
       end
end

% Filter
f = [-0.5 0.5];
ft = conv( t , ones(size(f)) , 'valid' ) / length(f) ;
fs = conv( s , fliplr(f)     , 'valid' ) ;

% Discrimination of peaks
dhi_max = max(dhi);
dlo_min = min(dlo);
delta_note2 = mean(dhi) - mean(dlo);

kd = zeros(1,length(kx));

% for j = 1:length(kx)
%     if (dhi(j) >= dhi_max/4) && (dlo(j) <= dlo_min/4)
%         kd(j)=kx(j);
%     end
% end

for j = 1:length(kx)
    if (dhi(j)-dlo(j) >= delta_note2)
        kd(j)=kx(j);
    end
end

kd_index = find(kd(1:end));
kd = kd(kd_index);        % indices of discriminated peaks 

t_d=tx(kd_index);                      

for l = 1:length(kd_index)
   s_d(l) = s(kd(l));     % value of discriminated peaks
end

% measure heart rate


% Plot
plot(t, s,'k-'...               % siganl s
    ,t_spl, s_spl,'bo--'...     % sampled signal s_n
    ,td,d,'gx--'...             % d(s)
    ,tx,sx_d,'cd' ...             % s_max
    ,tx,dhi,'c^' ...            % d_max_l
    ,tx,dlo,'cv' ...            % d_max_r
    ,kron(tx,[1 1 1]) , kron(sx,[0 1 nan]) , 'c-' ...                              % link note_1
    ,kron(tx,[1 1 1]) , kron(dlo,[1 0 nan]) + kron(dhi,[0 1 nan]) , 'c-' ...       % link note_2
    ,ft,fs,'m+'...              % filter
    ,t_d,s_d,'rp' ...
    );
xlabel('Time (sec)');


%% Discrimination with sampled signal

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

% Derivative, local maxima s_max, maximum slope around s_max
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
           if d_spl(kx(k)-i) >= d_spl(kx(k))
               dhi(k) = d_spl(kx(k)- i) ;
           end
           if d_spl(kx(k) +i) <= d_spl(kx(k))
               dlo(k) = d_spl(kx(k)+i);
           end
         i = i+1;
    end
       
       for i_ = 1:floor(0.5/dt)       % search for maxima 1 s around s_max
           if s_spl(kx(k)-i_) >= s_spl(kx(k))
               sx_d(k) = s_spl(kx(k)- i_) ;
           end
           if s(kx(k) +i_) <= d(kx(k))
               sx_d(k) = s_spl(kx(k)+ i_);
           end
         i_ = i_+1;
       end
end

% Filter
f = [-0.5 0.5];
ft = conv( t_spl , ones(size(f)) , 'valid' ) / length(f) ;
fs = conv( s_spl , fliplr(f)     , 'valid' ) ;

% Discrimination of peaks

dhi_max = max(dhi);
dlo_min = min(dlo);

delta_note2 = dhi - dlo;                        
normalisation = normlist(delta_note2);      % standard score of delta_note2 

kd = zeros(1,length(kx));

for j = 1:length(kx)
    if (abs(normalisation(j)) <= 1)
    kd(j)=kx(j);
    end
end

kd_index = find(kd(1:end));
kd = kd(kd_index);        % indices of discriminated peaks 

t_d=tx(kd_index);                      

for l = 1:length(kd_index)
   s_d(l) = s_spl(kd(l));     % value of discriminated peaks
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
    ,tx,sx,'cd' ...             % s_max
    ,tx,dhi,'c^' ...            % d_max_l
    ,tx,dlo,'cv' ...            % d_max_r
    ,t_d,s_d,'rp' ...           % detected peaks
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
      ,'s_{max}'...
      ,'maximum positive slope around s_{max}'...
      ,'maximum negative slope around s_{max}'...
      ,'detected peaks'...
      ,'Location','northeastoutside');


