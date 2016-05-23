%function plotATM(Name)
% Name = '3801060_0007m';
% Name = '3899985_0005m';
Name = 'a07m';


% Taken from plotATM

load(strcat(Name, '.mat'));
fid = fopen(strcat(Name, '.info'), 'rt');
fgetl(fid); fgetl(fid); fgetl(fid);
[interval] = sscanf(fgetl(fid), 'Sampling frequency: %f Hz  Sampling interval: %f sec');
interval = interval(2);
% fgetl(fid);
% for i = 1:size(val, 1)
%   [row(i), signal(i), gain(i), base(i), units(i)]=strread(fgetl(fid),'%d%s%f%f%s','delimiter','\t');
% end
fclose(fid);

% val(val==-32768) = NaN;
t = find(val~=val(end));
val = val(1:t(end));

t = (1:length(val)) * interval;


% Quantisize

dt = .05; quant = 10;
subels = 1:round(dt/interval):length(val);
subt = t(subels); subval = val(subels);
subval = quant*floor(subval/quant);

plot(t', val', subt', subval', 'o--');

xlabel('Time (sec)');
% grid on


%%

DT = 1;  Ds = round(DT/interval);

kmm = 1:(floor(length(val)/Ds)-1);
vmin  = zeros(size(kmm)); vmax  = vmin;
vmint = vmin; vmaxt = vmin;
for k = kmm
    [vmin(k),vmint(k)] = min(val( (k*Ds) : ((k+1)*Ds) ));
    [vmax(k),vmaxt(k)] = max(val( (k*Ds) : ((k+1)*Ds) ));
end
vmint = (vmint + kmm*Ds) * interval;
vmaxt = (vmaxt + kmm*Ds) * interval;

plot(t', val', vmint', vmin', 'o--', vmaxt', vmax', 'o--');

%end


%%

dval  = conv(val,  .5*[1 -1], 'valid');
ddval = conv(dval, .5*[1 -1], 'valid');

val_   = (val  - mean(val))/(max(val)-min(val));
dval_  =  dval / (max( dval) - min( dval));
ddval_ = ddval / (max(ddval) - min(ddval));

plot(t', val_', (t(1:end-1)+0.5*interval)', dval_', t(2:end-1)', ddval_');

%%

plot(subt', subval' - mean(subval), conv(subt,[.5 .5], 'valid')', conv(subval,[1 -1], 'valid'));

%%

db4_lpf = [1  3  3  1] + sqrt(3)*[ 1  1 -1 -1];
db4_hpf = [1 -3  3 -1] + sqrt(3)*[-1  1  1 -1];
w_ = dwt(val, db4_lpf, db4_hpf, 5);
t_ = { (t(1:2:end-1)+t(2:2:end))/2 };
for k = 2:5
    t_{k} = (t_{k-1}(1:2:end-1)+t_{k-1}(2:2:end))/2;
end
for k = 1:5
    t_{k} = t_{k}(1:length(w_{k}));
end
t_{6} = t_{5};

plot(t, (val - mean(val))/sqrt(var(val)))
hold on
for k = 1:6
    [ size(t_{k}) ; size(w_{k}) ]
    plot(t_{k}, (w_{k} - mean(w_{k}))/sqrt(var(w_{k})))
end

legend({'PPG','DB1','DB2','DB3','DB4','DB5 bf','DB5 hf'})

%%

db4_lpf = [1  3  3  1] + sqrt(3)*[ 1  1 -1 -1];
db4_hpf = [1 -3  3 -1] + sqrt(3)*[-1  1  1 -1];
w_ = dwt(subval, db4_lpf, db4_hpf, 5);
t_ = { (subt(1:2:end-1)+subt(2:2:end))/2 };
for k = 2:5
    t_{k} = (t_{k-1}(1:2:end-1)+t_{k-1}(2:2:end))/2;
end
for k = 1:5
    t_{k} = t_{k}(1:length(w_{k}));
end
t_{6} = t_{5};

plot(subt, (subval - mean(subval))/sqrt(var(subval)))
hold on
for k = 1:6
    plot(t_{k}, (w_{k} - mean(w_{k}))/sqrt(var(w_{k})))
end

legend({'PPG','DB1','DB2','DB3','DB4','DB5 bf','DB5 hf'})

%%

db4_lpf = [1  3  3  1] + sqrt(3)*[ 1  1 -1 -1];
db4_hpf = [1 -3  3 -1] + sqrt(3)*[-1  1  1 -1];
w_ = wt(subval, db4_lpf, db4_hpf, 5);
t_ = { conv( subt ,[1 1 1 1]/4,'valid') };
for k = 2:5
    t_{k} = conv( t_{k-1} ,[1 1 1 1]/4,'valid');
end
t_{6} = t_{5};

plot(subt, (subval - mean(subval))/sqrt(var(subval)))
hold on
for k = 1:6
    [ size(t_{k}) ; size(w_{k}) ]
    plot(t_{k}, (w_{k} - mean(w_{k}))/sqrt(var(w_{k})))
end

legend({'PPG','DB1','DB2','DB3','DB4','DB5 bf','DB5 hf'})





%%

Name = 'a07m';
% Name = 'a00m';
load(strcat(Name, '.mat'));
fid = fopen(strcat(Name, '.info'), 'rt');
fgetl(fid); fgetl(fid); fgetl(fid);
[interval] = sscanf(fgetl(fid), 'Sampling frequency: %f Hz  Sampling interval: %f sec');
interval = interval(2);
fclose(fid);

val(val==-32768) = NaN;

t = (1:length(val)) * interval;


% Quantisize

dt = .05; quant = 1; del = round(dt/interval);
subels = 1:del:length(val);
subt = t(subels); subval = val(subels);
subval = quant*floor(subval/quant);

% Low-pass

subhpf = fliplr([1 -3  3 -1] + sqrt(3)*[-1  1  1 -1]);
subtf  = ones(size(subhpf))/length(subhpf);
hpf = fliplr([-1 -2 -1 0 1 2 1]);
% hpf = interp1( 0:(length(subhpf)-1), subhpf, 0:(1/del):(length(subhpf)-1), 'spline');
tf  = ones(size(hpf))/length(hpf);

hf  = conv( val ,hpf,'valid');
hft = conv( t   ,tf,  'valid');
subhf  = conv( subval ,subhpf,'valid');
subhft = conv( subt   ,subtf, 'valid');

plot(t, normlist(val), 'x-', subt, normlist(subval), 'o', hft, normlist(hf), 'x-', subhft, normlist(subhf), 'x--');


%%
plot(t, normlist(val), 'x-', subt, normlist(subval), 'o', subhft, normlist(subhf), 'x-', ...
     conv(subhft,[1 1 1]/3,'valid'), normlist(conv(subhf,[-1 2 -1],'valid')), 'x-');

 %%

subd   = conv(subval,[ 1 -1],'valid');
subdd  = conv(subd  ,[ 1 -1],'valid');
subdt  = conv(subt  ,[.5 .5],'valid');
subddt = conv(subdt ,[.5 .5],'valid');
 
plot( t, normlist(val), 'x-', subt, normlist(subval), 'o' ...
    , subdt,  normlist(subd ), 'x-' ...
    , subddt, normlist(subdd), 'x-' ...
    );

 legend({'Signal','Sampled','D','DD'})

%%


hold on
for k = 50:5:100
    plot(xcorr(subval(k-49:k),subval(k+1:k+20)))
end




%%

fl = what;
fl_ = {};
for f = fl.mat'
    if exist([f{1}(1:end-3) 'info'],'file'); fl_{length(fl_)+1} = f{1}; end
end
fl_




%%  Define input

fname = 'a07m';
dt0 = .008;
load([fname '.mat']);               % Read signal
s0 = val(19438:19813);
s0 = (s0 - mean(s0))/sqrt(var(s0)); % Normalize
t0 = (0:length(s0)-1) * dt0;      % Define timeline

dt = 1/20; dNds = .05;
t = t0(1):dt:t0(end);
s = dNds * floor( interp1(t0,s0,t) / dNds );

%%  

%%

d1 =  s(2:end) -  s(1:end-1);   t1 = (  t(2:end) +  t(1:end-1) ) / 2;
d2 = d1(2:end) - d1(1:end-1);   t2 = ( t1(2:end) + t1(1:end-1) ) / 2;
d3 = d2(2:end) - d2(1:end-1);   t3 = ( t2(2:end) + t2(1:end-1) ) / 2;

ms  = d1 > 0;
ms  = find(ms(1:end-1) & ~ms(2:end));
xd1 = t1(ms) + (t1(ms+1)-t1(ms)) .* d1(ms)./(d1(ms)-d1(ms+1));
ms  = ms+1;

md1 = d2 > 0;
md1 = find(md1(1:end-1) & ~md1(2:end));
xd2 = t2(md1) + (t2(md1+1)-t2(md1)) .* d2(md1)./(d2(md1)-d2(md1+1));

Maxd2 = d2(md1);
Mind2 = d2(md1+1);
for k = 1:length(md1)
    i = md1(k)-1;   while i > 0          && d2(i) >= Maxd2(k); Maxd2(k) = d2(i); i = i-1; end
    i = md1(k)+2;   while i < length(d2) && d2(i) <= Mind2(k); Mind2(k) = d2(i); i = i+1; end
end

md1 = md1+1;


plot( t,s,'x--', t1,d1,'x--', t2,d2,'x--' ...
    , xd1,s(ms),'o' ...
    , xd2,d1(md1),'o' , xd2,Maxd2,'^' , xd2,Mind2,'v' ...
    , [t(1) t(end)],[0 0],':k')

%%

xx = 3;
while xx>0; xx = xx-1; end
xx



%%

Name = '3900497mB';
load(strcat(Name, '.mat'));
fid = fopen(strcat(Name, '.info'), 'rt');
fgetl(fid); fgetl(fid); fgetl(fid);
[interval] = sscanf(fgetl(fid), 'Sampling frequency: %f Hz  Sampling interval: %f sec');
interval = interval(2);
fgetl(fid);
for i = 1:size(val, 1)
    [~,signal,~,~,~] = strread(fgetl(fid),'%d%s%f%f%s','delimiter','\t');
    if signal == 'II'; ecg = i; end
    
        
%   [row(i), signal(i), gain(i), base(i), units(i)]=strread(fgetl(fid),'%d%s%f%f%s','delimiter','\t')
end
fclose(fid);
