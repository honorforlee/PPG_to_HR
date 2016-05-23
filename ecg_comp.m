%%

xx = {};
xx = [xx 'foo'];
xx = [xx 'bar'];
xx{1}




%%
% fname = 'a00m';
fname = '3900497mB';
% fname = '3900679m';
% Read sampling frequency
load([fname '.mat']);
fid = fopen([fname '.info'], 'rt');
for i = 1:6; s = fgetl(fid); end
while s
    s
    [~,signal,~,~,~] = strread(s,'%d%s%f%f%s','delimiter','\t');
    if strcmp(signal,'II'); fprintf('II'); break; end
    s = fgetl(fid);
end
fclose(fid);

%%

s = {};
s{length(s)+1} = 'foo';
s{length(s)+1} = 'bar';
for x = s
    x{1}
end



%% IN FILE

fname = '3900497m';
fname = '3900679m';
% Read sampling frequency
load([fname '.mat']);
fid = fopen([fname '.info'], 'rt');
fgetl(fid); fgetl(fid); fgetl(fid);
dt0 = sscanf(fgetl(fid), 'Sampling frequency: %f Hz  Sampling interval: %f sec');
dt0 = dt0(2);
fgetl(fid);
for i = 1:size(val, 1)
    [~,signal,~,~,~] = strread(fgetl(fid),'%d%s%f%f%s','delimiter','\t');
    if strcmp(signal,'II');    ecg = val(i,:); end
    if strcmp(signal,'PLETH'); s0  = val(i,:); end
end
fclose(fid);
% Remove tail
k = find(s0 ~= s0(end),1,'last'); s0 = s0(1:k); ecg = ecg(1:k);
% Normalize
s0  = (s0  - mean(s0 ))/sqrt(var(s0 ));
ecg = (ecg - mean(ecg))/sqrt(var(ecg));
% Define timeline
t0 = (0:length(s0)-1) * dt0;

val = [];

%% SAMPLE

dt = 1/20; dNds = .05;
t = t0(1):dt:t0(end);
s = dNds * floor( interp1(t0,s0,t) / dNds );

%%

tecg = (ecg(2:end) - ecg(1:end-1)) > 0;
tecg = find(tecg(1:end-1) & ~tecg(2:end) & ecg(2:end-1)>1) + 1;
tecg = t0(tecg);
tecg_x = kron(tecg,[1 1]) ;
tecg_y = 10*(2-abs(3 - 2*mod(0:length(tecg_x)-1,4))) ;


%%



%%

plot(t0,s0,t,s,t0,ecg,.5*(t0(1:end-1)+t0(2:end)),ecg(2:end)-ecg(1:end-1),tecg_x,tecg_y,':k')
