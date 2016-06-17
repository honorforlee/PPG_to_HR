%   - Load file and data -
%Name = '3987834m';     % BPM = 78
Name = '3900497m'; % BPM = 95

load(strcat(Name, '.mat'));
fid = fopen(strcat(Name, '.info'), 'rt');
fgetl(fid);
fgetl(fid);
fgetl(fid);
[interval] = sscanf(fgetl(fid), 'Sampling frequency: %f Hz  Sampling interval: %f sec');
interval = interval(2);              % data acquisition rate (interval = 1/f_spl_u = 0.5903 ms in practice)

fclose(fid);

t = (1:length(val)) * interval;            % timeline
s = val(6,1:length(val));
s  = (s  - mean(s ))/sqrt(var(s ));          % rescale s on 0 (standard score of signal)

%   - Timeline, noise, integration, quantization -
dt = 1/10;                           % sampling time: dt >> interval
t_int = dt * (1/3);                  % integration time: interval <= t_int < dt
quant = 1e-4;                        % LSB: vertical step

t_spl = t(1):dt:t(end);    

% % Noise
% frameNoise = (0:round(dt/interval))';
% frameNoise = bsxfun(@minus, subels, frameNoise);
% frameNoise_zero = find (frameNoise <= 0);
% frameNoise(frameNoise_zero) = 1;
% 
% noise= random('Normal',mean(s(frameNoise)),std(s(frameNoise)),1,length(subels));                     % Gaussian distribution (model thermal noise of finite BW)

% % Integration
% frameInteg = (0:t_int)';
% frameInteg = bsxfun(@minus, t_spl, frameInteg);
% frameInteg_zero = find (frameInteg <= 0);
% frameInteg(frameInteg_zero) = 1;                       % t_int < dt
% 
% %s_spl = mean( vertcat(s(frameInteg),noise) );          % sampled signal = average of Nint last values + noise during dt
% s_spl = mean( interp1(t,s,frameInteg) );                          % signal with no noise

% s_spl(1) = s(1);
% for k = 2:length(t_spl) 
%    
%     s_spl(k) = mean ( [s(floor( (t_spl(k-1)+ dt-t_int)/interval ): s(floor( (t_spl(k-1)+ dt-t_int)/interval ) :s(floor(dt/interval))) ] );
%             
% end 

%s_spl = interp1(t,s,t_spl); 
%s_spl = quant*floor(s_spl/quant);                      % quantization

% plot(t, s,'k-','MarkerSize',12,'LineWidth',1);               % siganl s
% hold on
% plot(t_spl, s_spl,'ko--','MarkerSize',10,'LineWidth',1);     % sampled signal s_n
% hold off
%%

for k = 1:length(t_spl)-1
frameNoise (:,k) = [ floor(t_spl(k)/interval) :  floor(t_spl(k)/interval) + floor(dt/interval) ]; 
end
noise= random('Normal',mean(s(frameNoise)),std(s(frameNoise)),1,length(frameNoise));                     % Gaussian distribution (model thermal noise of finite BW)

for k = 2:length(t_spl)
index = min (length(    [floor( (t_spl(k-1)+ dt-t_int)/interval ): floor( t_spl(k)/ interval ) ]    ));
end
index = index-1;
for k = 2:length(t_spl)
frameInteg(:,k-1) = [ floor( t_spl(k)/ interval ) - index : floor( t_spl(k)/ interval ) ];    
frameInteg_(:,k-1)= s(frameInteg(:,k-1)) ;
end
frameInteg_ = vertcat(frameInteg_,noise);

s_spl(1) = s(1);
for k = 2:length(t_spl)
    s_spl(k) = mean(frameInteg_(:,k-1));
end    

