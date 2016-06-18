% Ivan Ny Hanitra - Master thesis
%       -- Simulate integrating ADC: add thermal noise, integrate during t_int (share of t_sample), quantize the output --

function [t,s] = integration(t0,s0,dt0,dt,t_int,quant)
t = t0(1):dt:t0(end);                                   % timeline with new sampling frequency

% Noise
for k = 1:length(t)-1
    frameNoise (:,k) = [ floor(t(k)/dt0) :  floor(t(k)/dt0) + floor(dt/dt0) ];
end
noise= random('Normal',mean(s0(frameNoise)),std(s0(frameNoise)),1,length(frameNoise));                     % Gaussian distribution (model thermal noise of finite BW)

% Integration
if t_int ~=0
    
    for k = 2:length(t)
        index = min (length(    [floor( (t(k-1)+ dt - t_int)/dt0 ): floor( t(k)/ dt0 ) ]    ));
    end
    index = index-1;
    for k = 2:length(t)
        frameInteg(:,k-1) = [ floor( t(k)/ dt0 ) - index : floor( t(k)/ dt0 ) ];
        frameInteg_(:,k-1)= s0(frameInteg(:,k-1)) ;
    end
    
    frameInteg_ = vertcat(frameInteg_,noise);      % add Gaussian noise before integration
    
    s(1) = s0(1);
    for k = 2:length(t)
        s(k) = mean(frameInteg_(:,k-1));           % sampled signal = average of Nint last values + noise during dt
    end
    
else
    s = zeros(1,length(t));
end

s = quant * floor( s / quant );                % quantization: quant = LSB

%   - with time subdivided in indexes -
% subels = (1:round(dt/interval):length(t));
% t_spl = t(subels);                           % sample timeline
%
% % Noise
% frameNoise = (0:round(dt/interval))';
% frameNoise = bsxfun(@minus, subels, frameNoise);
% frameNoise_zero = find (frameNoise <= 0);
% frameNoise(frameNoise_zero) = 1;
%
% noise = random('Normal',mean(s(frameNoise)),std(s(frameNoise)),1,length(subels));      % Gaussian distribution (model thermal noise of finite BW)
%
% % Integration
% frameInteg = (0:round(t_int/interval))';
% frameInteg = bsxfun(@minus, subels, frameInteg);
% frameInteg_zero = find (frameInteg <= 0);
% frameInteg(frameInteg_zero) = 1;                       % t_int < dt
%
% s_spl = mean( vertcat(s(frameInteg), noise) );         % sampled signal = average of Nint last values + noise during dt
%
% s_spl = quant*floor(s_spl/quant);                      % quantization