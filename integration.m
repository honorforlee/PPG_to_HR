function [t,s] = integration(t0,s0,dt0,dt,t_int,quant)
t = t0(1):dt:t0(end);                                   % timeline with new sampling frequency

% Noise
for k = 1:length(t)-1
frameNoise (:,k) = [ floor(t(k)/dt0) :  floor(t(k)/dt0) + floor(dt/dt0) ]; 
end
noise= random('Normal',mean(s0(frameNoise)),std(s0(frameNoise)),1,length(frameNoise));                     % Gaussian distribution (model thermal noise of finite BW)

% Integration
for k = 2:length(t)
    index = min (length(    [floor( (t(k-1)+ dt - t_int)/dt0 ): floor( t(k)/ dt0 ) ]    ));
end
index = index-1;
for k = 2:length(t)
    frameInteg(:,k-1) = [ floor( t(k)/ dt0 ) - index : floor( t(k)/ dt0 ) ];
    frameInteg_(:,k-1)= s0(frameInteg(:,k-1)) ;
end
frameInteg_ = vertcat(frameInteg_,noise);

s(1) = s0(1);
for k = 2:length(t)
    s(k) = mean(frameInteg_(:,k-1));               % sampled signal = average of Nint last values + noise during dt
end
 s = quant * floor( s / quant );                % quantization