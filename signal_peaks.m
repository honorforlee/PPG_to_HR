% Ivan Ny Hanitra - Master thesis
%       -- Compute signal peaks--

function [kx,tx,sx, dhi,dlo, td,d, kx_n,tx_N,sx_N, note_x] = signal_peaks(t,s,signal,varargin)
%   - Derivative -
d = s(2:end) -  s(1:end-1);
%td = (  t(2:end) +  t(1:end-1) ) / 2;      % timeline of derivative shifted by t_sample/2
td = t(2:end);

kx = d > 0;
kx = find(kx(1:end-1) & ~kx(2:end));       % k_{x}:index where d > 0; d( k_{x} + 1 ) <= 0

%   - Local maxima sx, maximum slope around sx -
[tx,sx, dhi,dlo, kx_n,tx_N,sx_N, note_1,note_2, note_3] = peaks_processing(t,s,td,d,kx);

if nargin == 2
    %   - Events rating -
    note_x = 0.1*note_1 + 0.1*note_2 + 0.8*note_3;
    
else
    if signal == 'ecg'
        note_x = 0.5*note_2+0.5*note_3;
    end
end