function w = wt(s, lpf, hpf, max_level)
    % s   = input signal
    % lpf = low-pass filter,  or "sacling"   (for db4: lpf = [1  3  3  1] + sqrt(3)*[ 1  1 -1 -1])
    % hpf = high-pass filter, or "wavelet"   (for db4: hpf = [1 -3  3 -1] + sqrt(3)*[-1  1  1 -1])
    % max_level = max level of the discrete wavelet transform
    % returns  { hpf level 1, hpf level 2, ..., hpf level max_level, lpf level max_level }

    lpf = fliplr(lpf);
    hpf = fliplr(hpf);
    w = {s};
    for l = 1:max_level
        w{l+1} = conv( w{l} ,lpf,'valid');
        w{l}   = conv( w{l} ,hpf,'valid');
    end
    
end
