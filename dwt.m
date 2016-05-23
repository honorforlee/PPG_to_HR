function w = dwt(s, lpf, hpf, max_level)
    % s   = input signal
    % lpf = low-pass filter,  or "sacling"   (for db4: lpf = [1  3  3  1] + sqrt(3)*[ 1  1 -1 -1])
    % hpf = high-pass filter, or "wavelet"   (for db4: hpf = [1 -3  3 -1] + sqrt(3)*[-1  1  1 -1])
    % max_level = max level of the discrete wavelet transform
    % returns  { hpf level 1, hpf level 2, ..., hpf level max_level, lpf level max_level }

    lpf1 = fliplr(lpf); lpf2 = lpf1(1:2:end); lpf1 = lpf1(2:2:end);
    hpf1 = fliplr(hpf); hpf2 = hpf1(1:2:end); hpf1 = hpf1(2:2:end);
    
    w = {s};
    for l = 1:max_level
        w{l+1} = conv( w{l}(1:2:end-1) ,lpf1,'valid') + conv( w{l}(2:2:end) ,lpf2,'valid');
        w{l}   = conv( w{l}(1:2:end-1) ,hpf1,'valid') + conv( w{l}(2:2:end) ,hpf2,'valid');
    end
    
end
