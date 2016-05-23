function y = normlist(x)
    y = (x-mean(x))/sqrt(var(x));
end
