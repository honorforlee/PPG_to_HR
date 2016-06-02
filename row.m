%%
s = [2 0 1 -5 2 3 -8 0 3 2 5 1 -5 1 3];
d = s(2:end) -  s(1:end-1);

t = (1:length(s)) * 1/(125);    %timeline 
td = (  t(2:end) +  t(1:end-1) ) / 2;       

kx = d > 0;                                 
kx = find(kx(1:end-1) & ~kx(2:end));

sx = s(kx+1);                               % asign s_max to next index
tx = td(kx) + (td(kx+1)-td(kx)) .* d(kx)./(d(kx)-d(kx+1));
% test = (s(kx+1) .*td(kx+1) - s(kx) .*td(kx+1)+s(kx+1) .*td(kx)-s(kx+2) .*td(kx))./(2 .*s(kx+1)-s(kx)-s(kx+2));


dhi = d(kx);
dlo = d(kx+1);
for k = 1:length(kx)
    i = kx(k)-1;   while i > 0         && d(i) >= dhi(k); dhi(k) = d(i); i = i-1; end    % search for local maxima at left
    i = kx(k)+2;   while i < length(d) && d(i) <= dlo(k); dlo(k) = d(i); i = i+1; end    % search for local minima at right
end