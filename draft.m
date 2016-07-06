tx_major = [0:1:20];
sx_major = (5-0).*rand(1,21) + 0;

zeros = [4 7 12 18];
T_temp = .5;

insert = @(a, x, n)cat(2,  x(1:n), a, x(n+1:end));

for k = 1:length(zeros)
        tx_major = insert(tx_major(zeros(k)) + T_temp,tx_major, zeros(k));
        sx_major = insert(sx_major(zeros(k)), sx_major, zeros(k) );        % inset one nan value at added time sampled
        zeros = bsxfun(@plus ,zeros,ones(1,length(zeros)));
end
