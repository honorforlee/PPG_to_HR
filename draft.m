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


%%
%%
% %   - Cluster notation: size, mean note, periodicity -
% L = size(idx);
% for k = 1:L(2)
%     clust_size(k) = nnz(idx(:,k)) - sum(isnan(idx(:,k)));
%     if clust_size(k) > 2
%         NAN_idx = ~isnan(per(:,k));         % extract non NAN value of per, note
%         NAN_per = per(NAN_idx,k);
%         NAN_note = note(NAN_idx,k);
%
%         per_ = NAN_per(find(NAN_per));      % extract non zero value of per, note
%         note_ = NAN_note(find(NAN_note));
%
%         [PER_T(k),PER_eps(k), PER_R(k)] = periodicity(per_);    % note
%
%         if PER_eps(k) <= 0.01
%             PER_eps(k) = 0.01;
%         end
%         NOTE(k) = mean(note_);
%         SIZE(k) = clust_size(k);
%
%         clust_note(k) = (0.2 * NOTE(k) + 0.4 * SIZE(k)) / (PER_eps(k)/0.4);
%
%         clear NAN_idx NAN_per NAN_note per_ note_
%     else
%         NAN_idx = ~isnan(per(:,k));
%         NAN_note = note(NAN_idx,k);
%         note_ = NAN_note(find(NAN_note));
%
%         PER_T(k) = 0; PER_eps(k) = 0; PER_R(k) = 0;
%         NOTE(k) = mean(note_);
%         SIZE(k) = clust_size(k);
%         clust_note(k) = 0;
%
%         clear NAN_idx NAN_note note_
%     end
% end
%
% tbl_note = table([1:L(2)]', PER_T',PER_eps', PER_R', NOTE', SIZE', clust_note','VariableNames',{'Cluster','T','eps','R','Note_x','Size','Cluster_note'})

note_major(1) = max(clust_note);

j=1;
for k = 1:L(2)
    
    if var([note_major(1) clust_note(k)],1) < 5*eps
        clust_major(j) = k;
        j = j+1;
    end
    
end

kx_major = kx(clust_major);


%%
% plot note_note_x
for i = 1 : kmax
    figure(2);
    subplot(2,1,1);
    plot(note_note_x{i,kmax} , '.');
    hold on
end
Legend=cell(kmax,1);
for iter=1:kmax
    Legend{iter}=strcat('noteer ', num2str(iter));
end
legend(Legend);

hold off
subplot(2,1,2);
plot( note_x, '.');

base_array = cellfun(@length,note_tx);
base_max = max (base_array(:,kmax));
base = [1:base_max];

% plot note_periodicity
data = nan(base_max,kmax);

for i = 1 : kmax
    
    data(1:base_array(i,kmax),i) = note_tx{i,kmax};
    figure(3);
    tx_disp(i) = plot(base,data(:,i),'.');
    hold on
end

Legend=cell(kmax,1);
for iter=1:kmax
    Legend{iter}=strcat('noteer ', num2str(iter),': T = ', num2str(note_periodicity{iter,kmax}(1)), '; eps = ', num2str(note_periodicity{iter,kmax}(2)), '; R = ', num2str(note_periodicity{iter,kmax}(3)));
end
legend(Legend);

title('Linear regression of t_{x,k}');
xlabel('k');
ylabel('t_{x,k}, s');
hold off