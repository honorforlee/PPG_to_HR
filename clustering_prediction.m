% Ivan NY HANITRA - Master thesis
%       -- Peaks clustering and frequency prediction --

%   - Load file and data -
%Name = '3987834m';     % BPM = 78
Name = '3801060_0007m'; % BPM = 95

load(strcat(Name, '.mat'));
fid = fopen(strcat(Name, '.info'), 'rt');
fgetl(fid);
fgetl(fid);
fgetl(fid);
[interval] = sscanf(fgetl(fid), 'Sampling frequency: %f Hz  Sampling interval: %f sec');
interval = interval(2);              % data acquisition rate (interval = 1/f_spl_u = 0.5903 ms in practice)

fclose(fid);

%   - Create subdata of 5s -
range = (1 : (5/interval)) * interval;

for k = 0 : length(val) / length(range) - 1
    val_div(k+1,:)= val(  k*(length(range)) + 1 : (k+1)*(length(range)) ) ;
end

% % - Create subdata of 5s to get initial PPG frequency -
% range = (1 : (5/interval)) * interval;
% val_init(1,:)= val( 1 : length(range)) ;        % 5s data subset
% val_predict = val( length(range)+1 : end);      % remaining data

%   - Predicion and note for next data subset
freq_init = clustering_function(val_div(1,:), 8e-3, 0.1, 0.1/3, 1e-4);   % PPG average frequency of first 5s subdata

freq(1,1) = freq_init;       % PPG frequencies with initial parameters
freq_ref(1,1) = freq_init;   % reference frequency

ratio(1,1) = 1;              % rapport of frequency from one peak to next one
note_P(1,1) = 0;             % note

freq_ppg (1) = freq_init;

for j = 1 : length(val) / length(range) - 1
    dt = 0.1;
    warning (j) = 0;                                                            % count # major peaks with unexpected frequency
    [~,~, ~, ~, tx_major, ~, ~] = clustering_function(val_div(j+1,:), 8e-3, dt, dt/3, 1e-4);
    
    for k = 1 : length(tx_major) - 1
        
        freq(j,k+1) = 1./(tx_major(k+1) - tx_major(k));
        ratio(j,k+1) = freq(j,k+1)/freq_ref(j,k);
        
        freq_ref(j,k+1) = freq(j,k+1);
        note_P(j,k+1) = ratio(j,k+1);
        
        if abs(ratio(j,k+1) - 1) <= 0.15
            
        else
            warning(j) = warning(j) + 1;
        end
    end
    
    if warning > 0
        dt = dt + 0.1*dt;
        freq_ppg (j+1)= freq_ppg(j);
    else
        freq_ppg(j+1) = mean( freq(j,:) );
    end
    
    freq(j+1,1) = freq(j,length(tx_major));
    freq_ref(j+1,1) = freq(j,length(tx_major));
    ratio(j+1,1) = 1;
    note_P(j+1,1) = 0;
    
end


%%
%   - Plots -
figure (1);
plot(X(idx==1,1),X(idx==1,2),'r.','MarkerSize',12)  % cluster 1 correspoding sx,delta_note2 (minor)
hold on
plot(X(idx==2,1),X(idx==2,2),'b.','MarkerSize',12)  % cluster 2 correspoding sx,delta_note2 (major)
plot(C(:,1),C(:,2),'kx',...                         % plot centroids
    'MarkerSize',15,'LineWidth',3)

title 'Cluster Assignments and Centroids'
legend('Cluster1','Cluster2','Centroids',...
    'Location','NW')
xlabel ('Peak amplitude, a.u')
ylabel ('delta_{note2}, a.u');
hold off

figure(2);
plot(t, s,'k-','MarkerSize',12,'LineWidth',1);               % siganl s
hold on
plot(t_spl, s_spl,'ko--','MarkerSize',10,'LineWidth',1);     % sampled signal s_n
plot(td_spl, d_spl,'g--','MarkerSize',12,'LineWidth',1);     % derivative of s_n
plot(tx_major,sx_major,'rp','MarkerSize',15,'LineWidth',3);  % major peaks

title('Peaks discrimination for heart rate monitoring');
xlabel('Time, s');
ylabel('Arbitrary units');
legend('s: original signal','s_n: sampled signal','d_n: derivative of s_n','Major peaks','Location','northeastoutside');
hold off

