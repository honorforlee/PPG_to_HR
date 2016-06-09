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

% %   - Create subdata of 5s -
% range = (1 : (5/interval)) * interval;      
% 
% for k = 0 : length(val) / length(range) - 1   
%     val_div(k+1,:)= val(  k*(length(range)) + 1 : (k+1)*(length(range)) ) ;
% end

% - Create subdata of 5s to get initial PPG frequency -
range = (1 : (5/interval)) * interval; 
val_init(1,:)= val( 1 : length(range)) ;        % 5s data subset
val_predict = val( length(range)+1 : end);      % remaining data

%   - Predicion and note for next data subset -
freq_init = clustering_function(val_init, 8e-3, 0.1, 0.1/3, 1e-4);

[~,~,~, kx, tx, sx] = clustering_function(val_predict, 8e-3, 0.1, 0.1/3, 1e-4);

1./((kx(2) - kx(1)) * (tx(2) - tx(1))) 
    

% for k = 2:length(val) / length(range)
%     clustering_function(val_div(k,:), 8e-3, 0.1, 0.1/3, 1e-4);
%         
% end



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

