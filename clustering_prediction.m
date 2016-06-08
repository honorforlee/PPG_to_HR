% Ivan NY HANITRA - Master thesis
%       -- Discrimination with sampled signal --

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

%   - Create subdata -

range = (1 : (5/interval)) * interval;
L = length(range);

%val_div = zeros(12,length(range));

val_div(1,:) = val(1:length(range));


for k = 1 : length(val) / length(range) - 1   
    val_div(k+1,:)= val(  k*(length(range)) + 1 : (k+1)*(length(range)) ) ;
    
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
xlabel ('sx')
ylabel ('delta_{note2}');
hold off

figure(2);
plot(t, s,'k-'...               % siganl s
    ,t_spl, s_spl,'ko--'...     % sampled signal s_n
    ,td_spl, d_spl,'g--'...     % derivative of s_n 
    );
hold on

plot(tx_major_c,sx_major_c,'rp','MarkerSize',15,'LineWidth',3);

title('Peaks discrimination for heart rate monitoring');
xlabel('Time, s');
ylabel('Arbitrary units');
legend('s: original signal'...
    ,'s_n: sampled signal'...
    ,'derivative of s_n'...
    ,'Location','northeastoutside');
hold off
