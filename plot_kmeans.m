% X = [ note_1(:),note_2(:),delta(:) ];                        % data
% 
% [idx,C] = kmeans(X,2,'Distance','cityblock',...     % 2 clusters created: minor/major peaks
%     'Replicates',5,'Start','plus','Options',statset('Display','final'));  % initialize the replicates 5 times, separately using k-means++ algorithm, choose best arrangement and display final output
% 
% plot3(X(idx==1,1),X(idx==1,2),X(idx==1,3),'r.','MarkerSize',15);  % cluster 1 correspoding sx,delta_note2 (minor)
% hold on
% plot3(X(idx==2,1),X(idx==2,2),X(idx==2,3),'b.','MarkerSize',15);  % cluster 2 correspoding sx,delta_note2 (major)
% plot3(C(:,1),C(:,2),C(:,3),'kx','MarkerSize',20,'LineWidth',3);
% 
% title ('Cluster Assignments and Centroids','FontSize',30,'FontWeight','Bold');
% legend('Cluster1','Cluster2','Centroids','Location','NW');
% xlabel('note_1','FontSize',25,'FontWeight','Bold');
% ylabel('note_2','FontSize',25,'FontWeight','Bold');
% zlabel('note_3','FontSize',25,'FontWeight','Bold');
% grid on
% hold off

X = [note_x(:)];
[idx,C] = kmeans(X,2,'Distance','cityblock',...     % 2 clusters created: minor/major peaks
     'Replicates',5,'Start','plus','Options',statset('Display','final'));  % initialize the replicates 5 times, separately using k-means++ algorithm, choose best arrangement and display final output 

hold on
xlim([-1,5]);ylim([-.2,2])

%arrow
[arrowX,arrowY]=dsxy2figxy([-1,5],[0,0]);
annotation('arrow',arrowX,arrowY)

%crosses

plot(X(idx==1),0,'or','markersize',10);
hold on
plot(X(idx==2),0,'ob','markersize',10);
plot(C,0,'kx','MarkerSize',20,'LineWidth',3);
xlabel ('note_x','FontSize',25,'FontWeight','Bold');

%pipes
p=[0.5,0.65];
text(p,[0,0],'$$\vert$$','interpreter','latex')

grid on
axis on
print('-depsc','arrowFigure')
%%
cluster1 = find(idx==1)';
cluster2 = find(idx==2)';

if sx(cluster1(1)) > sx(cluster2(1))                % assign major peak cluster
    major_index = cluster1;
else
    major_index = cluster2;
end

tx_major = tx(major_index);
sx_major = sx(major_index);
plot(t,s,'o:k');
hold on
plot(tx_major,sx_major,'pr','MarkerSize',15);
 
%%

plot(X(idx==1,1),X(idx==1,2),X(idx==1,3),'r.','MarkerSize',15)  % cluster 1 correspoding sx,delta_note2 (minor)
hold on
plot(X(idx==2,1),X(idx==2,2),X(idx==2,3),'b.','MarkerSize',15)  % cluster 2 correspoding sx,delta_note2 (major)
plot(C(:,1),C(:,2),'kx',...                         % plot centroids
    'MarkerSize',20,'LineWidth',3)

title ('Cluster Assignments and Centroids','FontSize',30,'FontWeight','Bold');
legend('Cluster1','Cluster2','Centroids','Location','NW');
xlabel ('note_1','FontSize',25,'FontWeight','Bold');
ylabel ('note_3','FontSize',25,'FontWeight','Bold');
hold off

cluster1 = find(idx==1)';
cluster2 = find(idx==2)';

if sx(cluster1(1)) > sx(cluster2(1))                % assign major peak cluster
    major_index = cluster1;
else
    major_index = cluster2;
end

for k = 1:length(major_index)-1
    
    freq = mean( major_index(k+1) - major_index);
    
end
%%
X = [ note_1(:),note_2(:),delta(:) ];                        % data

c = clusterdata(X,'linkage','ward','savememory','on','maxclust',2);
scatter3(X(:,1),X(:,2),X(:,3),36,c);


%%
plot(t0,s0,'-k');
hold on
plot(t,s,'ok');
plot(tx,sx,'pb');