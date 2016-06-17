%X = [ note_1(:),note_2(:),note_P(:) ];  
rng default;  % For reproducibility
X = [gallery('uniformdata',[10 3],12);...
    gallery('uniformdata',[10 3],13)+1.2;...
    gallery('uniformdata',[10 3],14)+2.5];

T = clusterdata(X,'distance','cityblock','maxclust',3);
find(T==5)

scatter3(X(:,1),X(:,2),X(:,3),100,T,'filled')
