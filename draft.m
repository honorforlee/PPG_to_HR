%%
%X = [ note_1(:),note_2(:),note_P(:) ];  
rng default;  % For reproducibility
X = [gallery('uniformdata',[10 3],12);...
    gallery('uniformdata',[10 3],13)+1.2;...
    gallery('uniformdata',[10 3],14)+2.5];

T = clusterdata(X,'distance','cityblock','maxclust',3);
find(T==5)

scatter3(X(:,1),X(:,2),X(:,3),100,T,'filled')


%%
% tbl = table([1:length(tx)]', tx','VariableNames',{'k','tx'});
% mdl = fitlm(tbl,'tx~k');
% F_stat = anova(mdl);                        % analyse of variance
% F = F_stat.F(1);                            % F = MeanSq(xi)/MeanSq(Error) with MeanSq = SumSq/DF) (DF(xi)=1 , DF(error)=length(kx)-2)
tbl = table([1:length(kx)]', note_1', note_2', note_3',note_x','VariableNames',{'k','note_1','note_2','note_3','note_x'});

plot(tbl.k,tbl.note_1,'.r'...
    , tbl.k,tbl.note_2,'.b'...
    , tbl.k,tbl.note_3,'.g'...
    , tbl.k,tbl.note_x,'xk','MarkerSize',12 ...
    );
legend('note_1','note_2','note_3','note_x');
