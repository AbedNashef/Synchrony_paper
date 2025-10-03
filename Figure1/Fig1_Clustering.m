clear
bin2use=1;
tbf=1; taf=.75;
win2corr=.2; d2corr= .020; ax2corr=-tbf:d2corr:taf-win2corr;
ax=-tbf*1000:bin2use:taf*1000; ax=ax(1:end-1); ax=ax./1000;
c_all=0;
c_allv=0;
c_alla=0;
n=50;
move_thresh=0.2;
iterations=1;
lags=[-.2:.01:.2];
bin=10/1000;
win=200/1000;
%%
nax= [-0.5:1/1000:0.5];nax=nax(1:end-1);
t1= find(nax<-0.25); t2= find(nax>0 & nax<.25);
load ../Data/Figure1/fr2cluster.mat;

%% new clustering score silluhette
pcn=3;
[coeff,score,latent,tsquared] = pca(all_fr');
nc=15;
avg_distance_per_c_num=1;
sem_distance_per_c_num=0;
for i=2:nc
    clusters= get_activity_clusters(all_fr,ax,1:length(ax),i);
    close all
    D_this=[];
    D_others=[];
    for ii=1:i
        ix_this= find(clusters==ii);
        X_this= coeff(ix_this,1:pcn);
        D_this= [D_this; sum(sqrt((X_this-mean(X_this,1)).^2),2)];
    end
    for ii=1:i
        for iii=ii+1:i
            if ii~=iii
                ix1=  find(clusters==ii);
                ix2=  find(clusters==iii);
                X1=mean(coeff(ix1,1:pcn),1);
                X2=mean(coeff(ix2,1:pcn),1);
                D_others=[D_others; sum(sqrt((X1-X2).^2))];
            end
        end
    end
    avg_distance_per_c_num(i)=mean((D_this./mean(D_others)));
    sem_distance_per_c_num(i)=std((D_this./mean(D_others)))./sqrt(length(D_this));
    Nshort(i)= length(find((D_this./mean(D_others))<1))/length(D_this);
    if i==4
        clustering_score=(D_this./mean(D_others));
    end
end
i=4;
clusters= get_activity_clusters(all_fr,ax,1:length(ax),i);
% close all


subplot(2,3,4)
errorbar(2:nc,1./avg_distance_per_c_num(2:end),sem_distance_per_c_num(2:end),'k','CapSize',0)
xlabel('Cluster #');
ylabel('Dist(cluster)/Dist(others)')
set(gca,'TickDir','out'); box off

subplot(2,3,5)
plot(2:nc,Nshort(2:end),'k','LineWidth',2)
xlabel('Cluster #');
ylabel('fraction <1');
set(gca,'TickDir','out'); box off


%
subplot(2,3,6)
histogram(clustering_score,'BinWidth',.1,'FaceColor','k','EdgeColor','none')
hold on
a=axis;
plot([0 0],a(3:4),'-k','LineWidth',1)
plot([mean(clustering_score) mean(clustering_score)],a(3:4),'-r','LineWidth',2)
hold off
set(gca,'TickDir','out'); box off
xlabel('clustering score')
ylabel('#')

% savefig('/results/ClusteringGoodness.fig')