% model
clear
addpath C:\Users\PersonLab\Desktop\'synchrony paper'\Dylans_data\helper_functions\
C=89; vr=-60; vt=-40;
k=0.056; % k>0.056 and k<.2 gives ~90sp/s when no inhibition injected (Raman 2000)
a=.03; b=-2; c=-50; d=2;
gmax=5; %nS
vpeak= 35; % spike cutoff
Er=-75;
ex_input=600;
T=100;%s preceded by a brief increase in activity, presumably due to excitation of the DCN neuron by col0;
dt=1; % time span and strp (ms)
n=round(T/dt); % number of simulation steps
tau=2.5; % based on Person and Raman, 2012
tau_rise=.05;
Cm =1; %better as tau
Gsyn=0;
rp=2;
vplot=zeros(2,2);
win=50;
c_all=0;
ex_raster=[];
addpath ../../Dylans_data/helper_functions/
ntr=15;
nbin=1;
%%

load C:\Users\PersonLab\'OneDrive - The University of Colorado Denver'\Desktop\'synchrony paper'\Dylans_data\HH_model\FR4allClusters.mat;
load C:\Users\PersonLab\'OneDrive - The University of Colorado Denver'\Desktop\'synchrony paper'\Dylans_data\HH_model\inter_unit_xcorr_dist.mat;
load C:\Users\PersonLab\'OneDrive - The University of Colorado Denver'\Desktop\'synchrony paper'\Dylans_data\HH_model\ccg4model.mat;
load E:\Dylans_data\jPSTHs\Cluster_jPETHs\jPSTHs2model.mat;
nPC=40;
this_xcorr=all_ccg(find(sum(all_clusters,2)==3),find(xcorr_ax==0));
% this_xcorr=this_xcorr(find(this_xcorr>median(this_xcorr)));
% this_xcorr=xcorr_dist{1};
% this_xcorr=this_xcorr(find(this_xcorr>0 & this_xcorr<20));
this_xcorr=this_xcorr(this_xcorr>=quantile(this_xcorr,.8));
% this_xcorr(find(this_xcorr<0))=0;

xcorr_1_1=4; %4 sp/s cross-corr
n2sync= round((length(ax)/1000)*xcorr_1_1);
iterations=500;
jitter=1;%:size(all_perc,2);
t=(ax(end)-ax(1));

ccg_observed=[];
avg_ccg_observed=[];
dNuc=[];
%% high Sync pairs
index= find(cov_bfr_move(:,1)<=quantile(cov_bfr_move(:,1),1));
% index= find(cov_bfr_move(:,1)>quantile(cov_bfr_move(:,1),.75));
iter=0;
ccg_observed=[];
avg_ccg_observed=[];
dNuc=[];
while iter<iterations
    %
    [~,rix]= sort(rand(1,length(index)));
    rix=index(rix(1:nPC/2));
    clear MAT;
    try
    for i=1:length(rix)
        [~,tr_ix]= sort(rand(size(all_FR{rix(i),1},1),1));
        tr_ix=tr_ix(1:ntr);
        MAT(i,1,:,:)=all_FR{rix(i),1}(tr_ix,:);
        MAT(i,2,:,:)=all_FR{rix(i),2}(tr_ix,:);
    end
    added_sync=MAT;
    %
    iter=iter+1;
    catch me
        continue
    end
   
    for i=1:size(added_sync,1)
        cell1=squeeze(added_sync(i,1,:,:));
        cell2=squeeze(added_sync(i,2,:,:));
        [raw, shift_predict, pred, surprise, std, ~, ~] = my_JPSTH(cell1', cell2', nbin,1);
        ijpsth_syc(i,:,:)=raw-pred;
        
        cell1=squeeze(added_sync(i,1,:,find(ax>-250 & ax<100)));
        cell2=squeeze(added_sync(i,2,:,find(ax>-250 & ax<100)));

        for ii=1:size(cell1,1)
            ccg(ii,:)=get_ccg_whole([cell1(ii,:); cell2(ii,:)]);
        end
        this_ccg(i,:)=nanmean(ccg,1);

    end
    all_jpsth(iter,:,:)=squeeze(nanmean(ijpsth_syc,1));
    all_CCG(iter,:)=nanmean(this_ccg,1);
    for i=1:size(ijpsth_syc,1)
        mat=imgaussfilt(squeeze(ijpsth_syc(i,find(ax>-250 & ax<100),find(ax>-250 & ax<100))),1);
        all_diag(i,:)=diag(mat);
    end
    all_dMeans(:,iter)=nanmean(all_diag,2);
    all_CCGpeaks(:,iter)=this_ccg(:,ceil(size(this_ccg,2)/2));
end
%%



nbin=11;
X=all_dMeans(:)./ntr;
Y=all_CCGpeaks(:);
index= find(~isoutlier(X) & ~isoutlier(Y));
X=X(index);
Y=Y(index);
xmin=quantile(X,.025);
xmax=quantile(X,.975);
% xmin=-0.9975;  xmax= 67.4927;
dax=linspace(xmin,xmax,nbin+1);
d=median(diff(dax));
avg_ccg=[];std_ccg=[];
for i=1:length(dax)
    ix= find(X>=dax(i) & X<=dax(i)+d);
    avg_ccg(i)=nanmean(Y(ix));
    std_ccg(i)=nanstd(Y(ix))./sqrt(length(ix));
end

subplot(2,2,1)
scatter(X,Y,'.k')
lsline
set(gca,'TickDir','out'); box off
xlabel('1-ms Cofiring sp2/s2)')
ylabel('0-lag CCG')

subplot(2,2,2)
errorbar(dax,avg_ccg,std_ccg)
hold on
scatter(dax,avg_ccg);
lsline
hold off
h=fitlm(dax,avg_ccg);
title(['R2=' num2str(h.Rsquared.Adjusted)]);
set(gca,'TickDir','out'); box off
xlabel('1-ms Cofiring sp2/s2)')
ylabel('0-lag CCG')

%%
for i=1:size(all_jpsth,1)
    mat=squeeze(all_jpsth(i,:,:));
    all_da_diag(i,:)=diag(mat);
end
subplot(2,2,3)
scatter(nanmean(all_da_diag(:,find(ax>-250 & ax<100)),2)./ntr,all_CCG(:,ceil(size(this_ccg,2)/2)),30,'ok')
lsline
x=nanmean(all_da_diag(:,find(ax>-250 & ax<100)),2)./ntr;
y=all_CCG(:,ceil(size(this_ccg,2)/2));
set(gca,'TickDir','out'); box off
xlabel('1-ms Cofiring sp2/s2)')
ylabel('0-lag CCG')


nbin=7;
xmin=quantile(x,.05);
xmax=quantile(x,.95);
dax=linspace(xmin,xmax,nbin+1);
d=median(diff(dax));
avg_ccg=[];std_ccg=[];
for i=1:length(dax)
    ix= find(x>=dax(i) & x<=dax(i)+d);
    avg_ccg(i)=nanmean(y(ix));
    std_ccg(i)=nanstd(y(ix))./sqrt(length(ix));
end
subplot(2,2,4)
errorbar(dax,avg_ccg,std_ccg)
set(gca,'TickDir','out'); box off
xlabel('1-ms Cofiring sp2/s2)')
ylabel('0-lag CCG')

%%
y=all_CCG(:,ceil(size(this_ccg,2)/2));
d=quantile(y,0:.25:1);
for i=1:length(d)-1
    ix= find(y>=d(i) & y<=d(i+1));
    avg_diag_quant(i,:)=smooth(nanmean(all_da_diag(ix,:),1),11)./ntr;
end
subplot(2,2,4)
plot(ax,avg_diag_quant,'LineWidth',2)
xlim([-250 100])
set(gca,'TickDir','out'); box off
xlabel('Time (ms)')
ylabel('cofiring (sp2/s2)')