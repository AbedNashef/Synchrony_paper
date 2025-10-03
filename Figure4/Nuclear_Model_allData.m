% model
clear
C=89; vr=-60; vt=-40;
k=1; % not relevant here
a=.03; b=-2; c=-50; d=2;
gmax=5; %nS
vpeak= 35; % spike cutoff
Er=-75;
ex_input=800;
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
%%
load FR4allClusters.mat;
load inter_unit_xcorr_dist.mat;
load ccg4model.mat;
nPC=40;
this_xcorr=all_ccg(find(sum(all_clusters,2)==3),find(xcorr_ax==0));
% this_xcorr=this_xcorr(find(this_xcorr>median(this_xcorr)));
% this_xcorr=xcorr_dist{1};
% this_xcorr=this_xcorr(find(this_xcorr>0 & this_xcorr<20));
this_xcorr=this_xcorr(this_xcorr>=quantile(this_xcorr,.8));
% this_xcorr(find(this_xcorr<0))=0;

xcorr_1_1=4; %4 sp/s cross-corr
n2sync= round((length(ax)/1000)*xcorr_1_1);
iterations=1000;
jitter=1;%:size(all_perc,2);
t=(ax(end)-ax(1))*1000;

ccg_observed=[];
avg_ccg_observed=[];
dNuc=[];
%%
iter=0;
while iter<iterations
    c1_ix= find([all_FR.cluster]<=2);

    %
    [~,rix]= sort(rand(1,length(c1_ix)));
    rix= c1_ix(rix(1:nPC));
    for i=1:nPC
        nullFR(i,:)= mean(all_FR(rix(i)).null,1);
    end
    sim_FR=poissrnd(nullFR);
    sim_FR(find(sim_FR>1))=1;
    no_additional_sync=sim_FR;
    L1=sum(sim_FR(:,find(ax<-.5)),2);
    all_n2sync=nan(nPC,nPC);
    [~,rix]=sort(rand(1,nPC));
    pairs=reshape(rix(:),[],2);
    for cc=1:size(pairs,1)
        i=pairs(cc,1);
        % for i=1:size(sim_FR,1)
        % for ii=i+1:size(sim_FR,1)
        ii=pairs(cc,2);
        xcorr_1_1=this_xcorr(randi(length(this_xcorr),[1 1]))*1000; %4 sp/s cross-corr
        n2sync= round((length(ax)/1000)*xcorr_1_1);
        all_n2sync(i,ii)=n2sync;
        if n2sync>0 & xcorr_1_1>0
            ix1= find(sim_FR(i,:)==1 & sim_FR(ii,:)==0);
            ix2= find(sim_FR(ii,:)==1);
            try
                [~,tmp]=sort(rand(1,length(ix1))); ix1=ix1(tmp(1:n2sync));
                [~,tmp]=sort(rand(1,length(ix2))); ix2=ix2(tmp(1:n2sync));
                sim_FR(ii,ix2)=0;
                sim_FR(ii,ix1)=1;
            catch me
                n2sync=min([length(ix1) length(ix2)]);
                [~,tmp]=sort(rand(1,length(ix1))); ix1=ix1(tmp(1:n2sync));
                [~,tmp]=sort(rand(1,length(ix2))); ix2=ix2(tmp(1:n2sync));
                sim_FR(ii,ix2)=0;
                sim_FR(ii,ix1)=1;
            end
        else
            ix1= find(sim_FR(i,:)==1 & sim_FR(ii,:)==1);
            ix0= find(sim_FR(i,:)==0 & sim_FR(ii,:)==0);
            [~,rix]=sort(rand(1,length(ix0)));
            rix=rix(1:min([round(abs(n2sync)) length(ix1)]));
            [~,rix1]=sort(rand(1,length(ix1)));
            rix1=rix1(1:min([round(abs(n2sync)) length(ix1)]));
            sim_FR(ii,ix1(rix1))=0;
            sim_FR(ii,ix0(rix))=1;
        end
        % end

    end

    sim_FR(find(sim_FR==2))=1;
    L2=sum(sim_FR(:,find(ax<-.5)),2);

    [~,p]= ttest(L1,L2);
    if p>.05
        iter=iter+1;
    else
        continue
    end
    Pop_sync_noEnforced(iter,:)=nanmean(no_additional_sync,1);

    % if L2<L1
    %     index2add=randi
    % end
    Pop_sync_Enforced(iter,:)=nanmean(sim_FR,1);
    added_sync=sim_FR;


    PC_trc=no_additional_sync;
    % PC_trc=zeros(size(no_additional_sync));
    v=vr*ones(1,t+1); u=0*v; % initial values
    I=ones(1,t+1).*ex_input;
    for i=1:nPC
        ix= find(PC_trc(i,:));
        G_all=zeros(length(ix),t+1);
        for j=1:length(ix)
            for ii=ix(j):ix(j)+n
                tpeak=ix(j)+(tau*tau_rise/(tau-tau_rise))*log(tau/tau_rise);
                f=1/(-1*exp(-1*(tpeak-ix(j))/tau_rise)+exp(-1*(tpeak-ix(j))/tau));
                G_all(j,ii)=gmax*f*(exp(-1*(ii-ix(j))/tau)-exp(-1*((ii-ix(j))/tau_rise)));
            end
        end
        G_all=G_all(:,1:length(ax));
        G_this(i,:)= sum(G_all,1);
        %             G_this(i,:)= (G_this(i,:)./max(G_this(i,:))).*gmax;
    end
    for i=1:size(PC_trc,2)
        for ii=1:size(G_this,1)
            I(i)= I(i) + Cm*G_this(ii,i)*(Er-v(i));
        end
        v(i+1) = v(i)+tau*(k*(v(i)-vr)*(v(i)-vt)-u(i)+I(i))/C;
        u(i+1) = u(i)+tau*a*(b*(v(i)-vr)-u(i));
        if v(i+1)>=vpeak
            v(i)=vpeak; % a spike is fired
            nFR(i)=1;
            v(i+1)=c; % membrane voltage reset
            u(i+1)=u(i+1)+d; % recovery variable update
        else
            nFR(i)=0;
        end
    end

    OVERALL_NUC_G(iter,1)=nanmean(sum(G_this,2));
    OVERALL_NUC_V(iter,1)=nanmean(v(:));
    OVERALL_NUC_FR(iter,1)=sum(nFR)/(length(ax)/1000);
    OVERALL_PC_FR(iter,1)=mean(sum(added_sync,2))/(length(ax)/1000);
    fr1=sum(nFR(find(ax>-.3 & ax<.2)));
    fr0=sum(nFR(find(ax<-.5)));
    modulation_index(iter,1)= (fr1-fr0)/(fr1+fr0);
    Nuclear_FR(iter,1,:)=nFR;
    all_VOLTAGES(iter,1,:)=v;

    clear G_this v nFR;

    PC_trc=added_sync;
    v=vr*ones(1,t+1); u=0*v; % initial values
    I=ones(1,t+1).*ex_input;
    for i=1:nPC
        ix= find(PC_trc(i,:));
        G_all=zeros(length(ix),t+1);
        for j=1:length(ix)
            for ii=ix(j):ix(j)+n
                tpeak=ix(j)+(tau*tau_rise/(tau-tau_rise))*log(tau/tau_rise);
                f=1/(-1*exp(-1*(tpeak-ix(j))/tau_rise)+exp(-1*(tpeak-ix(j))/tau));
                G_all(j,ii)=gmax*f*(exp(-1*(ii-ix(j))/tau)-exp(-1*((ii-ix(j))/tau_rise)));
            end
        end
        G_all=G_all(:,1:length(ax));
        G_this(i,:)= sum(G_all,1);
        %             G_this(i,:)= (G_this(i,:)./max(G_this(i,:))).*gmax;
    end
    for i=1:size(PC_trc,2)
        for ii=1:size(G_this,1)
            I(i)= I(i) + Cm*G_this(ii,i)*(Er-v(i));
        end
        v(i+1) = v(i)+tau*(k*(v(i)-vr)*(v(i)-vt)-u(i)+I(i))/C;
        u(i+1) = u(i)+tau*a*(b*(v(i)-vr)-u(i));
        if v(i+1)>=vpeak
            v(i)=vpeak; % a spike is fired
            nFR(i)=1;
            v(i+1)=c; % membrane voltage reset
            u(i+1)=u(i+1)+d; % recovery variable update
        else
            nFR(i)=0;
        end
    end

    OVERALL_NUC_G(iter,2)=nanmean(sum(G_this,2));
    OVERALL_NUC_V(iter,2)=nanmean(v(:));
    OVERALL_NUC_FR(iter,2)=sum(nFR)/(length(ax)/1000);
    OVERALL_PC_FR(iter,2)=mean(sum(added_sync,2))/(length(ax)/1000);

    if OVERALL_NUC_FR(iter,2)-OVERALL_NUC_FR(iter,1)>5
        ex_raster(end+1,1,:,:)=no_additional_sync;
        ex_raster(end,2,:,:)=added_sync;
        dNuc(end+1,:)=[OVERALL_NUC_FR(iter,1) OVERALL_NUC_FR(iter,2)];
    end

    fr1=sum(nFR(find(ax>-.3 & ax<.2)));
    fr0=sum(nFR(find(ax<-.5)));
    modulation_index(iter,2)= (fr1-fr0)/(fr1+fr0);

    PC_FR(iter,1,:)=nanmean(no_additional_sync,1);
    PC_FR(iter,2,:)=nanmean(added_sync,1);

    PC_ISI(iter,1)=nanmean(diff(find(mean(no_additional_sync,1)>0)));
    PC_ISI(iter,2)=nanmean(diff(find(mean(added_sync,1)>0)));
    Nuclear_FR(iter,2,:)=nFR;
    all_VOLTAGES(iter,2,:)=v;

    ccg_observed=[ccg_observed; get_ccg(no_additional_sync) get_ccg(added_sync)];
    avg_ccg_observed=[avg_ccg_observed; mean(get_ccg(no_additional_sync)) mean(get_ccg(added_sync))];
    % all_ccg(iter,1)=get_ccg(no_additional_sync);
    % all_ccg(iter,2)=get_ccg(added_sync);
    all_ccg(iter)=nanmean(all_n2sync(:))/(length(ax)/1000);


    clear PC_trc nFR G_all G_this v XX
end
%%
ex2plot=3;
[~,rix]=sort(rand(1,size(ex_raster,1)));
[~,rix]=sort(diff(dNuc,[],2),'descend');
rix=rix(1:min([length(rix) 3]));% examples
for i=1:length(rix)
    subplot(ex2plot*3,4,(i-1)*12+1)
    plot(squeeze(all_VOLTAGES(rix(i),1,:)))
    set(gca,'TickDir','out'); box off
    subplot(ex2plot*3,4,[(i-1)*12+5 (i-1)*12+9])
    plot_raster(squeeze(ex_raster(rix(i),1,:,:)),ax)
    yyaxis right
    plot(ax,smooth(nanmean(squeeze(ex_raster(rix(i),1,:,:)),1)*1000,21))

    subplot(ex2plot*3,4,(i-1)*12+2)
    plot(squeeze(all_VOLTAGES(rix(i),2,:)))
    set(gca,'TickDir','out'); box off
    subplot(ex2plot*3,4,[(i-1)*12+6 (i-1)*12+10])
    plot_raster(squeeze(ex_raster(rix(i),2,:,:)),ax)
    yyaxis right
    plot(ax,smooth(nanmean(squeeze(ex_raster(rix(i),2,:,:)),1)*1000,21))

    title(diff(dNuc(rix(i),:)))
end
%%
subplot(4,4,3)
plot([1 2],OVERALL_PC_FR,'k','LineWidth',1)
hold on
errorbar([.95 2.05],nanmean(OVERALL_PC_FR,1),nanstd(OVERALL_PC_FR,[],1),'color','r','LineWidth',2,'CapSize',0);
hold off
set(gca,'TickDir','out'); box off
set(gca,'XTick',[1 2],'XTickLabel',{'No added sync' 'added sync'})
ylabel('Avg. PC firing')

subplot(4,4,7)
plot([1 2],PC_ISI,'k','LineWidth',1)
hold on
errorbar([.95 2.05],nanmean(PC_ISI,1),nanstd(PC_ISI,[],1),'color','r','LineWidth',2,'CapSize',0);
hold off
set(gca,'TickDir','out'); box off
set(gca,'XTick',[1 2],'XTickLabel',{'No added sync' 'added sync'})
ylabel('Avg ISI (ms)')


subplot(4,4,11)
plot([1 2],OVERALL_NUC_V,'k','LineWidth',1)
hold on
errorbar([.95 2.05],nanmean(OVERALL_NUC_V,1),nanstd(OVERALL_NUC_V,[],1),'color','r','LineWidth',2,'CapSize',0);
hold off
set(gca,'TickDir','out'); box off
set(gca,'XTick',[1 2],'XTickLabel',{'No added sync' 'added sync'})
ylabel('Avg. nuclear voltage')

subplot(4,4,15)
plot([1 2],OVERALL_NUC_FR,'k','LineWidth',1)
hold on
errorbar([.95 2.05],nanmean(OVERALL_NUC_FR,1),nanstd(OVERALL_NUC_FR,[],1),'color','r','LineWidth',2,'CapSize',0);
hold off
set(gca,'TickDir','out'); box off
set(gca,'XTick',[1 2],'XTickLabel',{'No added sync' 'added sync'})
ylabel('Avg. nuclear firing')


tmp_sync=squeeze(mean(PC_FR(:,2,1:500),3));
tmp_nosync=squeeze(mean(PC_FR(:,1,1:500),3));
stmp_nosync=squeeze(std(PC_FR(:,1,1:500),[],3))./sqrt(500);
ix = find(tmp_sync>=tmp_nosync-stmp_nosync & tmp_sync<=tmp_nosync+stmp_nosync);


subplot(5,4,4)
plot(ax,squeeze(nanmean(PC_FR,1)).*1000)
hold on
plot(ax,smooth(squeeze(nanmean(PC_FR(ix,2,:),1)).*1000,21),'k')
set(gca,'TickDir','out'); box off
xlabel('Time to endpoint (s)')
ylabel('PC firing rate (sp/s)')

subplot(5,4,8)
plot([ax ax(end)+1/1000],squeeze(nanmean(all_VOLTAGES(:,1,:),1)))
hold on
plot([ax ax(end)+1/1000],squeeze(nanmean(all_VOLTAGES(:,2,:),1)))
hold off
set(gca,'TickDir','out'); box off
xlabel('Time to endpoint (s)')
ylabel('Voltage')

subplot(5,4,12)
histogram(avg_ccg_observed(:,1),'BinWidth',1e-4,'EdgeColor','none','FaceColor','k')
hold on
histogram(avg_ccg_observed(:,2),'BinWidth',1e-4,'EdgeColor','none','FaceColor','g')
hold off
set(gca,'TickDir','out'); box off
xlabel('CCG')
ylabel('#')

subplot(5,4,16)
scatter(diff(avg_ccg_observed,[],2),diff(OVERALL_NUC_FR,[],2),'.k')
lsline
[r,p]=corr(diff(avg_ccg_observed,[],2),diff(OVERALL_NUC_FR,[],2));
title([num2str(r) '; ' num2str(p)])
xlabel('CCG')
ylabel('\DeltaFR')

subplot(5,4,20)
histogram(modulation_index(:,1),'BinWidth',.05,'FaceColor','b')
hold on
histogram(modulation_index(:,2),'BinWidth',.05,'FaceColor','r')
hold off
xlabel('modulation index')
set(gca,'TickDir','out'); box off

% subplot(5,4,20)
% plot([1 2],PC_ISI,'-k','LineWidth',1)
% hold on
% errorbar([.95 2.05],nanmean(PC_ISI,1),nanstd(PC_ISI,[],1),'color','r','LineWidth',2,'CapSize',0);
% hold off
% set(gca,'TickDir','out'); box off
% set(gca,'XTick',[1 2],'XTickLabel',{'No added sync' 'added sync'})
% ylabel('Inter-IPSC interval')

%%
gwin=gausswin(21); gwin= gwin./sum(gwin);
for i=1:size(Nuclear_FR,1)
    for ii=1:size(Nuclear_FR,2)
        this_fr=squeeze(Nuclear_FR(i,ii,:));
        dix= diff(find(this_fr==1));
        insta_fr(i,ii,:)=conv2(this_fr,gwin,'same');
        all_csum(i,ii,:)=cumsum(conv2(this_fr,gwin,'same'))./sum(conv2(this_fr,gwin,'same'));
        cv(i,ii)=std(dix)/mean(dix);
        dix=[dix(:) circshift(dix(:),1,1)];
        cv2(i,ii)=sum((abs(diff(dix,[],2))./sum(dix,2)))/size(dix,1);
    end
end

%%
gwin=gausswin(201); gwin= gwin./sum(gwin);
for i=1:size(Nuclear_FR,1)
    for ii=1:size(Nuclear_FR,2)
        x=squeeze(Nuclear_FR(i,ii,:));
        cx=conv2(x,gwin,'same');
        [peakx,ix]= max(cx);
        ix1= find(cx(1:ix)<peakx/2,1,'last');
        try
            edge1=interp1(cx(ix1:ix),ax(ix1:ix),peakx/2);
        catch me
            edge1=nan;
        end
        ix2= ix+find(cx(ix:end)<peakx/2,1,'first');
        try
            edge2=interp1(cx(ix:ix2),ax(ix:ix2),peakx/2);
        catch me
            edge2=nan;
        end
        HWHH(i,ii)=edge2-edge1;
    end
end
subplot(5,4,20)
plot([1 2],HWHH,'-k','LineWidth',1)
hold on
errorbar([.95 2.05],nanmean(HWHH,1),nanstd(HWHH,[],1),'color','r','LineWidth',2,'CapSize',0);
hold off
[~,p]= ttest(HWHH(:,1),HWHH(:,2));
title(p)
ylabel('HWHH (s)')
set(gca,'TickDir','out'); box off
set(gca,'XTick',[1 2],'XTickLabel',{'stochastic' 'physiological'})
%% For paper
subplot(2,2,1)
plot([1 2],OVERALL_PC_FR,'k','LineWidth',1)
hold on
errorbar([.9 2.1],nanmean(OVERALL_PC_FR,1),nanstd(OVERALL_PC_FR,1),'r','CapSize',0)
hold off
set(gca,'TickDir','out'); box off
ylabel('PC FR (sp/s)')

subplot(2,2,2)
plot([1 2],PC_ISI,'k','LineWidth',1)
hold on
errorbar([.9 2.1],nanmean(PC_ISI,1),nanstd(PC_ISI,1),'r','CapSize',0)
hold off
set(gca,'TickDir','out'); box off
ylabel('Average ISI (ms)')

