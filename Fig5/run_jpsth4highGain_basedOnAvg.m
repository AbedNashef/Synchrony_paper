% model
clear
C=89; vr=-60; vt=-40;
k=0.056;%1; % not relevant here
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
%%
load C:\Users\PersonLab\Desktop\'synchrony paper'\Dylans_data\HH_model\FR4allClusters.mat;
load C:\Users\PersonLab\Desktop\'synchrony paper'\Dylans_data\HH_model\inter_unit_xcorr_dist.mat;
load C:\Users\PersonLab\Desktop\'synchrony paper'\Dylans_data\HH_model\ccg4model.mat;
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
iterations=1;
jitter=1;%:size(all_perc,2);
t=(ax(end)-ax(1));

ccg_observed=[];
avg_ccg_observed=[];
dNuc=[];
all_dt=[];
maxshuffle=1:10:201;
%% Low gain Sync
index= find(cov_bfr_move(:,1)>60.7-10.7*3 & cov_bfr_move(:,1)<=60.7+10.7*3);

for ss=1:length(maxshuffle)
    ccg_observed=[];
    avg_ccg_observed=[];
    dNuc=[];
    iter=0;
    while iter<iterations
        %
        [~,rix]= sort(rand(1,length(index)));
        rix=index(rix(1:nPC/2));
        MAT=[];
        for i=1:length(rix)
            tr_ix= randi(size(all_FR{rix(i),1},1));
            MAT=[MAT; all_FR{rix(i),1}(tr_ix,:)];
            MAT=[MAT; all_FR{rix(i),2}(tr_ix,:)];
        end
        real_sync=MAT;

        newMAT=zeros(size(MAT));
        for i=1:size(MAT,1)
            ix= find(MAT(i,:)==1);
            tmp1=randn(1,length(ix));
            dtsign=(tmp1./abs(tmp1)); % shuffle of 1ms!
            dt=randi(maxshuffle(ss),1,length(ix)).*dtsign;
            rix= ix+dt;
            ix= find(rix<1 | rix>size(MAT,2));
            rix(ix)=randi(size(MAT,2),1,length(ix));
            rix=unique(rix);
            newMAT(i,rix)=1;
            spk2add=length(find(MAT(i,:)==1))-length(rix);
            ix= find(newMAT(i,:)==0);
            [~,tmp]= sort(rand(1,length(ix)));
            newMAT(i,ix(tmp(1:spk2add)))=1;
            all_dt=[all_dt; dt(:)];
        end
        desync_mat=newMAT;

        iter=iter+1;

        PC_trc=real_sync;
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
            if v(i+1)>=vt
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
        OVERALL_PC_FR(iter,1)=mean(sum(real_sync,2))/(length(ax)/1000);

        fr1=sum(nFR(find(ax>-.3 & ax<.2)));
        fr0=sum(nFR(find(ax<-.5)));
        modulation_index(iter)= (fr1-fr0)/(fr1+fr0);

        PC_FR(iter,1,:)=nanmean(real_sync,1);

        PC_ISI(iter,1)=nanmean(diff(find(mean(real_sync,1)>0)));
        Nuclear_FR(iter,1,:)=nFR;
        all_VOLTAGES(iter,1,:)=v;



        PC_trc=desync_mat;
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
            if v(i+1)>=vt
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
        OVERALL_PC_FR(iter,2)=mean(sum(desync_mat,2))/(length(ax)/1000);

        PC_FR(iter,2,:)=nanmean(desync_mat,1);
        PC_ISI(iter,2)=nanmean(diff(find(mean(desync_mat,1)>0)));
        Nuclear_FR(iter,2,:)=nFR;
        all_VOLTAGES(iter,2,:)=v;


        % ccg_observed=[ccg_observed; get_ccg(real_sync)];
        avg_ccg_observed=[avg_ccg_observed; mean(get_ccg(real_sync)) mean(get_ccg(desync_mat))];
        % all_ccg(iter,1)=get_ccg(no_additional_sync);
        % all_ccg(iter,2)=get_ccg(real_sync);

        clear PC_trc nFR G_all G_this v XX
    end
    all_CbNFR_low(ss,:,:)=OVERALL_NUC_FR;
    all_CCG_low(ss,:,:)=avg_ccg_observed;
    all_ISI_low(ss,:,:)=PC_ISI;
    clear OVERALL_NUC_FR avg_ccg_observed PC_ISI
end
%% high gain sync
index= find(cov_bfr_move(:,1)>123-10*3 & cov_bfr_move(:,1)<=123+10*3);

for ss=1:length(maxshuffle)
    ccg_observed=[];
    avg_ccg_observed=[];
    dNuc=[];
    iter=0;
    while iter<iterations
        %
        [~,rix]= sort(rand(1,length(index)));
        rix=index(rix(1:nPC/2));
        MAT=[];
        for i=1:length(rix)
            tr_ix= randi(size(all_FR{rix(i),1},1));
            MAT=[MAT; all_FR{rix(i),1}(tr_ix,:)];
            MAT=[MAT; all_FR{rix(i),2}(tr_ix,:)];
        end
        real_sync=MAT;

        newMAT=zeros(size(MAT));
        for i=1:size(MAT,1)
            ix= find(MAT(i,:)==1);
            tmp1=randn(1,length(ix));
            dtsign=(tmp1./abs(tmp1)); % shuffle of 1ms!
            dt=randi(maxshuffle(ss),1,length(ix)).*dtsign;
            rix= ix+dt;
            ix= find(rix<1 | rix>size(MAT,2));
            rix(ix)=randi(size(MAT,2),1,length(ix));
            rix=unique(rix);
            newMAT(i,rix)=1;
            spk2add=length(find(MAT(i,:)==1))-length(rix);
            ix= find(newMAT(i,:)==0);
            [~,tmp]= sort(rand(1,length(ix)));
            newMAT(i,ix(tmp(1:spk2add)))=1;
            all_dt=[all_dt; dt(:)];
        end
        desync_mat=newMAT;

        iter=iter+1;

        PC_trc=real_sync;
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
            if v(i+1)>=vt
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
        OVERALL_PC_FR(iter,1)=mean(sum(real_sync,2))/(length(ax)/1000);

        fr1=sum(nFR(find(ax>-.3 & ax<.2)));
        fr0=sum(nFR(find(ax<-.5)));
        modulation_index(iter)= (fr1-fr0)/(fr1+fr0);

        PC_FR(iter,1,:)=nanmean(real_sync,1);

        PC_ISI(iter,1)=nanmean(diff(find(mean(real_sync,1)>0)));
        Nuclear_FR(iter,1,:)=nFR;
        all_VOLTAGES(iter,1,:)=v;



        PC_trc=desync_mat;
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
            if v(i+1)>=vt
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
        OVERALL_PC_FR(iter,2)=mean(sum(desync_mat,2))/(length(ax)/1000);

        PC_FR(iter,2,:)=nanmean(desync_mat,1);
        PC_ISI(iter,2)=nanmean(diff(find(mean(desync_mat,1)>0)));
        Nuclear_FR(iter,2,:)=nFR;
        all_VOLTAGES(iter,2,:)=v;


        % ccg_observed=[ccg_observed; get_ccg(real_sync)];
        avg_ccg_observed=[avg_ccg_observed; mean(get_ccg(real_sync)) mean(get_ccg(desync_mat))];
        % all_ccg(iter,1)=get_ccg(no_additional_sync);
        % all_ccg(iter,2)=get_ccg(real_sync);

        clear PC_trc nFR G_all G_this v XX
    end
    all_CbNFR_high(ss,:,:)=OVERALL_NUC_FR;
    all_CCG_high(ss,:,:)=avg_ccg_observed;
    all_ISI_high(ss,:,:)=PC_ISI;
    clear OVERALL_NUC_FR avg_ccg_observed PC_ISI
end
%%
subplot(2,2,1)
errorbar([0 maxshuffle],[0; squeeze(mean(diff(all_CbNFR_high,[],3),2))],[0; squeeze(nanstd(diff(all_CbNFR_high,[],3),[],2))./sqrt(size(all_CbNFR_low,2))],'m','CapSize',0)
hold on
errorbar([0 maxshuffle],[0; squeeze(mean(diff(all_CbNFR_low,[],3),2))],[0; squeeze(nanstd(diff(all_CbNFR_low,[],3),[],2))./sqrt(size(all_CbNFR_low,2))],'g','CapSize',0)
hold off
set(gca,'TickDir','out'); box off
xlabel('jitter (ms)')
ylabel('\DeltaCbN FR')

subplot(2,2,2)
errorbar([0 maxshuffle],[0; squeeze(mean(diff(all_CbNFR_high,[],3)-diff(all_CbNFR_low,[],3),2))],[0; squeeze(nanstd(diff(all_CbNFR_high,[],3)-diff(all_CbNFR_low,[],3),[],2))./sqrt(size(all_CbNFR_low,2))],'k','CapSize',0)
set(gca,'TickDir','out'); box off
xlabel('jitter (ms)')
ylabel('High-Low')

subplot(2,2,3)
errorbar(maxshuffle,mean(all_CbNFR_high(:,:,1),2),nanstd(all_CbNFR_high(:,:,1),[],2)./sqrt(iterations),'m','CapSize',0)
hold on
errorbar(maxshuffle,mean(all_CbNFR_low(:,:,1),2),nanstd(all_CbNFR_low(:,:,1),[],2)./sqrt(iterations),'g','CapSize',0)
hold off
set(gca,'TickDir','out'); box off
xlabel('jitter (ms)')
ylabel('Baseline CbN FR')

subplot(2,2,4)
high2low=all_CbNFR_low-all_CbNFR_high;
errorbar(maxshuffle,squeeze(mean(diff(high2low,[],3),2))./abs(mean(high2low(:,:,1),2)),(squeeze(std(diff(high2low,[],3),[],2))./abs(mean(high2low(:,:,1),2)))./sqrt(iterations),'k','CapSize',0)
yyaxis right
plot(maxshuffle,cumsum([ diff([0; squeeze(mean(diff(high2low,[],3),2))])]),'--b')
ylabel('cumsum(diff(fraction))')

yyaxis left
set(gca,'TickDir','out'); box off
xlabel('jitter (ms)')
ylabel('Fraction explained')