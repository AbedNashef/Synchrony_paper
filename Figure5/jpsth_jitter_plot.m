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
load C:\Users\PersonLab\Desktop\'synchrony paper'\Dylans_data\HH_model\FR4allClusters.mat;
load C:\Users\PersonLab\Desktop\'synchrony paper'\Dylans_data\HH_model\inter_unit_xcorr_dist.mat;
load C:\Users\PersonLab\Desktop\'synchrony paper'\Dylans_data\HH_model\ccg4model.mat;
load E:\Dylans_data\jPSTHs\Cluster_jPETHs\jPSTHs2model.mat;
% load data\FR_gain.mat;
% all_FR=all_fr_high(find(sum(gain,2)==3),:);
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
max_jitter=[200];%:size(all_perc,2);
t=(ax(end)-ax(1)).*1000;

ccg_observed=[];
avg_ccg_observed=[];
dNuc=[];
%% high Sync pairs
% index= find(cov_bfr_move(:,1)<=quantile(cov_bfr_move(:,1),1));
index= find(cov_bfr_move(:,1)>123-10*3 & cov_bfr_move(:,1)<=123+10*3);
% index= find(cov_bfr_move(:,1)>60.7-10.7*3 & cov_bfr_move(:,1)<=60.7+10.7*3);
% index=1:size(all_FR,1);
iter=0;
ccg_observed=[];
avg_ccg_observed=[];
dNuc=[];
for jit=1:length(max_jitter)
    iter=0;

    while iter<iterations
        %
        this_dt=[];
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
            % remove 1-ms sync
            for i=1:size(MAT,1)
                for ii=1:size(MAT,3)
                    X1=squeeze(MAT(i,1,ii,:));
                    X2=squeeze(MAT(i,2,ii,:));
                    new_X1=zeros(1,length(X1));
                    new_X2=zeros(1,length(X2));

                    ix= find(X1==1);
                    tmp1=randn(1,length(ix));
                    dtsign=(tmp1./abs(tmp1)); % shuffle of 1ms!
                    dt=randi(max_jitter(jit),1,length(ix)).*dtsign;
                    rix= ix(:)+dt(:);
                    ix= find(rix<1 | rix>size(MAT,4));
                    rix(ix)=randi(length(X1),1,length(ix));
                    rix=unique(rix);
                    new_X1(rix)=1;
                    spk2add=length(find(X1==1))-length(rix);
                    ix= find(new_X1==0);
                    [~,tmp]= sort(rand(1,length(ix)));
                    new_X1(ix(tmp(1:spk2add)))=1;


                    ix= find(X2==1);
                    tmp1=randn(1,length(ix));
                    dtsign=(tmp1./abs(tmp1)); % shuffle of 1ms!
                    dt=randi(max_jitter(jit),1,length(ix)).*dtsign;
                    rix= ix(:)+dt(:);
                    ix= find(rix<1 | rix>size(MAT,4));
                    rix(ix)=randi(length(X2),1,length(ix));
                    rix=unique(rix);
                    new_X2(rix)=1;
                    spk2add=length(find(X2==1))-length(rix);
                    ix= find(new_X2==0);
                    [~,tmp]= sort(rand(1,length(ix)));
                    new_X2(ix(tmp(1:spk2add)))=1;


                    jittered_1ms(i,1,ii,:)=new_X1;
                    jittered_1ms(i,2,ii,:)=new_X2;
                    this_dt=[this_dt; dt(:)];
                end
            end
            %
            iter=iter+1;
        catch me
            continue
        end
        if ~isempty(dt)
            RANGE_dt(jit,iter)=range(dt);
            avg_dt(jit,iter)=mean(dt);
        end
        for tr=1:size(added_sync,3)
            PC_trc=[];
            for ii=1:size(added_sync,1)
                PC_trc=[PC_trc; squeeze(added_sync(ii,1,tr,:))'];
                PC_trc=[PC_trc; squeeze(added_sync(ii,2,tr,:))'];
            end
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

            OVERALL_NUC_FR(iter,tr,1)=sum(nFR)/(length(ax)/1000);
            OVERALL_PC_FR(iter,tr,1)=mean(sum(PC_trc,2))/(length(ax)/1000);


            PC_FR(iter,tr,1,:)=nanmean(PC_trc,1);

            PC_ISI(iter,tr,1)=nanmean(diff(find(mean(PC_trc,1)>0)));
            Nuclear_FR(iter,tr,1,:)=nFR;
            clear PC_trc nFR G_all G_this v XX

            % jittered
            PC_trc=[];
            for ii=1:size(jittered_1ms,1)
                PC_trc=[PC_trc; squeeze(jittered_1ms(ii,1,tr,:))'];
                PC_trc=[PC_trc; squeeze(jittered_1ms(ii,2,tr,:))'];
            end
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

            OVERALL_NUC_FR(iter,tr,2)=sum(nFR)/(length(ax)/1000);
            OVERALL_PC_FR(iter,tr,2)=mean(sum(PC_trc,2))/(length(ax)/1000);


            PC_FR(iter,tr,2,:)=nanmean(PC_trc,1);

            PC_ISI(iter,tr,2)=nanmean(diff(find(mean(PC_trc,1)>0)));
            Nuclear_FR(iter,i,2,:)=nFR;
            clear PC_trc nFR G_all G_this v XX

        end
        for i=1:size(added_sync,1)
            cell1=squeeze(added_sync(i,1,:,:));
            cell2=squeeze(added_sync(i,2,:,:));
            [raw, shift_predict, pred, surprise, std, ~, ~] = my_JPSTH(cell1', cell2', nbin,1);
            ijpsth_syc(i,:,:)=raw-pred;
            for ii=1:size(cell1,1)
                ccg(ii,:)=get_ccg_whole_with_time([cell1(ii,:); cell2(ii,:)],100);
            end
            this_ccg(i,:)=nanmean(ccg,1);

            cell1=squeeze(jittered_1ms(i,1,:,:));
            cell2=squeeze(jittered_1ms(i,2,:,:));
            [raw, shift_predict, pred, surprise, std, ~, ~] = my_JPSTH(cell1', cell2', nbin,1);
            ijpsth_jit(i,:,:)=raw-pred;
        end
        all_jpsth_real(iter,:,:)=squeeze(nanmean(ijpsth_syc,1));
        all_jpsth_jit(iter,:,:)=squeeze(nanmean(ijpsth_jit,1));
        diff1ms_sync(iter,:)=diag(squeeze(nanmean(ijpsth_jit-ijpsth_syc,1)));
        PC_overall_FR(iter,1,:)=squeeze(nanmean(nanmean(nanmean(added_sync,1),2),3));
        PC_overall_FR(iter,2,:)=squeeze(nanmean(nanmean(nanmean(jittered_1ms,1),2),3));
        CCG_all(iter,:)=nanmean(this_ccg,1);
        diag_real(iter,:)=diag(squeeze(nanmean(ijpsth_syc,1)));
    end
    %
    subplot(2,length(max_jitter)+1,1)
    pcolor(ax,ax,squeeze(nanmean(all_jpsth_real,1))./ntr);
    shading flat
    title('high; real')
    subplot(2,length(max_jitter)+1,jit+1)
    pcolor(ax,ax,squeeze(nanmean(all_jpsth_jit,1))./ntr);
    shading flat
    title(['high; jitter=' num2str(max_jitter(jit))])
    clear OVERALL_NUC_FR PC_ISI PC_FR  OVERALL_PC_FR all_jpsth_jit all_jpsth_real diff1ms_sync CCG_all
    % for i=1:size(all_jpsth_jit,1)
    % all_d(i,:)=diag(squeeze(all_jpsth_jit(i,:,:)-all_jpsth_real(i,:,:)));
    % end
    %
end

%% high Sync pairs
% index= find(cov_bfr_move(:,1)<=quantile(cov_bfr_move(:,1),1));
% index= find(cov_bfr_move(:,1)>123-10*3 & cov_bfr_move(:,1)<=123+10*3);
index= find(cov_bfr_move(:,1)>60.7-10.7*3 & cov_bfr_move(:,1)<=60.7+10.7*3);
% index=1:size(all_FR,1);
iter=0;
ccg_observed=[];
avg_ccg_observed=[];
dNuc=[];
for jit=1:length(max_jitter)
    iter=0;

    while iter<iterations
        %
        this_dt=[];
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
            % remove 1-ms sync
            for i=1:size(MAT,1)
                for ii=1:size(MAT,3)
                    X1=squeeze(MAT(i,1,ii,:));
                    X2=squeeze(MAT(i,2,ii,:));
                    new_X1=zeros(1,length(X1));
                    new_X2=zeros(1,length(X2));

                    ix= find(X1==1);
                    tmp1=randn(1,length(ix));
                    dtsign=(tmp1./abs(tmp1)); % shuffle of 1ms!
                    dt=randi(max_jitter(jit),1,length(ix)).*dtsign;
                    rix= ix(:)+dt(:);
                    ix= find(rix<1 | rix>size(MAT,4));
                    rix(ix)=randi(length(X1),1,length(ix));
                    rix=unique(rix);
                    new_X1(rix)=1;
                    spk2add=length(find(X1==1))-length(rix);
                    ix= find(new_X1==0);
                    [~,tmp]= sort(rand(1,length(ix)));
                    new_X1(ix(tmp(1:spk2add)))=1;


                    ix= find(X2==1);
                    tmp1=randn(1,length(ix));
                    dtsign=(tmp1./abs(tmp1)); % shuffle of 1ms!
                    dt=randi(max_jitter(jit),1,length(ix)).*dtsign;
                    rix= ix(:)+dt(:);
                    ix= find(rix<1 | rix>size(MAT,4));
                    rix(ix)=randi(length(X2),1,length(ix));
                    rix=unique(rix);
                    new_X2(rix)=1;
                    spk2add=length(find(X2==1))-length(rix);
                    ix= find(new_X2==0);
                    [~,tmp]= sort(rand(1,length(ix)));
                    new_X2(ix(tmp(1:spk2add)))=1;


                    jittered_1ms(i,1,ii,:)=new_X1;
                    jittered_1ms(i,2,ii,:)=new_X2;
                    this_dt=[this_dt; dt(:)];
                end
            end
            %
            iter=iter+1;
        catch me
            continue
        end
        if ~isempty(dt)
            RANGE_dt(jit,iter)=range(dt);
            avg_dt(jit,iter)=mean(dt);
        end
        for tr=1:size(added_sync,3)
            PC_trc=[];
            for ii=1:size(added_sync,1)
                PC_trc=[PC_trc; squeeze(added_sync(ii,1,tr,:))'];
                PC_trc=[PC_trc; squeeze(added_sync(ii,2,tr,:))'];
            end
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

            OVERALL_NUC_FR(iter,tr,1)=sum(nFR)/(length(ax)/1000);
            OVERALL_PC_FR(iter,tr,1)=mean(sum(PC_trc,2))/(length(ax)/1000);


            PC_FR(iter,tr,1,:)=nanmean(PC_trc,1);

            PC_ISI(iter,tr,1)=nanmean(diff(find(mean(PC_trc,1)>0)));
            Nuclear_FR(iter,tr,1,:)=nFR;
            clear PC_trc nFR G_all G_this v XX

            % jittered
            PC_trc=[];
            for ii=1:size(jittered_1ms,1)
                PC_trc=[PC_trc; squeeze(jittered_1ms(ii,1,tr,:))'];
                PC_trc=[PC_trc; squeeze(jittered_1ms(ii,2,tr,:))'];
            end
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

            OVERALL_NUC_FR(iter,tr,2)=sum(nFR)/(length(ax)/1000);
            OVERALL_PC_FR(iter,tr,2)=mean(sum(PC_trc,2))/(length(ax)/1000);


            PC_FR(iter,tr,2,:)=nanmean(PC_trc,1);

            PC_ISI(iter,tr,2)=nanmean(diff(find(mean(PC_trc,1)>0)));
            Nuclear_FR(iter,i,2,:)=nFR;
            clear PC_trc nFR G_all G_this v XX

        end
        for i=1:size(added_sync,1)
            cell1=squeeze(added_sync(i,1,:,:));
            cell2=squeeze(added_sync(i,2,:,:));
            [raw, shift_predict, pred, surprise, std, ~, ~] = my_JPSTH(cell1', cell2', nbin,1);
            ijpsth_syc(i,:,:)=raw-pred;
            for ii=1:size(cell1,1)
                ccg(ii,:)=get_ccg_whole_with_time([cell1(ii,:); cell2(ii,:)],100);
            end
            this_ccg(i,:)=nanmean(ccg,1);

            cell1=squeeze(jittered_1ms(i,1,:,:));
            cell2=squeeze(jittered_1ms(i,2,:,:));
            [raw, shift_predict, pred, surprise, std, ~, ~] = my_JPSTH(cell1', cell2', nbin,1);
            ijpsth_jit(i,:,:)=raw-pred;
        end
        all_jpsth_real(iter,:,:)=squeeze(nanmean(ijpsth_syc,1));
        all_jpsth_jit(iter,:,:)=squeeze(nanmean(ijpsth_jit,1));
        diff1ms_sync(iter,:)=diag(squeeze(nanmean(ijpsth_jit-ijpsth_syc,1)));
        PC_overall_FR(iter,1,:)=squeeze(nanmean(nanmean(nanmean(added_sync,1),2),3));
        PC_overall_FR(iter,2,:)=squeeze(nanmean(nanmean(nanmean(jittered_1ms,1),2),3));
        CCG_all(iter,:)=nanmean(this_ccg,1);
        diag_real(iter,:)=diag(squeeze(nanmean(ijpsth_syc,1)));
    end
    %
    subplot(2,length(max_jitter)+1,length(max_jitter)+2)
    pcolor(ax,ax,squeeze(nanmean(all_jpsth_real,1))./ntr);
    shading flat
    title('low; real')
    subplot(2,length(max_jitter)+1,jit+length(max_jitter)+2)
    pcolor(ax,ax,squeeze(nanmean(all_jpsth_jit,1))./ntr);
    shading flat
    title(['low; jitter=' num2str(max_jitter(jit))])
    clear OVERALL_NUC_FR PC_ISI PC_FR  OVERALL_PC_FR all_jpsth_jit all_jpsth_real diff1ms_sync CCG_all
    % for i=1:size(all_jpsth_jit,1)
    % all_d(i,:)=diag(squeeze(all_jpsth_jit(i,:,:)-all_jpsth_real(i,:,:)));
    % end
    %
end