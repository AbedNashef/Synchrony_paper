% model- panel k
% from Wu and Raman.
% MF-CbN convergence: 20-600 MF to CbN
% AMPA: Rise time: 0.28 ms, tau=1.06 ms; reversal potential: 11.5 mV;
% conductance: 1nS?
clear
C=89; vr=-60; vt=-40;
a=.03; b=-2; c=-50; d=2;
% a=3;
gmax_Inh=5; %nS
k=gmax_Inh/C; 
vpeak= 35; % spike cutoff
Er_Inh=-75;
T=100;%s preceded by a brief increase in activity, presumably due to excitation of the DCN neuron by col0;
dt=1; % time span and strp (ms)
n=round(T/dt); % number of simulation steps
tau_Inh=5;%2.5; % based on Person and Raman, 2012
tau_rise_Inh=.05;

nMF_range=[20 600];
Er_MF=11.5;
tau_MF=1.06;
tau_rise_MF=.28;
gmax_MF=1; %nS
MF_baselineFR=15; %sp/s

rp=2;
vplot=zeros(2,2);
win=50;
c_all=0;
ex_raster=[];
% mf_ex_raster=cell(0);
addpath ../../Dylans_data/helper_functions/
%%
load ../../Data/Figure4/V2model_alltrials.mat;
load ../../Data/Figure4/FR4allClusters.mat;
load ../../Data/Figure4/inter_unit_xcorr_dist.mat;
load ../../Data/Figure4/ccg4model.mat;
all_V=all_xv_allTrials;
nPC=40;
this_xcorr=all_ccg(find(sum(all_clusters,2)==3),find(xcorr_ax==0));
% this_xcorr=this_xcorr(find(this_xcorr>median(this_xcorr)));
% this_xcorr=xcorr_dist{1};
% this_xcorr=this_xcorr(find(this_xcorr>0 & this_xcorr<20));
this_xcorr=this_xcorr(this_xcorr>=quantile(this_xcorr,.8));
this_xcorr=this_xcorr.*5;
% this_xcorr(find(this_xcorr<0))=0;

xcorr_1_1=4; %4 sp/s cross-corr
n2sync= round((length(ax)/1000)*xcorr_1_1);
iterations=4000;
jitter=1;%:size(all_perc,2);
t=(ax(end)-ax(1))*1000;

ccg_observed=[];
avg_ccg_observed=[];
dNuc=[];

ix=(find(bhv_ax>-.3 & bhv_ax<.3));
strectched_V=all_V(:,ix);
new_bhv_ax=linspace(bhv_ax(1),bhv_ax(end),length(ix));
%%

iter=0;
while iter<iterations
    nMF=randi(nMF_range,1);
    baseline_fr=poissrnd(MF_baselineFR,nMF,1);
    clear MF_TRC
    % for imf=1:nMF
    %     ex_input=interp1(bhv_ax,all_V(randi(size(all_V,1)),:),ax);
    %     MF_TRC(imf,:) = simulateModulatedSpikes(ex_input./max(ex_input), 1000, poissrnd(MF_baselineFR),rand(1)-.2, 1, 1e-3);
    % end
    ex_input=interp1(new_bhv_ax,strectched_V(randi(size(all_V,1)),:),ax);
    MF_TRC = simulateModulatedSpikes(ex_input./max(ex_input), 1000, poissrnd(MF_baselineFR),rand(1), nMF, 1e-3);

    c1_ix= find([all_FR.cluster]<=2);
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
    I=zeros(1,t+1);
    for i=1:nMF
        ix= find(MF_TRC(i,:));
        G_MF=zeros(length(ix),t+1);
        for j=1:length(ix)
            for ii=ix(j):ix(j)+n
                tpeak=ix(j)+(tau_MF*tau_rise_MF/(tau_MF-tau_rise_MF))*log(tau_MF/tau_rise_MF);
                f=1/(-1*exp(-1*(tpeak-ix(j))/tau_rise_MF)+exp(-1*(tpeak-ix(j))/tau_MF));
                G_MF(j,ii)=gmax_MF*f*(exp(-1*(ii-ix(j))/tau_MF)-exp(-1*((ii-ix(j))/tau_rise_MF)));
            end
        end
        G_MF=G_MF(:,1:length(ax));
        G_MF_all(i,:)= sum(G_MF,1);
        %             G_this(i,:)= (G_this(i,:)./max(G_this(i,:))).*gmax;
    end
    for i=1:nPC
        ix= find(PC_trc(i,:));
        G_all=zeros(length(ix),t+1);
        for j=1:length(ix)
            for ii=ix(j):ix(j)+n
                tpeak=ix(j)+(tau_Inh*tau_rise_Inh/(tau_Inh-tau_rise_Inh))*log(tau_Inh/tau_rise_Inh);
                f=1/(-1*exp(-1*(tpeak-ix(j))/tau_rise_Inh)+exp(-1*(tpeak-ix(j))/tau_Inh));
                G_all(j,ii)=gmax_Inh*f*(exp(-1*(ii-ix(j))/tau_Inh)-exp(-1*((ii-ix(j))/tau_rise_Inh)));
            end
        end
        G_all=G_all(:,1:length(ax));
        G_this(i,:)= sum(G_all,1);
        %             G_this(i,:)= (G_this(i,:)./max(G_this(i,:))).*gmax;
    end
    %
    % Inhibitory current
    for i=1:size(PC_trc,2)
        for ii=1:size(G_this,1)
            I(i)= I(i) + G_this(ii,i)*(Er_Inh-v(i));
            I_inh(i)=G_this(ii,i)*(Er_Inh-v(i));
        end
        for ii=1:size(G_MF_all,1)
            I(i)= I(i)+ G_MF_all(ii,i)*(Er_MF-v(i));
            I_MF(i)=G_MF_all(ii,i)*(Er_MF-v(i));
        end
        v(i+1) = v(i)+(k*(v(i)-vr)*(v(i)-vt)-u(i)+I(i))/C+(k*(v(i)-vr)*(v(i)-vt)-u(i)+I(i))/C;
        % v(i+1) = v(i)+(k*(v(i)-vr)*(v(i)-vt)-u(i)+I(i))/C+(k*(v(i)-vr)*(v(i)-vt)-u(i)+I(i))/C;
        u(i+1) = u(i)+tau_Inh*a*(b*(v(i)-vr)-u(i))+tau_MF*a*(b*(v(i)-vr)-u(i));
        if v(i+1)>=vt%vpeak
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
    OVERALL_NUC_FR_move(iter,1)=sum(nFR(find(ax>-.2 & ax<.0)))./.2;
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
                tpeak=ix(j)+(tau_Inh*tau_rise_Inh/(tau_Inh-tau_rise_Inh))*log(tau_Inh/tau_rise_Inh);
                f=1/(-1*exp(-1*(tpeak-ix(j))/tau_rise_Inh)+exp(-1*(tpeak-ix(j))/tau_Inh));
                G_all(j,ii)=gmax_Inh*f*(exp(-1*(ii-ix(j))/tau_Inh)-exp(-1*((ii-ix(j))/tau_rise_Inh)));
            end
        end
        G_all=G_all(:,1:length(ax));
        G_this(i,:)= sum(G_all,1);
        %             G_this(i,:)= (G_this(i,:)./max(G_this(i,:))).*gmax;
    end
   for i=1:size(PC_trc,2)
        for ii=1:size(G_this,1)
            I(i)= I(i) + G_this(ii,i)*(Er_Inh-v(i));
            I_inh(i)=G_this(ii,i)*(Er_Inh-v(i));
        end
        for ii=1:size(G_MF_all,1)
            I(i)= I(i)+ G_MF_all(ii,i)*(Er_MF-v(i));
            I_MF(i)=G_MF_all(ii,i)*(Er_MF-v(i));
        end
        % v(i+1) = v(i)+tau_Inh*(k*(v(i)-vr)*(v(i)-vt)-u(i)+I(i))/C+tau_MF*(k*(v(i)-vr)*(v(i)-vt)-u(i)+I(i))/C;
        v(i+1) = v(i)+(k*(v(i)-vr)*(v(i)-vt)-u(i)+I(i))/C+(k*(v(i)-vr)*(v(i)-vt)-u(i)+I(i))/C;
        u(i+1) = u(i)+tau_Inh*a*(b*(v(i)-vr)-u(i))+tau_MF*a*(b*(v(i)-vr)-u(i));
        if v(i+1)>=vt%vpeak
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
    OVERALL_NUC_FR_move(iter,2)=sum(nFR(find(ax>-.2 & ax<.0)))./.2;
    OVERALL_PC_FR(iter,2)=mean(sum(added_sync,2))/(length(ax)/1000);

    if OVERALL_NUC_FR(iter,2)-OVERALL_NUC_FR(iter,1)>20
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

    
    for nn=1:size(added_sync,1)
        null_fr(nn,:)=get_null_fr(added_sync(nn,:),1);
    end
    cofiring(iter,2,:)=nanmean(added_sync-null_fr,1);
    for nn=1:size(added_sync,1)
        null_fr(nn,:)=get_null_fr(no_additional_sync(nn,:),1);
    end
    cofiring(iter,1,:)=nanmean(no_additional_sync-null_fr,1);
    ccg_observed=[ccg_observed; get_ccg(no_additional_sync) get_ccg(added_sync)];
    avg_ccg_observed=[avg_ccg_observed; mean(get_ccg(no_additional_sync)) mean(get_ccg(added_sync))];
    % all_ccg(iter,1)=get_ccg(no_additional_sync);
    % all_ccg(iter,2)=get_ccg(added_sync);
    all_ccg(iter)=nanmean(all_n2sync(:))/(length(ax)/1000);
    mf_ex_raster(iter,:)=nanmean(MF_TRC,1);
    avg_ex(iter)=mean(sum(G_MF_all,1));
    avg_ex_mean(iter)=nanmean(nanmean(G_MF_all,1));
    num_of_MF(iter)=nMF;
    clear PC_trc nFR G_all G_this v XX G_MF_all
end
%%
subplot(2,2,1)
scatter(avg_ex,OVERALL_NUC_FR(:,1),'.k')
hold on
scatter(avg_ex,OVERALL_NUC_FR(:,2),'.r')
lsline
hold off

bin=1;
dax=[round((min(avg_ex)/bin))*bin:bin:max(avg_ex)+bin];
clear std_fr avg_fr avg_dfr std_dfr  avg_decel std_decel avg_decelGain std_decelGain
for i=1:length(dax)
    ix= find(avg_ex>=dax(i) & avg_ex<dax(i)+bin);
    % if length(ix)<10
    %     OVERALL_NUC_FR(ix,:)=nan;
    % end
    avg_fr(i,:)=nanmean(OVERALL_NUC_FR(ix,:),1);
    std_fr(i,:)=nanstd(OVERALL_NUC_FR(ix,:),[],1)./sqrt(length(ix));

    avg_dfr(i,:)=nanmean(diff(OVERALL_NUC_FR(ix,:),[],2),1);
    std_dfr(i,:)=nanstd(diff(OVERALL_NUC_FR(ix,:),[],2),[],1)./sqrt(length(ix));

    avg_decel(i,:)=nanmean(diff(OVERALL_NUC_FR(ix,:),[],2).*1.2371+15.113,1); % based on Becker  & Person + my kinematic analysis
    std_decel(i,:)=nanstd(diff(OVERALL_NUC_FR(ix,:),[],2).*1.2371+15.113,[],1)./sqrt(length(ix));

    % avg_decelGain(i,:)=nanmean(diff(OVERALL_NUC_FR(ix,:),[],2).*0.41183+70.175,1); % based on Becker  & Person + my kinematic analysis // sorted by decel
    % std_decelGain(i,:)=nanstd(diff(OVERALL_NUC_FR(ix,:),[],2).*0.41183+70.175,[],1)./sqrt(length(ix));

    avg_decelGain(i,:)=nanmean(diff(OVERALL_NUC_FR(ix,:),[],2).*.63973+54.266,1); % based on Becker  & Person + my kinematic analysis // sorted by endpoint (hypo-hyper)
    std_decelGain(i,:)=nanstd(diff(OVERALL_NUC_FR(ix,:),[],2).*.63973+54.266,[],1)./sqrt(length(ix));

    avg_dfr_move(i,:)=nanmean(diff(OVERALL_NUC_FR_move(ix,:),[],2),1);
    std_dfr_move(i,:)=nanstd(diff(OVERALL_NUC_FR_move(ix,:),[],2),[],1)./sqrt(length(ix));
end
subplot(2,4,3)
errorbar(dax,avg_fr,std_fr)
set(gca,'TickDir','out'); box off

subplot(2,4,4)
errorbar(dax,avg_dfr,std_dfr,'LineWidth',2,'Color','k','CapSize',0)
% hold on
% errorbar(dax,avg_dfr_move,std_dfr_move,'LineWidth',2,'Color','g','CapSize',0)
% hold off
set(gca,'TickDir','out'); box off

subplot(2,3,4)
scatter(diff(avg_ccg_observed,[],2),diff(OVERALL_NUC_FR,[],2),'.b')
hold on
a=axis;
plot(a(1:2),[0 0],'--k')
lsline
xlabel('CCG')
ylabel('\DeltaFR')
h=fitlm(diff(avg_ccg_observed,[],2),diff(OVERALL_NUC_FR,[],2));
title(h.Coefficients.Estimate(2))
set(gca,'TickDir','out'); box off
%%
% Sample data

% Fit a linear regression model
mdl = fitlm(avg_ex, OVERALL_NUC_FR(:,1));

% Get predictions and confidence intervals
[x_sorted, idx] = sort(avg_ex(:)); % Sort x for smooth plotting
[y_pred, y_ci] = predict(mdl, x_sorted);

% Plot the data
subplot(2,3,5);
hold on;
% scatter(x, y, 'b', 'filled'); % Scatter plot of data
plot(x_sorted, y_pred, 'k-', 'LineWidth', 2); % Regression line

% Plot confidence intervals
fill([x_sorted; flipud(x_sorted)], [y_ci(:,1); flipud(y_ci(:,2))], ...
    'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none'); % Confidence bounds

mdl = fitlm(avg_ex, OVERALL_NUC_FR(:,2));

% Get predictions and confidence intervals
[x_sorted, idx] = sort(avg_ex(:)); % Sort x for smooth plotting
[y_pred, y_ci] = predict(mdl, x_sorted);

% scatter(x, y, 'b', 'filled'); % Scatter plot of data
plot(x_sorted, y_pred, 'b-', 'LineWidth', 2); % Regression line

% Plot confidence intervals
fill([x_sorted; flipud(x_sorted)], [y_ci(:,1); flipud(y_ci(:,2))], ...
    'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none'); % Confidence bounds


% Add labels and legend
set(gca,'TickDir','out'); box off
xlabel('Avergae excitation');
ylabel('CbN firing rate');
hold off;

%%
% subplot(2,4,7)
% scatter(randn(1,iter),diff(OVERALL_NUC_FR_move,[],2),10,'ob','fill')
% set(gca,'TickDir','out'); box off
subplot(2,3,6)
scatter(randn(1,iter),diff(OVERALL_NUC_FR,[],2),10,'ob','fill')
set(gca,'TickDir','out'); box off


%%
X=diff(squeeze(nanmean(cofiring,1)),1);
x1=squeeze(nanmean(cofiring,1));
s=diff(squeeze(nanstd(cofiring,1)),[],1)./sqrt(iter);
gwin=gausswin(151);
gwin=gwin./sum(gwin);
plot(ax,conv2(X,gwin','same'))
%%
a=1.1070; b=-148.6085; 
real_gain_range=[8.8695 93.6987];
median_gain_change=38.2600;
DEC_GAIN=OVERALL_NUC_FR_move*a+b;
% ix = find(DEC_GAIN(:,1)<0 & DEC_GAIN(:,2)>0)
ix=find(OVERALL_NUC_FR_move(:,1)>0 & OVERALL_NUC_FR_move(:,2)>0);% & OVERALL_NUC_FR_move(:,1)<400 & OVERALL_NUC_FR_move(:,2)<400);
figure
subplot 221
histogram(OVERALL_NUC_FR_move(ix,1).*a+b,'BinWidth',20,'FaceColor','k')
hold on
histogram(OVERALL_NUC_FR_move(ix,2).*a+b,'BinWidth',20,'FaceColor','b')
hold off
h=get(gca,'Children');
x1=h(1).BinEdges;
y1=h(1).BinCounts;
x2=h(2).BinEdges;
y2=h(2).BinCounts;

subplot 222
plot(x1(1:end-1),y1./sum(y1),'b')
hold on
plot(x2(1:end-1),y2./sum(y2),'k')
hold off

subplot 223
histogram(OVERALL_NUC_FR_move(ix,2).*a-OVERALL_NUC_FR_move(ix,1).*a,'BinWidth',10);
hold on
aa=axis;
plot(real_gain_range,[aa(4)-1 aa(4)-1],'k','LineWidth',2)
scatter(median_gain_change,aa(4),50,'or')

clear X S
dax=[0:5:floor(max(diff(OVERALL_NUC_FR_move,[],2))+5)];
for ii=1:length(dax)-1
    ix= find(diff(OVERALL_NUC_FR_move,[],2)>=dax(ii) & diff(OVERALL_NUC_FR_move,[],2)<=dax(ii+1));
    X(ii)= mean(OVERALL_NUC_FR_move(ix,2).*a-OVERALL_NUC_FR_move(ix,1).*a);
    S(ii)=std(OVERALL_NUC_FR_move(ix,2).*a-OVERALL_NUC_FR_move(ix,1).*a);%./sqrt(length(ix));
end
subplot 224
errorbar(dax(1:end-1),X,S)

clear X S C
diffCCG=diff(avg_ccg_observed,[],2);
bin=1;
dax=[round((min(avg_ex)/bin))-bin:bin:max(avg_ex)+bin];
dax=linspace(min(avg_ex),max(avg_ex),11);
for ii=1:length(dax)-1
    ix= find(avg_ex>=dax(ii) & avg_ex<dax(ii+1));
    if length(ix)>15
        X(ii)= mean(OVERALL_NUC_FR_move(ix,2).*a-OVERALL_NUC_FR_move(ix,1).*a);
        C(ii)=mean(num_of_MF(ix));
        S(ii)=std(OVERALL_NUC_FR_move(ix,2).*a-OVERALL_NUC_FR_move(ix,1).*a)./sqrt(length(ix));
    else
        X(ii)=nan;
        C(ii)=nan;
        S(ii)=nan;
    end
end
subplot 224
errorbar(dax(1:end-1),X,S)
hold on
scatter(dax(1:end-1),X,50,C,'fill')
hold off


clear X S C
diffCCG=diff(avg_ccg_observed,[],2);
bin=1;
dax=linspace(nMF_range(1),nMF_range(2),21);
for ii=1:length(dax)-1
    ix= find(num_of_MF>=dax(ii) & num_of_MF<dax(ii+1));
    if length(ix)>15
        X(ii)= mean(OVERALL_NUC_FR_move(ix,2).*a-OVERALL_NUC_FR_move(ix,1).*a);
        C(ii)=mean(avg_ex(ix));
        S(ii)=std(OVERALL_NUC_FR_move(ix,2).*a-OVERALL_NUC_FR_move(ix,1).*a)./sqrt(length(ix));
    else
        X(ii)=nan;
        C(ii)=nan;
        S(ii)=nan;
    end
end
subplot 224
errorbar(dax(1:end-1),X,S)
hold on
scatter(dax(1:end-1),X,50,C,'fill')
hold off

clear X S C fr_all fr_move sfr_all sfr_move
bin=.001;
dax=[round((min(diffCCG)))-bin:bin:max(diffCCG)+bin];
dax=linspace(min(diffCCG),max(diffCCG),11);
for ii=1:length(dax)-1
    ix= find(diffCCG>=dax(ii) & diffCCG<dax(ii+1));
    if length(ix)>15
        X(ii)= mean(OVERALL_NUC_FR_move(ix,2).*a-OVERALL_NUC_FR_move(ix,1).*a);
        S(ii)=std(OVERALL_NUC_FR_move(ix,2).*a-OVERALL_NUC_FR_move(ix,1).*a)./sqrt(length(ix));
        C(ii)=mean(avg_ex(ix));
        fr_move(ii)= mean(OVERALL_NUC_FR_move(ix,2)-OVERALL_NUC_FR_move(ix,1));
        fr_all(ii)= mean(OVERALL_NUC_FR(ix,2)-OVERALL_NUC_FR(ix,1));

        sfr_move(ii)= std(OVERALL_NUC_FR_move(ix,2)-OVERALL_NUC_FR_move(ix,1))./sqrt(length(ix));
        sfr_all(ii)= std(OVERALL_NUC_FR(ix,2)-OVERALL_NUC_FR(ix,1))./sqrt(length(ix));

    else
        X(ii)=nan;
        S(ii)=nan;
        C(ii)=nan;
        fr_move(ii)=nan;
        fr_all(ii)=nan;
                sfr_move(ii)=nan;
        sfr_all(ii)=nan;
    end
end
subplot 224
errorbar(dax(1:end-1),X,S)
hold on
scatter(dax(1:end-1),X,50,C,'fill')
hold off
mdl=fitlm(dax(1:end-1),X);


clear X S C fr_all fr_move
bin=.001;
decelGain=OVERALL_NUC_FR_move(:,2).*a-OVERALL_NUC_FR_move(:,1).*a;
dax=linspace(min(decelGain),max(decelGain),21);
for ii=1:length(dax)-1
    ix= find(decelGain>=dax(ii) & decelGain<dax(ii+1));
    if length(ix)>15
        X(ii)= mean(avg_ex(ix));
        S(ii)=std(avg_ex(ix))./sqrt(length(ix));
        C(ii)=mean(num_of_MF(ix));
        fr_move(ii)= mean(OVERALL_NUC_FR_move(ix,2)-OVERALL_NUC_FR_move(ix,1));
        fr_all(ii)= mean(OVERALL_NUC_FR(ix,2)-OVERALL_NUC_FR(ix,1));
    else
        X(ii)=nan;
        S(ii)=nan;
        C(ii)=nan;
        fr_move(ii)=nan;
        fr_all(ii)=nan;
    end
end
subplot 224
errorbar(dax(1:end-1),X,S)
hold on
scatter(dax(1:end-1),X,50,C,'fill')
hold off
mdl_2MFnum=fitlm(dax(1:end-1),X);
mdl_2Ex=fitlm(dax(1:end-1),C);
%%
%%
this_xcorr=this_xcorr./5;
%%
[x_sorted, idx] = sort(this_xcorr(:)); % Sort x for smooth plotting
[y_pred, y_ci] = predict(mdl, x_sorted);

% Plot the data
subplot(2,2,4);
hold on;
% scatter(x, y, 'b', 'filled'); % Scatter plot of data
plot(x_sorted, y_pred, 'k-', 'LineWidth', 2); % Regression line

% Plot confidence intervals
fill([x_sorted; flipud(x_sorted)], [y_ci(:,1); flipud(y_ci(:,2))], ...
    'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none'); % Confidence bounds

%%
figure
m=nanmean(mf_ex_raster,1).*1000;
s=(nanstd(mf_ex_raster,[],1).*1000)./sqrt(size(mf_ex_raster,1));
subplot 221
patch([ax fliplr(ax)],[m-s fliplr(m+s)],[.5 .5 .5],'EdgeColor','none')
hold on
plot(ax,m,'k','LineWidth',1)
hold off

subplot(2,2,2)
scatter(avg_ex_mean.*1000,OVERALL_NUC_FR(:,1),'.k')
lsline
%%
a=1.1070; b=-148.6085; 
real_gain_range=[8.8695 93.6987];
median_gain_change=38.2600;

ix_move=find(ax>-.3 & ax<.2);
for i=1:size(Nuclear_FR,1)
    for ii=1:size(Nuclear_FR,2)
        CbN_move(i,ii)=sum(squeeze(Nuclear_FR(i,ii,ix_move)))/(length(ix_move)/length(ax));
    end
end

clear X S C fr_all fr_move
bin=.001;
decelGain=CbN_move(:,2).*a-CbN_move(:,1).*a;
dax=linspace(min(decelGain),max(decelGain),21);
for ii=1:length(dax)-1
    ix= find(decelGain>=dax(ii) & decelGain<dax(ii+1));
    if length(ix)>15
        X(ii)= mean(avg_ex(ix));
        S(ii)=std(avg_ex(ix))./sqrt(length(ix));
        C(ii)=mean(num_of_MF(ix));
        fr_move(ii)= mean(CbN_move(ix,2)-CbN_move(ix,1));
        fr_all(ii)= mean(OVERALL_NUC_FR(ix,2)-OVERALL_NUC_FR(ix,1));
    else
        X(ii)=nan;
        S(ii)=nan;
        C(ii)=nan;
        fr_move(ii)=nan;
        fr_all(ii)=nan;
    end
end
subplot 221
errorbar(dax(1:end-1),X,S)
hold on
scatter(dax(1:end-1),X,50,C,'fill')
hold off
mdl_2MFnum=fitlm(dax(1:end-1),X);
mdl_2Ex=fitlm(dax(1:end-1),C);


clear X S C fr_all fr_move
dax=linspace(min(avg_ex),max(avg_ex),16);
for ii=1:length(dax)-1
    ix= find(avg_ex>=dax(ii) & avg_ex<dax(ii+1));
    if length(ix)>15
        X(ii)= mean(decelGain(ix));
        S(ii)=std(decelGain(ix))./sqrt(length(ix));
        C(ii)=mean(num_of_MF(ix));
        fr_move(ii)= mean(CbN_move(ix,2)-CbN_move(ix,1));
        fr_all(ii)= mean(OVERALL_NUC_FR(ix,2)-OVERALL_NUC_FR(ix,1));
    else
        X(ii)=nan;
        S(ii)=nan;
        C(ii)=nan;
        fr_move(ii)=nan;
        fr_all(ii)=nan;
    end
end
subplot 221
errorbar(dax(1:end-1),X,S)
hold on
scatter(dax(1:end-1),X,50,C,'fill')
hold off
mdl_2MFnum=fitlm(dax(1:end-1),X);
mdl_2Ex=fitlm(dax(1:end-1),C);



subplot 222
histogram(CbN_move(:,2).*a-CbN_move(:,1).*a,'BinWidth',10);
hold on
aa=axis;
plot(real_gain_range,[aa(4)-1 aa(4)-1],'k','LineWidth',2)
scatter(median_gain_change,aa(4),50,'or')




clear X S C fr_all fr_move sfr_all sfr_move
bin=.001;
dax=[round((min(diffCCG)))-bin:bin:max(diffCCG)+bin];
dax=linspace(min(diffCCG),max(diffCCG),11);
for ii=1:length(dax)-1
    ix= find(diffCCG>=dax(ii) & diffCCG<dax(ii+1));
    if length(ix)>15
        X(ii)= mean(CbN_move(ix,2).*a-CbN_move(ix,1).*a);
        S(ii)=std(CbN_move(ix,2).*a-CbN_move(ix,1).*a)./sqrt(length(ix));
        C(ii)=mean(avg_ex(ix));
        fr_move(ii)= mean(CbN_move(ix,2)-CbN_move(ix,1));
        fr_all(ii)= mean(OVERALL_NUC_FR(ix,2)-OVERALL_NUC_FR(ix,1));

        sfr_move(ii)= std(CbN_move(ix,2)-CbN_move(ix,1))./sqrt(length(ix));
        sfr_all(ii)= std(OVERALL_NUC_FR(ix,2)-OVERALL_NUC_FR(ix,1))./sqrt(length(ix));

    else
        X(ii)=nan;
        S(ii)=nan;
        C(ii)=nan;
        fr_move(ii)=nan;
        fr_all(ii)=nan;
        sfr_move(ii)=nan;
        sfr_all(ii)=nan;
    end
end
subplot 236
errorbar(dax(1:end-1),X,S)
hold on
scatter(dax(1:end-1),X,50,C,'fill')
hold off
mdl=fitlm(dax(1:end-1),X);

subplot(2,3,5)
errorbar(dax(1:end-1),fr_all,sfr_all)
hold on
errorbar(dax(1:end-1),fr_move,sfr_move)
scatter(dax(1:end-1),fr_all,'ok')
scatter(dax(1:end-1),fr_move,'ok')
hold off


clear X S C fr_all fr_move sfr_all sfr_move
bin=.001;
dax=[round((min(diffCCG)))-bin:bin:max(diffCCG)+bin];
dax=linspace(min(diffCCG),max(diffCCG),11);
for ii=1:length(dax)-1
    ix= find(diffCCG>=dax(ii) & diffCCG<dax(ii+1));
    if length(ix)>15
        X(ii)= mean(CbN_move(ix,2).*a-CbN_move(ix,1).*a);
        S(ii)=std(CbN_move(ix,2).*a-CbN_move(ix,1).*a)./sqrt(length(ix));
        C(ii)=mean(avg_ex(ix));
        fr_move(ii)= mean(CbN_move(ix,2)-CbN_move(ix,1));
        fr_all(ii)= mean(OVERALL_NUC_FR(ix,2)-OVERALL_NUC_FR(ix,1));

        sfr_move(ii)= std(CbN_move(ix,2)-CbN_move(ix,1))./sqrt(length(ix));
        sfr_all(ii)= std(OVERALL_NUC_FR(ix,2)-OVERALL_NUC_FR(ix,1))./sqrt(length(ix));

    else
        X(ii)=nan;
        S(ii)=nan;
        C(ii)=nan;
        fr_move(ii)=nan;
        fr_all(ii)=nan;
                sfr_move(ii)=nan;
        sfr_all(ii)=nan;
    end
end
subplot 234
errorbar(dax(1:end-1),X,S)
hold on
scatter(dax(1:end-1),X,50,C,'fill')
hold off
mdl=fitlm(dax(1:end-1),X);
%
[x_sorted, idx] =sort([0:0.005:max(this_xcorr)]');% sort(this_xcorr(:)); % Sort x for smooth plotting
[y_pred, y_ci] = predict(mdl, x_sorted);

% Plot the data
subplot(2,3,4);
hold on;
% scatter(x, y, 'b', 'filled'); % Scatter plot of data
plot(x_sorted, y_pred, 'k-', 'LineWidth', 2); % Regression line

% Plot confidence intervals
fill([x_sorted; flipud(x_sorted)], [y_ci(:,1); flipud(y_ci(:,2))], ...
    'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none'); % Confidence bounds
%% for figure in paper
a=1.1070; b=-148.6085; 
real_gain_range=[8.8695 93.6987];
median_gain_change=38.2600;

ix_move=find(ax>-.3 & ax<.2);
for i=1:size(Nuclear_FR,1)
    for ii=1:size(Nuclear_FR,2)
        CbN_move(i,ii)=sum(squeeze(Nuclear_FR(i,ii,ix_move)))/(length(ix_move)/length(ax));
    end
end

close all
clear X S C fr_all fr_move sfr_all sfr_move
bin=.001;
dax=[round((min(diffCCG)))-bin:bin:max(diffCCG)+bin];
dax=linspace(min(diffCCG),max(diffCCG),11);
for ii=1:length(dax)-1
    ix= find(diffCCG>=dax(ii) & diffCCG<dax(ii+1));
    if length(ix)>15
        X(ii)= mean(CbN_move(ix,2).*a-CbN_move(ix,1).*a);
        S(ii)=std(CbN_move(ix,2).*a-CbN_move(ix,1).*a)./sqrt(length(ix));
        C(ii)=mean(avg_ex(ix));
        fr_move(ii)= mean(CbN_move(ix,2)-CbN_move(ix,1));
        fr_all(ii)= mean(OVERALL_NUC_FR(ix,2)-OVERALL_NUC_FR(ix,1));

        sfr_move(ii)= std(CbN_move(ix,2)-CbN_move(ix,1))./sqrt(length(ix));
        sfr_all(ii)= std(OVERALL_NUC_FR(ix,2)-OVERALL_NUC_FR(ix,1))./sqrt(length(ix));

    else
        X(ii)=nan;
        S(ii)=nan;
        C(ii)=nan;
        fr_move(ii)=nan;
        fr_all(ii)=nan;
        sfr_move(ii)=nan;
        sfr_all(ii)=nan;
    end
end
subplot(2,3,1)
errorbar(dax(1:end-1),fr_all,sfr_all)
hold on
errorbar(dax(1:end-1),fr_move,sfr_move)
scatter(dax(1:end-1),fr_all,'ok')
scatter(dax(1:end-1),fr_move,'ok')
hold off



clear X S C fr_all fr_move
dax=linspace(min(avg_ex),max(avg_ex),16);
for ii=1:length(dax)-1
    ix= find(avg_ex>=dax(ii) & avg_ex<dax(ii+1));
    if length(ix)>15
        X(ii)= mean(decelGain(ix));
        S(ii)=std(decelGain(ix))./sqrt(length(ix));
        C(ii)=mean(num_of_MF(ix));
        fr_move(ii)= mean(CbN_move(ix,2)-CbN_move(ix,1));
        fr_all(ii)= mean(OVERALL_NUC_FR(ix,2)-OVERALL_NUC_FR(ix,1));
    else
        X(ii)=nan;
        S(ii)=nan;
        C(ii)=nan;
        fr_move(ii)=nan;
        fr_all(ii)=nan;
    end
end

subplot(2,3,2)
histogram(CbN_move(:,2).*a-CbN_move(:,1).*a,'BinWidth',10);
hold on
aa=axis;
plot(real_gain_range,[aa(4)-1 aa(4)-1],'k','LineWidth',2)
scatter(median_gain_change,aa(4),50,'or')

clear X S C fr_all fr_move sfr_all sfr_move
bin=.001;
dax=[round((min(diffCCG)))-bin:bin:max(diffCCG)+bin];
dax=linspace(min(diffCCG),max(diffCCG),11);
for ii=1:length(dax)-1
    ix= find(diffCCG>=dax(ii) & diffCCG<dax(ii+1));
    if length(ix)>15
        X(ii)= mean(CbN_move(ix,2).*a-CbN_move(ix,1).*a);
        S(ii)=std(CbN_move(ix,2).*a-CbN_move(ix,1).*a)./sqrt(length(ix));
        C(ii)=mean(avg_ex(ix));
        fr_move(ii)= mean(CbN_move(ix,2)-CbN_move(ix,1));
        fr_all(ii)= mean(OVERALL_NUC_FR(ix,2)-OVERALL_NUC_FR(ix,1));

        sfr_move(ii)= std(CbN_move(ix,2)-CbN_move(ix,1))./sqrt(length(ix));
        sfr_all(ii)= std(OVERALL_NUC_FR(ix,2)-OVERALL_NUC_FR(ix,1))./sqrt(length(ix));

    else
        X(ii)=nan;
        S(ii)=nan;
        C(ii)=nan;
        fr_move(ii)=nan;
        fr_all(ii)=nan;
                sfr_move(ii)=nan;
        sfr_all(ii)=nan;
    end
end
subplot 233
errorbar(dax(1:end-1),X,S)
hold on
scatter(dax(1:end-1),X,50,C,'fill')
hold off
mdl=fitlm(dax(1:end-1),X);
%
[x_sorted, idx] =sort([0:0.005:max(this_xcorr)]');% sort(this_xcorr(:)); % Sort x for smooth plotting
[y_pred, y_ci] = predict(mdl, x_sorted);

% Plot the data
hold on;
% scatter(x, y, 'b', 'filled'); % Scatter plot of data
plot(x_sorted, y_pred, 'k-', 'LineWidth', 2); % Regression line

% Plot confidence intervals
fill([x_sorted; flipud(x_sorted)], [y_ci(:,1); flipud(y_ci(:,2))], ...
    'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none'); % Confidence bounds