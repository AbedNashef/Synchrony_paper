%
clear
udir= '../new_data/';
dates = dir([udir 'D*']);
bin2use=1;
addpath ../helper_functions/
tbf=1; taf=.5;
win2corr=.2; d2corr= .020; ax2corr=-tbf:d2corr:taf-win2corr;
ax=-tbf*1000:bin2use:taf*1000; ax=ax(1:end-1); ax=ax./1000;
c_all=0;
c_allv=0;
c_alla=0;
n=50;
nbin=1;
class2take={'PC'};
freq= 120;
L= ((taf+tbf)*freq);
move_thresh=0.3;
prob=1/4;
%%
nax= [-0.5:1/1000:0.5];nax=nax(1:end-1);
t1= find(nax<-0.25); t2= find(nax>0 & nax<.25);
bhv_ax=[-tbf:1/freq:taf];
all_acc=[];
all_vf=[];
all_v0=[];
all_dec=[];
all_SI=[];
all_mxV=[];
all_endpoint=[];
all_cell1=[];
all_cell2=[];
date_id=[];
gain1=[]; gain2=[];
cell1_id=[];
cell2_id=[];
all_mice_id=[];
% [ res shift_predict psth_pred surprise_mat std_mat NUM_TRIAL pval_mat] = my_JPSTH(cut1, cut2, BIN_SIZE,gauss_filt_std)
%%
for j=1:length(dates)
    load([udir dates(j).name]);
    c_trials=0;
    index=strfind(dates(j).name,'_');
    mouse_id=dates(j).name(index(end-1)+1:index(end)+1);
    for i=1:length(ReachS)
        times = ReachS(i).filt_kin(:,1);
        ix= find(ReachS(i).out(:,2)>move_thresh,1,'first');
        endpoint_time=ReachS(i).out(ix,1);
        endpoint_time=ReachS(i).out(end,1);
        if isempty(endpoint_time)
            ix=find(ReachS(i).filt_kin(:,2)>move_thresh,1,'first');
            endpoint_time=ReachS(i).filt_kin(ix,1);
        end
        times= ReachS(i).filt_kin(:,1);
        P1= repmat(ReachS(i).filt_kin(1,2:4),size(ReachS(i).filt_kin,1),1);
        tt=times(find(times>=endpoint_time-tbf & times<=endpoint_time+taf));
        
        ix= find(ReachS(i).filt_kin(:,1)>=endpoint_time-tbf & ReachS(i).filt_kin(:,1)<=endpoint_time+taf);
        dt=median(diff(ReachS(i).filt_kin(ix,1)));
        tq=-tbf:1/freq:taf;
        t=ReachS(i).filt_kin(ix,1)-endpoint_time;
        xV=interp1(t,ReachS(i).filt_kin(ix,6),tq);
        yV=interp1(t,ReachS(i).filt_kin(ix,7),tq);
        V=interp1(t,ReachS(i).filt_kin(ix,5),tq);
        tt=times(find(times>=endpoint_time-tbf & times<=endpoint_time+taf));
        
        vStimMode = isfield(ReachS(i),'stim');
        if vStimMode, vStim=ReachS(i).stim; else, vStim=0; end
        if isempty(vStim), vStim=0; end
        vExclude=ReachS(i).exclude;
        c_this=0;
        if ~vExclude & ~vStim
            c_trials=c_trials+1;
            for cc=1:length(cellData)
                c_this=c_this+1;
                this_trc=cellData(cc).Bin1;
                index= find(this_trc(:,1)>=endpoint_time-tbf & this_trc(:,1)<=endpoint_time+taf);
                FR(c_this,c_trials,:)= this_trc(index,2);
                
                endpoint_XY(c_trials,:)=[ReachS(i).out(end,2) ReachS(i).out(end,3)];
                endpoint_D(c_trials)= sqrt(ReachS(i).out(end,2)^2+ReachS(i).out(end,3)^2);
                
                T= ReachS(i).out(end,1);
                ix= find(ReachS(i).filt_kin(:,1)==T);
                index= find(ReachS(i).filt_kin(:,1)>=T-.3 & ReachS(i).filt_kin(:,1)<=T+.3);
                maxV(c_trials,:)= [max(ReachS(i).filt_kin(index,5)) max(ReachS(i).filt_kin(index,6)) max(ReachS(i).filt_kin(index,7)) max(ReachS(i).filt_kin(index,8))];
                all_gain(c_this)=cellData(cc).gain;
                
                [h,t]=findpeaks(xV);
                maxV(c_trials,1)=h(find(bhv_ax(t)>0,1,'first'));
                mix=t(find(bhv_ax(t)>0,1,'first'));
                time_max(c_trials,1)=bhv_ax(mix);
                h=find(islocalmin(xV));
                minix=h(find(bhv_ax(h)>time_max(c_trials,1),1,'first'));
                minV(c_trials,1)=xV(minix);
                time_min(c_trials,1)=bhv_ax(minix);
                xDec(c_trials)=(maxV(c_trials,1)-minV(c_trials,1))/(time_max(c_trials,1)-time_min(c_trials,1));
                
                [h,t]=findpeaks(yV);
                maxV(c_trials,2)=h(find(bhv_ax(t)>0,1,'first'));
                mix=t(find(bhv_ax(t)>0,1,'first'));
                time_max(c_trials,2)=bhv_ax(mix);
                h=find(islocalmin(yV));
                minix=h(find(bhv_ax(h)>time_max(c_trials,2),1,'first'));
                minV(c_trials,2)=yV(minix);
                time_min(c_trials,2)=bhv_ax(minix);
                yDec(c_trials)=(maxV(c_trials,2)-minV(c_trials,2))/(time_max(c_trials,2)-time_min(c_trials,2));
                
                V0(c_trials)=ReachS(i).out(1,6);
                tt= ReachS(i).out(end,1)-ReachS(i).out(:,1);
                try
                    Vf(c_trials)=interp1(tt,ReachS(i).out(:,6),.1);
                catch ME
                    Vf(c_trials)=nan;
                end
                [~,T]= max(ReachS(i).out(:,6));
                acc(c_trials)= (max(ReachS(i).out(:,6))-ReachS(i).out(1,6))./(ReachS(i).out(T,1)-ReachS(i).out(1,1));
                Chs{c_this}= cellData(cc).Channels;
            end
        end
    end
    endpoint_XY(:,1)=endpoint_XY(:,1);%./max(endpoint_XY(:,1));
    endpoint_XY(:,2)=endpoint_XY(:,2);%./max(endpoint_XY(:,2));
    endpoint_XY=endpoint_XY-repmat(nanmean(endpoint_XY,1),size(endpoint_XY,1),1);
    maxV=maxV./repmat(max(maxV,[],1),size(maxV,1),1);
    endpoint_D=endpoint_D-nanmean(endpoint_D);
    xDec=abs(xDec); xDec=xDec./max(xDec);
    yDec=abs(yDec); yDec=yDec./max(yDec);
    for i=1:size(FR,1)-1
        for ii=i+1:size(FR,1)
            ch1=Chs{i};
            ch2=Chs{ii};
            if ~isempty(ch1) & ~isempty(ch2)
                overlap= length(intersect(ch1,ch2))/length(unique([ch1,ch2]));
            else
                overlap=0;
            end
            if overlap==0 & (all_gain(i)<-.05 & all_gain(ii)<-.05)
                c_all=c_all+1;
                cell1= squeeze(FR(i,:,:));
                cell2= squeeze(FR(ii,:,:));
                
                tix= find(ax>-.3 & ax<.3);
                Pr_S1= nanmean(cell1(:,tix),1);
                Pr_S2= nanmean(cell2(:,tix),1);
                Pr_S1S2= nanmean(cell1(:,tix).*cell2(:,tix),1);
                sumSI= nansum((cell1(:,tix).*cell2(:,tix))./repmat(Pr_S1.*Pr_S2,size(all_gain,1),1),2);
                %             sumSI= max((cell1(:,tix).*cell2(:,tix))./repmat(Pr_S1.*Pr_S2,size(all_gain,1),1),[],2);
                %             tmp=(cell1(:,tix).*cell2(:,tix))./repmat(Pr_S1.*Pr_S2,size(all_gain,1),1);
                %             sumSI= tmp(:,find(ax==-.05));
                all_dec=[all_dec; xDec(:) yDec(:)];
                all_acc=[all_acc; acc(:)];
                all_vf=[all_vf; Vf(:)];
                all_v0=[all_v0; V0(:)];
                all_SI=[all_SI; sumSI./max(sumSI)];
                all_mxV=[all_mxV; maxV];
                all_endpoint=[all_endpoint; endpoint_XY endpoint_D(:)];
                tmp=(sum(cell1(:,tix),2)-sum(cell1(:,find(ax<-.5)),2));
                tmp=tmp;%./max(tmp);
                all_cell1=[all_cell1; tmp./max(abs(tmp))];
                tmp=(sum(cell2(:,tix),2)-sum(cell2(:,find(ax<-.5)),2));
                tmp=tmp;%./max(tmp);
                all_cell2=[all_cell2; tmp./max(abs(tmp))];
                date_id= [date_id; ones(size(sumSI)).*j];
                cell1_id=[cell1_id; ones(size(sumSI)).*i];
                cell2_id=[cell2_id; ones(size(sumSI)).*ii];
                
                gain1=[gain1; ones(length(sumSI),1).*all_gain(i)];
                gain2=[gain2; ones(length(sumSI),1).*all_gain(ii)];
                Dist4cells(c_all)= abs(cellData(i).depth-cellData(ii).depth);
                all_mice_id=[all_mice_id; repmat(mat2cell(mouse_id,1),length(sumSI),1)];
            end
        end
    end
    clear v0 vf acc FR FR2gain this_S nX nY nZ endpoint_XY endpoint_D maxV xDec yDec Chs
end
%%
all_mod= nanmean([all_cell1 all_cell2],2);
ix= find(gain1<-.05 & gain2<-.05);
pausers_SI= all_SI(ix);
pausers_mod= all_mod(ix);
pausers_SI=pausers_SI;%./max(pausers_SI);
pausers_mod=pausers_mod;%./max(pausers_mod);
pausers_endpoint= all_endpoint(ix,:);
pausers_V= all_mxV(ix);
pausers_dec= all_dec(ix,:);
pausers_days= date_id(ix);
pausers_cell1= cell1_id(ix);
pausers_cell2= cell2_id(ix);
pausers_mouse= all_mice_id(ix);
pausers_acc= all_acc(ix);
pausers_vf= all_vf(ix);
pausers_v0= all_v0(ix);
%%
top_qut=.99; bot_qut=.05;
nbin=10;
mod_lim= [quantile(pausers_mod,bot_qut) quantile(pausers_mod,top_qut)];
si_lim= [quantile(pausers_SI,bot_qut) quantile(pausers_SI,top_qut)];
dmod=mod_lim(1):range(mod_lim)/nbin:mod_lim(2);
dsi=si_lim(1):range(si_lim)/nbin:si_lim(2);
avg_v=[]; avg_end=[]; avg_dec=[]; avg_vf=[];  avg_acc=[];  avg_v0=[];
for i=1:length(dmod)-1
    ix= find(pausers_mod>=dmod(i) & pausers_mod<=dmod(i+1));
    avg_v(i)= nanmean(pausers_V(ix));
    avg_end(i)= nanmean(pausers_endpoint(ix));
    std_end(i)= nanstd(pausers_endpoint(ix))./sqrt(length(ix));
    avg_dec(i)= nanmean(pausers_dec(ix));
    
    avg_vf(i)=nanmean(pausers_vf(ix));
    avg_v0(i)=nanmean(pausers_v0(ix));
    avg_acc(i)=nanmean(pausers_acc(ix));
    
end

%%
si_avg_v=[]; si_avg_end=[]; si_avg_dec=[];
for i=1:length(dsi)-1
    ix= find(pausers_SI>=dsi(i) & pausers_SI<=dsi(i+1));
    si_avg_v(i)= nanmean(pausers_V(ix));
    si_avg_end(i)= nanmean(pausers_endpoint(ix));
    si_std_end(i)= nanstd(pausers_endpoint(ix))./sqrt(length(ix));
    si_avg_dec(i)= nanmean(pausers_dec(ix));
end
%%
simulated_end=(pausers_vf(:).^2-pausers_v0(:).^2)./(2*pausers_acc(:));
simulated_end=simulated_end-nanmean(simulated_end);
nbin=10;
avg_sim_end=[]; si_avg_sim_end=[];
for i=1:length(dmod)-1
    ix= find(pausers_mod>=dmod(i) & pausers_mod<=dmod(i+1));
    avg_sim_end(i)=nanmean(simulated_end(ix));
    std_sim_end(i)=nanstd(simulated_end(ix))./sqrt(length(ix));
    ix= find(pausers_SI>=dsi(i) & pausers_SI<=dsi(i+1));
    si_avg_sim_end(i)=nanmean(simulated_end(ix));
    si_std_sim_end(i)=nanstd(simulated_end(ix))./sqrt(length(ix));
end

%%
sim_end4both=[]; real_end4both=[];
for i=1:length(dsi)-1
    for ii=1:length(dmod)-1
        ix= find(pausers_mod>=dmod(ii) & pausers_mod<=dmod(ii+1) & pausers_SI>=dsi(i) & pausers_SI<=dsi(i+1));
        sim_end4both(i,ii)=nanmean(simulated_end(ix));
        real_end4both(i,ii)= nanmean(pausers_endpoint(ix));
    end
end
subplot(2,2,1)
pcolor(dmod(1:end-1),dsi(1:end-1),imgaussfilt(sim_end4both,1))
shading flat
colormap(flipud(brewermap([],'RdBu')))
xlabel('Norm. Mod');
ylabel('Norm. SI');

subplot(2,2,2)
pcolor(dmod(1:end-1),dsi(1:end-1),imgaussfilt(real_end4both,1))
shading flat
colormap(flipud(brewermap([],'RdBu')))
xlabel('Norm. Mod');
ylabel('Norm. SI');

subplot(2,2,3)
% scatter(dmod(1:end-1),smooth(avg_sim_end,3),40,'or','fill')
errorbar(dmod(1:end-1),smooth(nanmean(sim_end4both,1)),smooth(nanstd(sim_end4both,[],1))./sqrt(size(sim_end4both,1)),'r')
hold on
errorbar(dmod(1:end-1),smooth(nanmean(real_end4both,1)),smooth(nanstd(real_end4both,[],1))./sqrt(size(real_end4both,1)),'k')
% hold on
% scatter(dmod(1:end-1),smooth(avg_end,3),40,'ok','fill')
% errorbar(dmod(1:end-1),smooth(avg_sim_end,3),std_sim_end,'Color','r')
% errorbar(dmod(1:end-1),smooth(avg_end,3),std_end,'Color','k')
hold off
xlabel('Norm. Mod');
ylabel('\DeltaEnd');

subplot(2,2,4)
errorbar(dsi(1:end-1),smooth(nanmean(sim_end4both,2)),smooth(nanstd(sim_end4both,[],2))./sqrt(size(sim_end4both,1)),'r')
hold on
errorbar(dsi(1:end-1),smooth(nanmean(real_end4both,2)),smooth(nanstd(real_end4both,[],2))./sqrt(size(real_end4both,1)),'k')

% scatter(dsi(1:end-1),smooth(si_avg_sim_end,3),40,'or','fill')
% hold on
% scatter(dsi(1:end-1),smooth(si_avg_end,3),40,'ok','fill')
% errorbar(dsi(1:end-1),smooth(si_avg_sim_end,3),si_std_sim_end,'Color','r')
% errorbar(dsi(1:end-1),smooth(si_avg_end,3),si_std_end,'Color','k')
hold off
xlabel('Norm. SI');
ylabel('\DeltaEnd');
