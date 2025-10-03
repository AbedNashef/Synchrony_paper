%
clear
% % % % udir= 'E:/all_data/';
% % % % dates = dir([udir 'D*']);
% % % % bin2use=1;
% % % % addpath ../helper_functions/
% % % % tbf=0.5; taf=.5;
% % % % win2corr=.2; d2corr= .020; ax2corr=-tbf:d2corr:taf-win2corr;
% % % % c_all=0;
% % % % c_allv=0;
% % % % c_alla=0;
% % % % n=50;
% % % % nbin=1;
% % % % class2take={'PC'};
% % % % freq= 120;
% % % % L= ((taf+tbf)*freq);
% % % % move_thresh=0.2;
% % % % prob=1/4;
% % % % load C:/Users/PersonLab/Desktop/'synchrony paper'/Dylans_data/LASSO/FullData_withAutoClusters.mat;
% % % % ax=-tbf*1000:bin2use:taf*1000; ax=ax(1:end-1); ax=ax./1000;
% % % % %%
% % % % nax= [-tbf:1/1000:taf];nax=nax(1:end-1);
% % % % t1= find(nax<-0.25); t2= find(nax>0 & nax<.25);
% % % % bhv_ax=[-tbf:1/freq:taf];
% % % % 
% % % % % [ res shift_predict psth_pred surprise_mat std_mat NUM_TRIAL pval_mat] = my_JPSTH(cut1, cut2, BIN_SIZE,gauss_filt_std)
% % % % %%
% % % % for j=1:length(dates)
% % % %     load([udir dates(j).name]);
% % % %     c_trials=0;
% % % %     index=strfind(dates(j).name,'_');
% % % %     mouse_id=dates(j).name(index(end-1)+1:index(end)+1);
% % % %     for i=1:length(ReachS)
% % % %         times = ReachS(i).filt_kin(:,1);
% % % %         ix= find(ReachS(i).out(:,2)>move_thresh,1,'first');
% % % %         [~,ix]= max(ReachS(i).out(:,6));
% % % % 
% % % %         %         endpoint_time=ReachS(i).out(ix,1);
% % % %         endpoint_time=ReachS(i).out(end,1);
% % % %         if isempty(endpoint_time)
% % % %             ix=find(ReachS(i).filt_kin(:,2)>move_thresh,1,'first');
% % % %             endpoint_time=ReachS(i).filt_kin(ix,1);
% % % %         end
% % % %         ix= find(ReachS(i).filt_kin(:,1)>=endpoint_time-tbf & ReachS(i).filt_kin(:,1)<=endpoint_time+taf);
% % % %         dt=median(diff(ReachS(i).filt_kin(ix,1)));
% % % %         tq=-tbf:1/freq:taf;
% % % %         t=ReachS(i).filt_kin(ix,1)-endpoint_time;
% % % %         xV=interp1(ReachS(i).filt_kin(:,1),ReachS(i).filt_kin(:,6),tq+endpoint_time);
% % % %         yV=interp1(ReachS(i).filt_kin(:,1),ReachS(i).filt_kin(:,7),tq+endpoint_time);
% % % %         V=interp1(ReachS(i).filt_kin(:,1),ReachS(i).filt_kin(:,5),tq+endpoint_time);
% % % %         if ~isempty(find(xV>50))
% % % %             ix=find(xV<50);
% % % %             V=V(ix);
% % % %             yV=yV(ix);
% % % %             xV=xV(ix);
% % % %             V=interp1(tq(ix)+endpoint_time,V,tq+endpoint_time);
% % % %             yV=interp1(tq(ix)+endpoint_time,yV,tq+endpoint_time);
% % % %             xV=interp1(tq(ix)+endpoint_time,xV,tq+endpoint_time);
% % % %         end
% % % %         times= ReachS(i).filt_kin(:,1);
% % % %         P1= repmat(ReachS(i).filt_kin(1,2:4),size(ReachS(i).filt_kin,1),1);
% % % %         tt=times(find(times>=endpoint_time-tbf & times<=endpoint_time+taf));
% % % % 
% % % %         vStimMode = isfield(ReachS(i),'stim');
% % % %         if vStimMode, vStim=ReachS(i).stim; else, vStim=0; end
% % % %         if isempty(vStim), vStim=0; end
% % % %         vExclude=ReachS(i).exclude;
% % % %         c_this=0;
% % % %         if ~vExclude & ~vStim
% % % %             c_trials=c_trials+1;
% % % %             for cc=1:length(cellData)
% % % %                 c_this=c_this+1;
% % % %                 this_trc=cellData(cc).Bin1;
% % % %                 index= find(this_trc(:,1)>=endpoint_time-tbf & this_trc(:,1)<=endpoint_time+taf);
% % % %                 FR(c_this,c_trials,:)= this_trc(index,2);
% % % % 
% % % %                 index= find(this_trc(:,1)>=endpoint_time-0.5 & this_trc(:,1)<=endpoint_time+0.5);
% % % %                 FR2gain(c_this,c_trials,:)= this_trc(index,2);
% % % %                 %                 [h,t]=findpeaks(xV);
% % % %                 %                 maxV(c_trials,1)=h(find(bhv_ax(t)>0,1,'first'));
% % % %                 %                 mix=t(find(bhv_ax(t)>0,1,'first'));
% % % %                 %                 time_max(c_trials,1)=bhv_ax(mix);
% % % %                 %                 h=find(islocalmin(xV));
% % % %                 %                 minix=h(find(bhv_ax(h)>time_max(c_trials,1),1,'first'));
% % % %                 %                 minV(c_trials,1)=xV(minix);
% % % %                 %                 time_min(c_trials,1)=bhv_ax(minix);
% % % %                 %                 xDec(c_trials)=(maxV(c_trials,1)-minV(c_trials,1))/(time_max(c_trials,1)-time_min(c_trials,1));
% % % % 
% % % %                 endpoints(c_trials)=ReachS(i).out(end,2);
% % % %                 [h,t]=findpeaks(xV);
% % % %                 [~,mix]= max(xV(1:find(bhv_ax==0)));
% % % %                 maxV(c_trials,1)=xV(mix);
% % % %                 %                     mix=t(find(bhv_ax(t)>-.10,1,'first'));
% % % %                 time_max(c_trials,1)=bhv_ax(mix);
% % % %                 h=find(islocalmin(xV));
% % % %                 minix=h(find(bhv_ax(h)>time_max(c_trials,1),1,'first'));
% % % %                 minV(c_trials,1)=xV(minix);
% % % %                 time_min(c_trials,1)=bhv_ax(minix);
% % % %                 xDec(c_trials)=(maxV(c_trials,1)-minV(c_trials,1))/(time_max(c_trials,1)-time_min(c_trials,1));
% % % % 
% % % %                 [h,t]=findpeaks(yV);
% % % %                 maxV(c_trials,2)=h(find(bhv_ax(t)>0,1,'first'));
% % % %                 mix=t(find(bhv_ax(t)>0,1,'first'));
% % % %                 time_max(c_trials,2)=bhv_ax(mix);
% % % %                 h=find(islocalmin(yV));
% % % %                 minix=h(find(bhv_ax(h)>time_max(c_trials,2),1,'first'));
% % % %                 minV(c_trials,2)=yV(minix);
% % % %                 time_min(c_trials,2)=bhv_ax(minix);
% % % %                 yDec(c_trials)=(maxV(c_trials,2)-minV(c_trials,2))/(time_max(c_trials,2)-time_min(c_trials,2));
% % % %                 Chs{c_this}= cellData(cc).Channels;
% % % %                 all_gain(c_this)=cellData(cc).gain;
% % % %             end
% % % %         end
% % % %     end
% % % %     all_gain= AllData(j).Trial(1).clusters;
% % % %     overall_dec=sqrt(xDec.^2+yDec.^2);
% % % %     xDec=-1.*(xDec);
% % % %     %     xDec=abs(xDec);
% % % %     endpoints=endpoints-nanmean(endpoints);
% % % %     for i=1:size(FR,1)-1
% % % %         for ii=i+1:size(FR,1)
% % % %             ch1=Chs{i};
% % % %             ch2=Chs{ii};
% % % %             if ~isempty(ch1) & ~isempty(ch2)
% % % %                 overlap= length(intersect(ch1,ch2))/length(unique([ch1,ch2]));
% % % %             else
% % % %                 overlap=0;
% % % %             end
% % % %             if overlap==0 & ((all_gain(i)<=2 & all_gain(ii)<=2) | (all_gain(i)<=2 & all_gain(ii)<=2))
% % % %                 c_all=c_all+1;
% % % %                 cell1= squeeze(FR(i,:,:));
% % % %                 cell2= squeeze(FR(ii,:,:));
% % % %                 tix= find(ax>-.3 & ax<.0);
% % % % 
% % % %                 %                 sim_spikes=nansum((cell1(:,tix).*cell2(:,tix)),2);
% % % %                 sim_spikes=nanmean((cell1(:,tix).*cell2(:,tix)),2)-nanmean(nanmean((cell1(:,tix).*cell2(:,tix)),2));
% % % % 
% % % %                 [r_real(c_all),p_real(c_all)]= corr(sim_spikes,xDec(:));
% % % %                 h=fitlm(maxV(:,1),xDec);
% % % %                 se= h.Residuals.Raw;
% % % % 
% % % %                 [~,ix] = sort(sim_spikes);
% % % %                 low_ix=ix(1:round(length(ix)*prob));
% % % %                 high_ix=ix(end-round(length(ix)*prob):end);
% % % %                 dec_lowSI(1,c_all)=nanmean(xDec(low_ix));
% % % %                 dec_lowSI(2,c_all)=nanmean(yDec(low_ix));
% % % %                 dec_highSI(1,c_all)=nanmean(xDec(high_ix));
% % % %                 dec_highSI(2,c_all)=nanmean(yDec(high_ix));
% % % % 
% % % %                 [~,ix]= sort((se));
% % % % %                 [~,ix]=sort(overall_dec);
% % % %                 Pr_S1= cell1(ix(1:round(length(ix)*prob)),:);
% % % %                 Pr_S2= cell2(ix(1:round(length(ix)*prob)),:);
% % % %                 [raw, shift_predict, pred, surprise, std, ~, ~] = my_JPSTH(Pr_S1', Pr_S2', nbin,1);
% % % %                 all_sync_slow(c_all,:,:)=(raw-pred)./size(Pr_S1,1);%std;
% % % %                 all_fr_slow(c_all,:)=nanmean([nanmean(Pr_S1,1); nanmean(Pr_S2,1)]);
% % % % 
% % % %                 Pr_S1= cell1(ix(end-round(length(ix)*prob):end),:);
% % % %                 Pr_S2= cell2(ix(end-round(length(ix)*prob):end),:);
% % % %                 [raw, shift_predict, pred, surprise, std, ~, ~] = my_JPSTH(Pr_S1', Pr_S2', nbin,1);
% % % %                 all_sync_fast(c_all,:,:)=(raw-pred)./size(Pr_S1,1);%std;
% % % %                 all_fr_fast(c_all,:)=nanmean([nanmean(Pr_S1,1); nanmean(Pr_S2,1)]);
% % % %                 time_max_all(c_all,1)=mean(time_max(ix(1:round(length(ix)*prob)),1));
% % % %                 time_max_all(c_all,2)=mean(time_max(ix(end-round(length(ix)*prob):end),1));
% % % %                 maxV_all(c_all,1)=mean(maxV(ix(1:round(length(ix)*prob)),1));
% % % %                 maxV_all(c_all,2)=mean(maxV(ix(end-round(length(ix)*prob):end),1));
% % % %                 end_all(c_all,1)=mean(endpoints(ix(1:round(length(ix)*prob))));
% % % %                 end_all(c_all,2)=mean(endpoints(ix(end-round(length(ix)*prob):end)));
% % % %                 dec_all(c_all,1)=mean(xDec(ix(1:round(length(ix)*prob))));
% % % %                 dec_all(c_all,2)=mean(xDec(ix(end-round(length(ix)*prob):end)));
% % % % 
% % % %                 gain(c_all,1)=all_gain(i);
% % % %                 gain(c_all,2)=all_gain(ii);
% % % %                 date_id(c_all)= j;
% % % %                 Dist4cells(c_all)= abs(cellData(i).depth-cellData(ii).depth);
% % % %                 all_mice_id{c_all}=mouse_id;
% % % % 
% % % %             end
% % % %         end
% % % %     end
% % % %     clear FR FR2gain this_S nX nY nZ xVelocity yVelocity Speed xDec yDec endpoints
% % % %     clear time_min time_max minV maxV Chs
% % % % end
%%
load('E:\PAPER\Data2load\jPSTH_data\jPSTH2DecGain.mat');
%%
subplot(2,3,1)
pcolor(ax,ax,imgaussfilt(squeeze(nanmean(all_sync_fast,1)),3))
title('fast');
shading flat
colormap jet

subplot(2,3,2)
pcolor(ax,ax,imgaussfilt(squeeze(nanmean(all_sync_slow,1)),3))
title('slow');
shading flat
colormap jet

subplot(2,3,3)
pcolor(ax,ax,imgaussfilt(squeeze(nanmean(all_sync_slow-all_sync_fast,1)),3))
title('slow- fast');
shading flat
colormap jet

%%
sm=21;
gwin=gausswin(sm);
gwin=gwin./sum(gwin);
for i=1:size(all_sync_fast,1)
    diag_hypo(i,:)= conv2(diag(squeeze(all_sync_slow(i,:,:))),gwin,'same');
    diag_hyper(i,:)= conv2(diag(squeeze(all_sync_fast(i,:,:))),gwin,'same');
    %     diag_hypo(i,:)= smooth(diag(squeeze(all_sync_slow(i,:,:))),sm);
    % diag_hyper(i,:)= smooth(diag(squeeze(all_sync_fast(i,:,:))),sm);
end
ix= find(max(abs(diag_hyper),[],2)>1);
% btime= -.2;
% bix= find(ax>btime & ax<0);
bix= find(ax>-.15 & ax<.0);

subplot(2,3,4)
m=nanmean(diag_hypo,1);
s=nanstd(diag_hypo,[],1)./sqrt(size(diag_hypo,1));
patch([ax fliplr(ax)],[m-s fliplr(m+s)],[1 .7 1],'EdgeColor','none')
% errorbar(ax,nanmean(diag_hypo,1),nanstd(diag_hypo,[],1)./sqrt(size(diag_hypo,1)),'Color','m')
hold on
m=nanmean(diag_hyper,1);
s=nanstd(diag_hyper,[],1)./sqrt(size(diag_hyper,1));
patch([ax fliplr(ax)],[m-s fliplr(m+s)],[.5 1 .5],'EdgeColor','none')
plot(ax,nanmean(diag_hypo,1),'m','LineWidth',2)
plot(ax,nanmean(diag_hyper,1),'g','LineWidth',2)
% errorbar(ax,nanmean(diag_hyper,1),nanstd(diag_hyper,[],1)./sqrt(size(diag_hyper,1)),'Color','g')
a=axis;
plot([0 0],a(3:4),'--k')
plot(a(1:2),[0 0],'--k')
plot([ax(bix(1)) ax(bix(1)) ax(bix(end)) ax(bix(end)) ax(bix(1))],[a(3) a(4) a(4) a(3) a(3)],'Color',[.5 .5 .5]);
hold off

subplot(2,3,5)
histogram(nanmean(diag_hypo(:,bix),2),'BinWidth',.0200,'FaceColor','m');
hold on
histogram(nanmean(diag_hyper(:,bix),2),'BinWidth',.0200,'FaceColor','g');
a=axis;
plot([nanmean(nanmean(diag_hypo(:,bix),2)) nanmean(nanmean(diag_hypo(:,bix),2))],a(3:4),'m')
plot([nanmean(nanmean(diag_hyper(:,bix),2)) nanmean(nanmean(diag_hyper(:,bix),2))],a(3:4),'g')
hold off
[~,p]= ttest(nanmean(diag_hypo(:,bix),2),nanmean(diag_hyper(:,bix),2));
title(p);

subplot(2,3,6)
% histogram(r_real,'BinWidth',.02,'FaceColor','k','EdgeColor','none')
% hold on
% plot([nanmean(r_real) nanmean(r_real)],[0 25],'k','LineWidth',1)
% [~,p]= ttest(r_real,0)
% title(p)
% set(gca,'TickDir','out'); box off

[~,ix]=unique(maxV_all(:,1));
scatter(maxV_all(ix,1),maxV_all(ix,2),70,'ok','fill')
a=axis;
hold on
plot([min(a) max(a)],[min(a) max(a)],'-k')
hold off
[~,p]= ttest(maxV_all(ix,1),maxV_all(ix,2));
title(p)
set(gca,'TickDir','out'); box off
xlabel('max V; slow')
ylabel('max V; fast')