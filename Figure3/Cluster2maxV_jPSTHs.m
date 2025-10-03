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
% % % %                 endpoint_time=ReachS(i).out(ix,1);
% % % % %         endpoint_time=ReachS(i).out(end,1);
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
% % % %                 [h,t]=findpeaks(xV);
% % % %                 maxV(c_trials,1)=max(xV);%h(find(bhv_ax(t)>0,1,'first'));
% % % % 
% % % %                 Chs{c_this}= cellData(cc).Channels;
% % % %                 all_gain(c_this)=cellData(cc).gain;
% % % %             end
% % % %         end
% % % %     end
% % % %     all_gain= AllData(j).Trial(1).clusters;
% % % % %     xDec=maxV;
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
% % % %                 [r_real(c_all),p_real(c_all)]= corr(sim_spikes,maxV(:));
% % % % 
% % % %                 [~,ix] = sort(sim_spikes);
% % % %                 low_ix=ix(1:round(length(ix)*prob));
% % % %                 high_ix=ix(end-round(length(ix)*prob):end);
% % % %                 dec_lowSI(1,c_all)=nanmean(maxV(low_ix));
% % % %                 dec_highSI(1,c_all)=nanmean(maxV(high_ix));
% % % % 
% % % %                 [~,ix]= sort((maxV));
% % % %                 Pr_S1= cell1(ix(1:round(length(ix)*prob)),:);
% % % %                 Pr_S2= cell2(ix(1:round(length(ix)*prob)),:);
% % % %                 [raw, shift_predict, pred, surprise, std, ~, ~] = my_JPSTH(Pr_S1', Pr_S2', nbin,1);
% % % %                 all_sync_slow(c_all,:,:)=(raw-pred)./size(Pr_S1,1);%std;
% % % %                 all_fr_slow(c_all,:)=nanmean([nanmean(Pr_S1,1); nanmean(Pr_S2,1)]);
% % % %                 Pr_S1= cell1(ix(end-round(length(ix)*prob):end),:);
% % % %                 Pr_S2= cell2(ix(end-round(length(ix)*prob):end),:);
% % % %                 [raw, shift_predict, pred, surprise, std, ~, ~] = my_JPSTH(Pr_S1', Pr_S2', nbin,1);
% % % %                 all_sync_fast(c_all,:,:)=(raw-pred)./size(Pr_S1,1);%std;
% % % %                 all_fr_fast(c_all,:)=nanmean([nanmean(Pr_S1,1); nanmean(Pr_S2,1)]);
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
% % % %     clear FR FR2gain this_S nX nY nZ xVelocity yVelocity Speed xDec yDec
% % % %     clear time_min time_max minV maxV Chs
% % % % end
%%
load('E:\PAPER\Data2load\jPSTH_data\jPSTH2PeakV.mat');
%%
pair2take=find(gain(:,1) & gain(:,2)==1);
%%
subplot(2,3,1)
pcolor(ax,ax,imgaussfilt(squeeze(nanmean(all_sync_fast(pair2take,:,:),1)),3))
title('fast');
shading flat
colormap jet
axis([-.2 .1 -.2 .1]);
xlabel('Time to peak V (s)')
ylabel('Time to peak V (s)')
set(gca,'TickDir','out'); box off

subplot(2,3,2)
pcolor(ax,ax,imgaussfilt(squeeze(nanmean(all_sync_slow(pair2take,:,:),1)),3))
title('slow');
shading flat
colormap jet
axis([-.2 .1 -.2 .1]);
xlabel('Time to peak V (s)')
ylabel('Time to peak V (s)')
set(gca,'TickDir','out'); box off

subplot(2,3,3)
pcolor(ax,ax,imgaussfilt(squeeze(nanmean(all_sync_slow(pair2take,:,:)-all_sync_fast(pair2take,:,:),1)),3))
title('slow- fast');
shading flat
colormap jet
axis([-.2 .1 -.2 .1]);
caxis([-200 200])
xlabel('Time to peak V (s)')
ylabel('Time to peak V (s)')
set(gca,'TickDir','out'); box off
%%
sm=21;
gwin=gausswin(sm);
gwin=gwin./sum(gwin);
for i=1:length(pair2take)
    diag_slow(i,:)= conv2(diag(squeeze(all_sync_slow(pair2take(i),:,:))),gwin,'same');
    diag_fast(i,:)= conv2(diag(squeeze(all_sync_fast(pair2take(i),:,:))),gwin,'same');
end
%%
ix= find(max(abs(diag_fast),[],2)>1);
% btime= -.2;
% bix= find(ax>btime & ax<0);
bix= find(ax>-.05 & ax<.1);

subplot(2,3,4)
m=nanmean(diag_slow,1);
s=nanstd(diag_slow,[],1)./sqrt(size(diag_slow,1));
patch([ax fliplr(ax)],[m-s fliplr(m+s)],[1 .7 1],'EdgeColor','none')
% errorbar(ax,nanmean(diag_hypo,1),nanstd(diag_hypo,[],1)./sqrt(size(diag_hypo,1)),'Color','m')
hold on
m=nanmean(diag_fast,1);
s=nanstd(diag_fast,[],1)./sqrt(size(diag_fast,1));
patch([ax fliplr(ax)],[m-s fliplr(m+s)],[.5 1 .5],'EdgeColor','none')
plot(ax,nanmean(diag_slow,1),'m','LineWidth',2)
plot(ax,nanmean(diag_fast,1),'g','LineWidth',2)
% errorbar(ax,nanmean(diag_hyper,1),nanstd(diag_hyper,[],1)./sqrt(size(diag_hyper,1)),'Color','g')
a=axis;
plot([0 0],a(3:4),'--k')
plot(a(1:2),[0 0],'--k')
plot([ax(bix(1)) ax(bix(1)) ax(bix(end)) ax(bix(end)) ax(bix(1))],[a(3) a(4) a(4) a(3) a(3)],'Color',[.5 .5 .5]);
hold off
xlim([-.2 .1]);
xlabel('Time to peak V (s)')
ylabel('Cofiring (sp^2/s^2)')
set(gca,'TickDir','out'); box off

subplot(2,3,5)
histogram(nanmean(diag_slow(:,bix),2),'BinWidth',20,'FaceColor','m');
hold on
histogram(nanmean(diag_fast(:,bix),2),'BinWidth',20,'FaceColor','g');
a=axis;
plot([nanmean(nanmean(diag_slow(:,bix),2)) nanmean(nanmean(diag_slow(:,bix),2))],a(3:4),'m')
plot([nanmean(nanmean(diag_fast(:,bix),2)) nanmean(nanmean(diag_fast(:,bix),2))],a(3:4),'g')
hold off
[~,p]= ttest(nanmean(diag_slow(:,bix),2),nanmean(diag_fast(:,bix),2));
title(p);
xlabel('Cofiring (sp^2/s^2)')
ylabel('# of pairs')
set(gca,'TickDir','out'); box off

subplot(2,3,6)
plot(ax,smooth(nanmean(all_fr_slow(pair2take,:),1),sm),'g')
hold on
plot(ax,smooth(nanmean(all_fr_slow(pair2take,:),1)-nanstd(all_fr_slow(pair2take,:),[],1)./sqrt(length(pair2take)),sm),'--g')
plot(ax,smooth(nanmean(all_fr_slow(pair2take,:),1)+nanstd(all_fr_slow(pair2take,:),[],1)./sqrt(length(pair2take)),sm),'--g')
plot(ax,smooth(nanmean(all_fr_fast(pair2take,:),1),sm),'m')
plot(ax,smooth(nanmean(all_fr_fast(pair2take,:),1)-nanstd(all_fr_fast(pair2take,:),[],1)./sqrt(length(pair2take)),sm),'--m')
plot(ax,smooth(nanmean(all_fr_fast(pair2take,:),1)+nanstd(all_fr_fast(pair2take,:),[],1)./sqrt(length(pair2take)),sm),'--m')
a=axis;
plot([0 0],a(3:4),'--k')
hold off
xlim([-.2 .1])
xlabel('Time to peak V (s)')
ylabel('firing rate (sp/s)')
set(gca,'TickDir','out'); box off