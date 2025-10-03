%
clear
% % % % % % udir= 'E:/all_data/';
% % % % % % dates = dir([udir 'D*']);
% % % % % % bin2use=1;
% % % % % % addpath ../helper_functions/
% % % % % % tbf=0.5; taf=.5;
% % % % % % win2corr=.2; d2corr= .020; ax2corr=-tbf:d2corr:taf-win2corr;
% % % % % % c_all=0;
% % % % % % c_allv=0;
% % % % % % c_alla=0;
% % % % % % n=50;
% % % % % % nbin=1;
% % % % % % class2take={'PC'};
% % % % % % freq= 120;
% % % % % % L= ((taf+tbf)*freq);
% % % % % % move_thresh=0.2;
% % % % % % prob=1/4;
% % % % % % load C:/Users/PersonLab/Desktop/'synchrony paper'/Dylans_data/LASSO/FullData_withAutoClusters.mat;
% % % % % % ax=-tbf*1000:bin2use:taf*1000; ax=ax(1:end-1); ax=ax./1000;
% % % % % % 
% % % % % % %%
% % % % % % nax= [-0.5:1/1000:0.5];nax=nax(1:end-1);
% % % % % % t1= find(nax<-0.25); t2= find(nax>0 & nax<.25);
% % % % % % 
% % % % % % % [ res shift_predict psth_pred surprise_mat std_mat NUM_TRIAL pval_mat] = my_JPSTH(cut1, cut2, BIN_SIZE,gauss_filt_std)
% % % % % % %%
% % % % % % for j=1:length(dates)
% % % % % %     load([udir dates(j).name]);
% % % % % %     c_trials=0;
% % % % % %     index=strfind(dates(j).name,'_');
% % % % % %     mouse_id=dates(j).name(index(end-1)+1:index(end)+1);
% % % % % %     for i=1:length(ReachS)
% % % % % %         times = ReachS(i).filt_kin(:,1);
% % % % % %         ix= find(ReachS(i).out(:,2)>move_thresh,1,'first');
% % % % % %         [~,ix]= max(ReachS(i).out(:,6));
% % % % % %         endpoint_time=ReachS(i).out(ix,1);
% % % % % %         endpoint_time=ReachS(i).out(end,1);
% % % % % %         if isempty(endpoint_time)
% % % % % %             ix=find(ReachS(i).filt_kin(:,2)>move_thresh,1,'first');
% % % % % %             endpoint_time=ReachS(i).filt_kin(ix,1);
% % % % % %         end
% % % % % %         times= ReachS(i).filt_kin(:,1);
% % % % % %         P1= repmat(ReachS(i).filt_kin(1,2:4),size(ReachS(i).filt_kin,1),1);
% % % % % %         tt=times(find(times>=endpoint_time-tbf & times<=endpoint_time+taf));
% % % % % % 
% % % % % %         vStimMode = isfield(ReachS(i),'stim');
% % % % % %         if vStimMode, vStim=ReachS(i).stim; else, vStim=0; end
% % % % % %         if isempty(vStim), vStim=0; end
% % % % % %         vExclude=ReachS(i).exclude;
% % % % % %         c_this=0;
% % % % % %         if ~vExclude & ~vStim
% % % % % %             c_trials=c_trials+1;
% % % % % %             for cc=1:length(cellData)
% % % % % %                 c_this=c_this+1;
% % % % % %                 this_trc=cellData(cc).Bin1;
% % % % % %                 index= find(this_trc(:,1)>=endpoint_time-tbf & this_trc(:,1)<=endpoint_time+taf);
% % % % % %                 FR(c_this,c_trials,:)= this_trc(index,2);
% % % % % % 
% % % % % %                 index= find(this_trc(:,1)>=endpoint_time-0.5 & this_trc(:,1)<=endpoint_time+0.5);
% % % % % %                 FR2gain(c_this,c_trials,:)= this_trc(index,2);
% % % % % %                 endpoint_XY(c_trials,:)=[ReachS(i).out(end,2) ReachS(i).out(end,3)];
% % % % % %                 maxV(c_trials)=max(ReachS(i).out(:,6));
% % % % % %                 Chs{c_this}= cellData(cc).Channels;
% % % % % %                 all_gain(c_this)=cellData(cc).gain;
% % % % % %             end
% % % % % %         end
% % % % % %     end
% % % % % %     all_gain= AllData(j).Trial(1).clusters;
% % % % % %     endpoint_XY(:,3)=sqrt(endpoint_XY(:,1).^2 + endpoint_XY(:,2).^2);
% % % % % %     endpoint_XY=endpoint_XY-repmat(nanmean(endpoint_XY,1),size(endpoint_XY,1),1);
% % % % % %     for i=1:size(FR,1)-1
% % % % % %         for ii=i+1:size(FR,1)
% % % % % %             ch1=Chs{i};
% % % % % %             ch2=Chs{ii};
% % % % % %             if ~isempty(ch1) & ~isempty(ch2)
% % % % % %                 overlap= length(intersect(ch1,ch2))/length(unique([ch1,ch2]));
% % % % % %             else
% % % % % %                 overlap=0;
% % % % % %             end
% % % % % %             if overlap==0 & ((all_gain(i)<=2 & all_gain(ii)<=2) | (all_gain(i)<=2 & all_gain(ii)<=2))
% % % % % %                 c_all=c_all+1;
% % % % % %                 cell1= squeeze(FR(i,:,:));
% % % % % %                 cell2= squeeze(FR(ii,:,:));
% % % % % % 
% % % % % %                 tix= find(ax>-.15 & ax<0);
% % % % % %                 %                 sim_spikes=nansum((cell1(:,tix).*cell2(:,tix)),2);
% % % % % %                 sim_spikes=nanmean((cell1(:,tix).*cell2(:,tix)),2)-nanmean(nanmean((cell1(:,tix).*cell2(:,tix)),2));
% % % % % % 
% % % % % %                 [r_real(c_all),p_real(c_all)]= corr(sim_spikes,endpoint_XY(:,1));
% % % % % % 
% % % % % %                 [~,ix] = sort(sim_spikes);
% % % % % %                 low_ix=ix(1:round(length(ix)*prob));
% % % % % %                 endpoint_high2lowSI(1,c_all,1)=nanmean(endpoint_XY(low_ix,1));
% % % % % %                 endpoint_high2lowSI(2,c_all,1)=nanmean(endpoint_XY(low_ix,2));
% % % % % %                 endpoint_high2lowSI(3,c_all,1)=nanmean(endpoint_XY(low_ix,3));
% % % % % %                 high_ix=ix(end-round(length(ix)*prob):end);
% % % % % %                 endpoint_high2lowSI(1,c_all,2)=nanmean(endpoint_XY(high_ix,1));
% % % % % %                 endpoint_high2lowSI(2,c_all,2)=nanmean(endpoint_XY(high_ix,2));
% % % % % %                 endpoint_high2lowSI(3,c_all,2)=nanmean(endpoint_XY(high_ix,3));
% % % % % % 
% % % % % % 
% % % % % %                 [~,endix]=sort(endpoint_XY(:,1));
% % % % % %                 Pr_S1= cell1(endix(1:round(length(endix)*prob)),:);
% % % % % %                 Pr_S2= cell2(endix(1:round(length(endix)*prob)),:);
% % % % % %                 [raw, shift_predict, pred, surprise, std, ~, ~] = my_JPSTH(Pr_S1', Pr_S2', nbin,1);
% % % % % %                 jPSTH_hypo(c_all,:,:)=(raw-pred)./size(Pr_S1,1);%./std;%
% % % % % %                 all_fr_hypo(c_all,:)=nanmean([nanmean(Pr_S1,1); nanmean(Pr_S2,1)]);
% % % % % %                 maxV_hypo(c_all)=mean(maxV(endix(1:round(length(endix)*prob))));
% % % % % % 
% % % % % %                 Pr_S1= cell1(endix(end-round(length(endix)*prob):end),:);
% % % % % %                 Pr_S2= cell2(endix(end-round(length(endix)*prob):end),:);
% % % % % %                 [raw, shift_predict, pred, surprise, std, ~, ~] = my_JPSTH(Pr_S1', Pr_S2', nbin,1);
% % % % % %                 jPSTH_hyper(c_all,:,:)=(raw-pred)./size(Pr_S1,1);%std;%
% % % % % %                 maxV_hyper(c_all)=mean(maxV(endix(end-round(length(endix)*prob):end)));
% % % % % %                 all_fr_hyper(c_all,:)=nanmean([nanmean(Pr_S1,1); nanmean(Pr_S2,1)]);
% % % % % % 
% % % % % %                 gain(c_all,1)=all_gain(i);
% % % % % %                 gain(c_all,2)=all_gain(ii);
% % % % % %                 date_id(c_all)= j;
% % % % % %                 all_mice_id{c_all}=mouse_id;
% % % % % %             end
% % % % % %         end
% % % % % %     end
% % % % % %     clear FR FR2gain this_S nX nY nZ endpoint_XY Chs maxV
% % % % % % end
%%
load('E:\PAPER\Data2load\jPSTH_data\jPSTH2Endpoint.mat');
%%
pair2take=find(gain(:,1) & gain(:,2)==1);
%%
s=15;
gwin=gausswin(s); gwin=gwin./sum(gwin);
for i=1:size(all_fr_hyper,1)
    x= squeeze(all_fr_hyper(i,:));
    x=conv2(x,gwin,'same');
    smoothed_fr_hyper(i,:)=x;

    x=x-mean(x(find(ax<-.25)));
    x=conv2(x,gwin,'same');
    norm_fr_hyper(i,:)=x;


    x= squeeze(all_fr_hypo(i,:));
    x=conv2(x,gwin,'same');
    smoothed_fr_hypo(i,:)=x;

    x=x-mean(x(find(ax<-.25)));
    x=conv2(x,gwin,'same');
    norm_fr_hypo(i,:)=x;
end
%%
subplot(2,3,1)
pcolor(ax,ax,imgaussfilt(squeeze(nanmean(jPSTH_hyper,1)),3))
title('Hypermetric');
shading flat
colormap jet

subplot(2,3,2)
pcolor(ax,ax,imgaussfilt(squeeze(nanmean(jPSTH_hypo,1)),3))
title('Hypometric');
shading flat
colormap jet

subplot(2,3,3)
pcolor(ax,ax,imgaussfilt(squeeze(nanmean(jPSTH_hypo-jPSTH_hyper,1)),3))
title('Hypo- Hyper');
shading flat
colormap jet

%%
sm=21;
gwin=gausswin(sm);
gwin=gwin./sum(gwin);
for i=1:size(jPSTH_hyper,1)
    % diag_hypo(i,:)= smooth(diag(squeeze(jPSTH_hypo(i,:,:))),sm);
    % diag_hyper(i,:)= smooth(diag(squeeze(jPSTH_hyper(i,:,:))),sm);
    diag_hypo(i,:)= conv2(diag(squeeze(jPSTH_hypo(i,:,:))),gwin,'same');
    diag_hyper(i,:)= conv2(diag(squeeze(jPSTH_hyper(i,:,:))),gwin,'same');
end
ix= find(max(abs(diag_hyper),[],2)>1);
btime= -.2;
bix= find(ax>btime & ax<0);
bix= find(ax>-.15 & ax<.0);

subplot(2,3,4)
errorbar(ax,nanmean(diag_hypo,1),nanstd(diag_hypo,[],1)./sqrt(size(diag_hypo,1)),'Color','m')
hold on
errorbar(ax,nanmean(diag_hyper,1),nanstd(diag_hyper,[],1)./sqrt(size(diag_hyper,1)),'Color','g')
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
%
% histogram(r_real,'BinWidth',.02,'FaceColor','k','EdgeColor','none')
% hold on
% plot([nanmean(r_real) nanmean(r_real)],[0 25],'k','LineWidth',1)
[~,ix]=unique(maxV_hyper);
scatter(maxV_hyper(ix),maxV_hypo(ix),70,'ok','fill')
a=axis;
hold on
plot([min(a) max(a)],[min(a) max(a)],'-k')
hold off
[~,p]= ttest(maxV_hyper(ix),maxV_hypo(ix));
title(p)
set(gca,'TickDir','out'); box off
xlabel('max V; Hyper')
ylabel('max V; Hypo')