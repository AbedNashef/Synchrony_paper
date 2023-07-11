%
clear
udir= '../new_data/';
dates = dir([udir 'D*']);
bin2use=1;
addpath ../helper_functions/
tbf=0.5; taf=.75;
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
move_thresh=0.2;
%%
nax= [-0.5:1/1000:0.5];nax=nax(1:end-1);
t1= find(nax<-0.25); t2= find(nax>0 & nax<.25);

% [ res shift_predict psth_pred surprise_mat std_mat NUM_TRIAL pval_mat] = my_JPSTH(cut1, cut2, BIN_SIZE,gauss_filt_std)
%%
for j=1:length(dates)
    load([udir dates(j).name]);
    c_trials=0;
    for i=1:length(ReachS)
        times = ReachS(i).filt_kin(:,1);
        ix= find(ReachS(i).out(:,2)>move_thresh,1,'first');
        endpoint_time=ReachS(i).out(ix,1);
        %         endpoint_time=ReachS(i).out(end,1);
        if isempty(endpoint_time)
            ix=find(ReachS(i).filt_kin(:,2)>move_thresh,1,'first');
            endpoint_time=ReachS(i).filt_kin(ix,1);
        end
        times= ReachS(i).filt_kin(:,1);
        P1= repmat(ReachS(i).filt_kin(1,2:4),size(ReachS(i).filt_kin,1),1);
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
                
                index= find(this_trc(:,1)>=endpoint_time-0.5 & this_trc(:,1)<=endpoint_time+0.5);
                FR2gain(c_this,c_trials,:)= this_trc(index,2);
                Chs{c_this}= cellData(cc).Channels;
                
            end
        end
    end
    
    for i=1:size(FR,1)-1
        for ii=i+1:size(FR,1)
            ch1=Chs{i};
            ch2=Chs{ii};
            if ~isempty(ch1) & ~isempty(ch2)
                overlap= length(intersect(ch1,ch2))/length(unique([ch1,ch2]));
            else
                overlap=0;
            end
            if overlap==0
                c_all=c_all+1;
                cell1= squeeze(FR(i,:,:));
                cell2= squeeze(FR(ii,:,:));
                
                Pr_S1= nanmean(cell1,1);
                Pr_S2= nanmean(cell2,1);
                Pr_S1S2= nanmean(cell1.*cell2,1);
                
                all_sync(c_all,:)= Pr_S1S2;
                
                
                Pr_S1= nanmean(cell1,1);
                Pr_S2= nanmean(cell2,1);
                [~,ix1]=sort(rand(1,size(cell1,1)));
                [~,ix2]=sort(rand(1,size(cell2,1)));
                Pr_S1S2= nanmean(cell1(ix1,:).*cell2(ix2,:),1);
                
                all_sync_random(c_all,:)= Pr_S1S2;
                
                fr1= nanmean(nanmean(squeeze(FR2gain(i,:,t1)))); fr2= nanmean(nanmean(squeeze(FR2gain(i,:,t2))));
                gain(c_all,1)=(fr2-fr1)/(fr1+fr2);
                fr1= nanmean(nanmean(squeeze(FR2gain(ii,:,t1)))); fr2= nanmean(nanmean(squeeze(FR2gain(ii,:,t2))));
                gain(c_all,2)=(fr2-fr1)/(fr1+fr2);
                date_id(c_all)= j;
                Dist4cells(c_all)= abs(cellData(i).depth-cellData(ii).depth);
            end
        end
    end
    
    clear FR FR2gain this_S nX nY nZ
end
%%
%% keep original data
orig_all_sync=all_sync;
orig_all_sync_random=all_sync_random;

%% smooth data
sm=11;
for i=1:size(all_sync,1)
    all_sync_random(i,:)=smooth(all_sync_random(i,:),sm);
    all_sync(i,:)=smooth(all_sync(i,:),sm);
end
%%
sm=1;
axx=[-tbf:nbin/1000:taf]; axx=axx(1:end-1);
subplot(2,2,1)
ix= find(gain(:,1)<-.05 & gain(:,2)<-.05);
m=smooth(nanmean(all_sync(ix,:),1),sm); m=m(:)';
s=smooth(nanstd(all_sync(ix,:),1),sm)./sqrt(length(ix)); s=s(:)';
xx=[ax ax(end) ax(end) fliplr(ax) ax(1) ax(1)];
yy=[m-s m(end)-s(end) m(end)+s(end) fliplr(m+s) m(1)+s(1) m(1)-s(1)];
patch(xx,yy,[0.5 0.5 1],'EdgeColor','none');
hold on
m=smooth(nanmean(all_sync_random(ix,:),1),sm); m=m(:)';
s=smooth(nanstd(all_sync_random(ix,:),1),sm)./sqrt(length(ix)); s=s(:)';
xx=[ax ax(end) ax(end) fliplr(ax) ax(1) ax(1)];
yy=[m-s m(end)-s(end) m(end)+s(end) fliplr(m+s) m(1)+s(1) m(1)-s(1)];
patch(xx,yy,[0.5 0.5 0.5],'EdgeColor','none');
plot(ax,smooth(nanmean(all_sync(ix,:),1),sm),'Color','b','LineWidth',1)
plot(ax,smooth(nanmean(all_sync_random(ix,:),1),sm),'Color','k','LineWidth',1)
a=axis;
plot([0 0],[a(3) a(4)],'--k','LineWidth',2)
[h,p]= ttest(all_sync(ix,:),all_sync_random(ix,:));
scatter(ax,ones(1,length(ax)).*a(3),20,-log(p),'s')
hold off
set(gca,'TickDir','out'); box off
ylabel('SI'); xlabel('Time (s)')
title('Pausers')

subplot(2,2,2)
ix= find(gain(:,1)>.05 & gain(:,2)>.05);
m=smooth(nanmean(all_sync(ix,:),1),sm); m=m(:)';
s=smooth(nanstd(all_sync(ix,:),1),sm)./sqrt(length(ix)); s=s(:)';
xx=[ax ax(end) ax(end) fliplr(ax) ax(1) ax(1)];
yy=[m-s m(end)-s(end) m(end)+s(end) fliplr(m+s) m(1)+s(1) m(1)-s(1)];
patch(xx,yy,[1 0.5 0.5],'EdgeColor','none');
hold on
m=smooth(nanmean(all_sync_random(ix,:),1),sm); m=m(:)';
s=smooth(nanstd(all_sync_random(ix,:),1),sm)./sqrt(length(ix)); s=s(:)';
xx=[ax ax(end) ax(end) fliplr(ax) ax(1) ax(1)];
yy=[m-s m(end)-s(end) m(end)+s(end) fliplr(m+s) m(1)+s(1) m(1)-s(1)];
patch(xx,yy,[0.5 0.5 0.5],'EdgeColor','none');
plot(ax,smooth(nanmean(all_sync(ix,:),1),sm),'Color','r','LineWidth',1)
plot(ax,smooth(nanmean(all_sync_random(ix,:),1),sm),'Color','k','LineWidth',1)
a=axis;
plot([0 0],[a(3) a(4)],'--k','LineWidth',2)
[h,p]= ttest(all_sync(ix,:),all_sync_random(ix,:));
scatter(ax,ones(1,length(ax)).*a(3),20,-log(p),'s')
hold off
set(gca,'TickDir','out'); box off
ylabel('SI'); xlabel('Time (s)')
title('Bursters')

subplot(2,2,3)
ix= find((gain(:,1)>.05 & gain(:,2)<-.05) | (gain(:,1)<-.05 & gain(:,2)>.05));
m=smooth(nanmean(all_sync(ix,:),1),sm); m=m(:)';
s=smooth(nanstd(all_sync(ix,:),1),sm)./sqrt(length(ix)); s=s(:)';
xx=[ax ax(end) ax(end) fliplr(ax) ax(1) ax(1)];
yy=[m-s m(end)-s(end) m(end)+s(end) fliplr(m+s) m(1)+s(1) m(1)-s(1)];
patch(xx,yy,[0.5 1 0.5],'EdgeColor','none');
hold on
m=smooth(nanmean(all_sync_random(ix,:),1),sm); m=m(:)';
s=smooth(nanstd(all_sync_random(ix,:),1),sm)./sqrt(length(ix)); s=s(:)';
xx=[ax ax(end) ax(end) fliplr(ax) ax(1) ax(1)];
yy=[m-s m(end)-s(end) m(end)+s(end) fliplr(m+s) m(1)+s(1) m(1)-s(1)];
patch(xx,yy,[0.5 0.5 0.5],'EdgeColor','none');
plot(ax,smooth(nanmean(all_sync(ix,:),1),sm),'Color','g','LineWidth',1)
plot(ax,smooth(nanmean(all_sync_random(ix,:),1),sm),'Color','k','LineWidth',1)
a=axis;
plot([0 0],[a(3) a(4)],'--k','LineWidth',2)
[h,p]= ttest(all_sync(ix,:),all_sync_random(ix,:));
scatter(ax,ones(1,length(ax)).*a(3),20,-log(p),'s')
hold off
set(gca,'TickDir','out'); box off
ylabel('SI'); xlabel('Time (s)')
title('Bursters-Pausers')

mat1=[]; mat2=[];
for i=1:size(all_sync,1)
    mat1(i,:)=nanmean(reshape(all_sync(i,:),10,[]),1);
    mat2(i,:)=nanmean(reshape(all_sync_random(i,:),10,[]),1);
end
nax= reshape(ax,10,[]);
nax=nax(1,:);
dec_all_sync=mat1;
dec_all_sync_rand=mat2;
Same_prob=dec_all_sync-dec_all_sync_rand~=0;
inc_prob=dec_all_sync-dec_all_sync_rand>0;
dec_prob=dec_all_sync-dec_all_sync_rand<0;
subplot(4,4,[11 12])
ix= find(gain(:,1)<-.05 & gain(:,2)<-.05);
plot(nax,nanmean(inc_prob(ix,:),1)./nanmean(Same_prob(ix,:),1),'b');
ix= find(gain(:,1)>.05 & gain(:,2)>.05);
hold on
plot(nax,nanmean(inc_prob(ix,:),1)./nanmean(Same_prob(ix,:),1),'r');
ix= find((gain(:,1)>.05 & gain(:,2)<-.05) | (gain(:,1)<-.05 & gain(:,2)>.05));
plot(nax,nanmean(inc_prob(ix,:),1)./nanmean(Same_prob(ix,:),1),'m');
a=axis;
plot(a(1:2),[.5 .5],'--k');
plot([0 0],a(3:4),'--k');
hold off
set(gca,'TickDir','out'); box off
xlabel('Time (s)')
ylabel('%increase')

% subplot(4,4,12)
% ix= find(gain(:,1)<-.05 & gain(:,2)<-.05);
% plot(ax,smooth(nanmean(all_sync(ix,:)-all_sync_random(ix,:),1),11),'b');
% ix= find(gain(:,1)>.05 & gain(:,2)>.05);
% hold on
% plot(ax,smooth(nanmean(all_sync(ix,:)-all_sync_random(ix,:),1),11),'r');
% ix= find((gain(:,1)>.05 & gain(:,2)<-.05) | (gain(:,1)<-.05 & gain(:,2)>.05));
% plot(ax,smooth(nanmean(all_sync(ix,:)-all_sync_random(ix,:),1),11),'m');
% a=axis;
% plot(a(1:2),[0 0],'--k');
% plot([0 0],a(3:4),'--k');
% hold off
% set(gca,'TickDir','out'); box off
% xlabel('Time (s)')
% ylabel('\DeltaSI')

inc_prob=all_sync>all_sync_random;
[~,PP]= ttest(all_sync',all_sync_random');
MM=nanmean(all_sync,2)>nanmean(all_sync_random,2);
subplot(4,4,15)
ix= find(gain(:,1)<-.05 & gain(:,2)<-.05);
l=length(find(MM(ix) & PP(ix)'<.05));
bar(1,l/length(ix),'FaceColor','b','EdgeColor','k');
ix= find(gain(:,1)>.05 & gain(:,2)>.05);
l=length(find(MM(ix) & PP(ix)'<.05));
hold on
bar(2,l/length(ix),'FaceColor','r','EdgeColor','k');
ix= find((gain(:,1)>.05 & gain(:,2)<-.05) | (gain(:,1)<-.05 & gain(:,2)>.05));
l=length(find(MM(ix) & PP(ix)'<.05));
bar(3,l/length(ix),'FaceColor','m','EdgeColor','k');
a=axis;
plot(a(1:2),[.5 .5],'--k')
hold off
set(gca,'TickDir','out'); box off
set(gca,'XTick',1:3,'XTickLabel',{'Pausers' 'Bursters' 'Mixed'})
ylabel('%Sig. Higher SI');
title('Bursters-Pausers')

inc_prob=all_sync>all_sync_random;
ix=find(ax>0 & ax<.5);
[~,PP]= ttest(all_sync(:,ix)',all_sync_random(:,ix)');
MM=nanmean(all_sync(:,ix),2)>nanmean(all_sync_random(:,ix),2);
subplot(4,4,16)
ix= find(gain(:,1)<-.05 & gain(:,2)<-.05);
l=length(find(MM(ix) & PP(ix)'<.05));
bar(1,l/length(ix),'FaceColor','b','EdgeColor','k');
ix= find(gain(:,1)>.05 & gain(:,2)>.05);
l=length(find(MM(ix) & PP(ix)'<.05));
hold on
bar(2,l/length(ix),'FaceColor','r','EdgeColor','k');
ix= find((gain(:,1)>.05 & gain(:,2)<-.05) | (gain(:,1)<-.05 & gain(:,2)>.05));
l=length(find(MM(ix) & PP(ix)'<.05));
bar(3,l/length(ix),'FaceColor','m','EdgeColor','k');
a=axis;
plot(a(1:2),[.5 .5],'--k')
hold off
set(gca,'TickDir','out'); box off
set(gca,'XTick',1:3,'XTickLabel',{'Pausers' 'Bursters' 'Mixed'})
ylabel('%Sig. Higher SI');
title('Bursters-Pausers')