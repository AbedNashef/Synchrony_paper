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
move_thresh=0.3;
all_CS_times=[];
iter=100;
%%
nax= [-0.5:1/1000:0.5];nax=nax(1:end-1);
t1= find(nax<-0.25); t2= find(nax>0 & nax<.25);

% [ res shift_predict psth_pred surprise_mat std_mat NUM_TRIAL pval_mat] = my_JPSTH(cut1, cut2, BIN_SIZE,gauss_filt_std)
%%
for j=1:length(dates)
    load([udir dates(j).name]);
    c_trials=0;
    if isfield(cellData,'CS_Bin1')
        for i=1:length(ReachS)
            times = ReachS(i).filt_kin(:,1);
            ix= find(ReachS(i).out(:,2)>move_thresh,1,'first');
            mo_time=ReachS(i).out(ix,1);
            endpoint_time= ReachS(i).out(ix,1);
            %             endpoint_time=ReachS(i).out(end,1);
            if isempty(mo_time)
                ix=find(ReachS(i).filt_kin(:,2)>move_thresh,1,'first');
                mo_time=ReachS(i).filt_kin(ix,1);
            end
            times= ReachS(i).filt_kin(:,1);
            P1= repmat(ReachS(i).filt_kin(1,2:4),size(ReachS(i).filt_kin,1),1);
            tt=times(find(times>=mo_time-tbf & times<=mo_time+taf));
            
            vStimMode = isfield(ReachS(i),'stim');
            if vStimMode, vStim=ReachS(i).stim; else, vStim=0; end
            if isempty(vStim), vStim=0; end
            vExclude=ReachS(i).exclude;
            c_this=0;
            if ~vExclude & ~vStim
                c_trials=c_trials+1;
                for cc=1:length(cellData)
                    if ~isempty(cellData(cc).CS_Bin1)
                        c_this=c_this+1;
                        this_trc=cellData(cc).Bin1;
                        index= find(this_trc(:,1)>=mo_time-tbf & this_trc(:,1)<=mo_time+taf);
                        FR(c_this,c_trials,:)= this_trc(index,2);
                        
                        this_trc=cellData(cc).CS_Bin1;
                        index= find(this_trc(:,1)>=mo_time & this_trc(:,1)<=mo_time+0.3);
                        CS_on(c_this,c_trials)=length(find(cellData(cc).CS_Bin1(index,2)))>0;
                        nCS(c_this,c_trials)=length(find(cellData(cc).CS_Bin1(index,2)));
                        gain(c_this)=cellData(cc).gain;
                        ix=find(cellData(cc).CS_Bin1(index,2));
                        all_CS_times=[all_CS_times; cellData(cc).CS_Bin1(index(ix),1)-mo_time];
                        Chs{c_this}= cellData(cc).Channels;
                        
                    end
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
                    ix= find(CS_on(i,:)==1 | CS_on(ii,:)==1);
                    Pr_S1= nanmean(cell1(ix,:),1);
                    Pr_S2= nanmean(cell2(ix,:),1);
                    Pr_S1S2= nanmean(cell1(ix,:).*cell2(ix,:),1);
                    
                    ix_on=ix;
                    all_sync_CSon(c_all,:)= Pr_S1S2./(Pr_S1.*Pr_S2);
                    
                    ix= find(CS_on(i,:)==0 & CS_on(ii,:)==0);
                    Pr_S1= nanmean(cell1(ix,:),1);
                    Pr_S2= nanmean(cell2(ix,:),1);
                    Pr_S1S2= nanmean(cell1(ix,:).*cell2(ix,:),1);
                    clear tmp dtmp;
                    all_sync_CSoff(c_all,:)= Pr_S1S2./(Pr_S1.*Pr_S2);
                    ix_off=ix;
                    if length(ix_off)>length(ix_on)
                        for rep=1:iter
                            [~,ix] = sort(rand(1,length(ix_off)));
                            ix=ix_off(ix(1:length(ix_on)));
                            Pr_S1= nanmean(cell1(ix,:),1);
                            Pr_S2= nanmean(cell2(ix,:),1);
                            Pr_S1S2= nanmean(cell1(ix,:).*cell2(ix,:),1);
                            tmp(rep,:)=Pr_S1S2./(Pr_S1.*Pr_S2);
                            dtmp(rep,:)=all_sync_CSoff(c_all,:)-tmp(rep,:);
                        end
                        all_sync_CSoff_matched(c_all,:)= nanmean(tmp,1);
                        all_dsync_matched(c_all,:)= nanmean(dtmp,1);
                    else
                        all_sync_CSoff_matched(c_all,:)= all_sync_CSoff(c_all,:);
                    end
                    clear tmp dtmp;
                    if length(ix_on)>length(ix_off)
                        for rep=1:iter
                            [~,ix] = sort(rand(1,length(ix_on)));
                            ix=ix_on(ix(1:length(ix_off)));
                            Pr_S1= nanmean(cell1(ix,:),1);
                            Pr_S2= nanmean(cell2(ix,:),1);
                            Pr_S1S2= nanmean(cell1(ix,:).*cell2(ix,:),1);
                            tmp(rep,:)=Pr_S1S2./(Pr_S1.*Pr_S2);
                            dtmp(rep,:)=tmp(rep,:)-all_sync_CSon(c_all,:);
                        end
                        all_sync_CSon_matched(c_all,:)=nanmean(tmp,1);
                        all_dsync_matched(c_all,:)= nanmean(dtmp,1);
                    else
                        all_sync_CSon_matched(c_all,:)= all_sync_CSon(c_all,:);
                    end
                    clear tmp;
                    if length(ix_on)<3 | length(ix_off)<3
                        all_sync_CSon_matched(c_all,:)=nan(1,length(ax));
                        all_sync_CSoff_matched(c_all,:)=nan(1,length(ax));
                        all_dsync_matched(c_all,:)=nan(1,length(ax));
                    end
                    
                    ix= find(CS_on(i,:)==1 & CS_on(ii,:)==1);
                    if ~isempty(ix)
                        Pr_S1= nanmean(cell1(ix,:),1);
                        Pr_S2= nanmean(cell2(ix,:),1);
                        Pr_S1S2= nanmean(cell1(ix,:).*cell2(ix,:),1);
                        ix_on=ix;
                        all_sync_ShareCS(c_all,:)= Pr_S1S2./(Pr_S1.*Pr_S2);
                    else
                        all_sync_ShareCS(c_all,:)= nan(1,length(ax));
                    end
                    % trial matching 
                    ix_both=ix;
                    if length(ix_off)>length(ix_both)
                        for rep=1:iter
                            [~,ix] = sort(rand(1,length(ix_off)));
                            ix=ix_off(ix(1:length(ix_both)));
                            Pr_S1= nanmean(cell1(ix,:),1);
                            Pr_S2= nanmean(cell2(ix,:),1);
                            Pr_S1S2= nanmean(cell1(ix,:).*cell2(ix,:),1);
                            tmp(rep,:)=Pr_S1S2./(Pr_S1.*Pr_S2);
                            dtmp(rep,:)=all_sync_CSoff(c_all,:)-tmp(rep,:);
                        end
                        all_sync_CSoff_matched4shared(c_all,:)= nanmean(tmp,1);
                        all_dsync_matched4shared(c_all,:)= nanmean(dtmp,1);
                    else
                        all_sync_CSoff_matched4shared(c_all,:)= all_sync_CSoff(c_all,:);
                    end
                    clear tmp dtmp;
                    if length(ix_both)>length(ix_off)
                        for rep=1:iter
                            [~,ix] = sort(rand(1,length(ix_both)));
                            ix=ix_both(ix(1:length(ix_off)));
                            Pr_S1= nanmean(cell1(ix,:),1);
                            Pr_S2= nanmean(cell2(ix,:),1);
                            Pr_S1S2= nanmean(cell1(ix,:).*cell2(ix,:),1);
                            tmp(rep,:)=Pr_S1S2./(Pr_S1.*Pr_S2);
                            dtmp(rep,:)=tmp(rep,:)-all_sync_CSon(c_all,:);
                        end
                        all_sync_CSon_matched4shared(c_all,:)=nanmean(tmp,1);
                        all_dsync_matched4shared(c_all,:)= nanmean(dtmp,1);
                    else
                        all_sync_CSon_matched4shared(c_all,:)= all_sync_CSon(c_all,:);
                    end
                    clear tmp;
                    if length(ix_both)<3 | length(ix_off)<3
                        all_sync_CSon_matched4shared(c_all,:)=nan(1,length(ax));
                        all_sync_CSoff_matched4shared(c_all,:)=nan(1,length(ax));
                        all_dsync_matched4shared(c_all,:)=nan(1,length(ax));
                    end
                    %
                    all_gain(c_all,:)=[gain(i) gain(ii)];
                    date_id(c_all)= j;
                    Dist4cells(c_all)= abs(cellData(i).depth-cellData(ii).depth);
                end
            end
        end
        clear FR FR2gain this_S nX nY nZ gain CS_on nCS Chs
    end
end
%% keep original data
orig_all_sync_CSon=all_sync_CSon;
orig_all_sync_CSon_matched=all_sync_CSon_matched;
orig_all_sync_CSoff=all_sync_CSoff;
orig_all_sync_CSoff_matched=all_sync_CSoff_matched;
orig_all_sync_ShareCS=all_sync_ShareCS;
orig_all_dsync_matched=all_dsync_matched;
orig_all_sync_CSoff_matched4shared=all_sync_CSoff_matched4shared;
orig_all_sync_CSon_matched4shared=all_sync_CSon_matched4shared;
orig_all_dsync_matched4shared=all_dsync_matched4shared;
%%
all_sync_CSon=orig_all_sync_CSon;
all_sync_CSon_matched=orig_all_sync_CSon_matched;
all_sync_CSoff=orig_all_sync_CSoff;
all_sync_CSoff_matched=orig_all_sync_CSoff_matched;
all_sync_ShareCS=orig_all_sync_ShareCS;
all_dsync_matched=orig_all_dsync_matched;
all_sync_CSoff_matched4shared=orig_all_sync_CSoff_matched4shared;
all_sync_CSon_matched4shared=orig_all_sync_CSon_matched4shared;
all_dsync_matched4shared=orig_all_dsync_matched4shared;
%% smooth data
sm=11;
for i=1:size(all_sync_CSon,1)
    all_sync_CSon(i,:)=smooth(all_sync_CSon(i,:),sm);
    all_sync_CSon_matched(i,:)=smooth(all_sync_CSon_matched(i,:),sm);
    all_sync_CSoff(i,:)=smooth(all_sync_CSoff(i,:),sm);
    all_sync_CSoff_matched(i,:)=smooth(all_sync_CSoff_matched(i,:),sm);
    all_sync_ShareCS(i,:)=smooth(all_sync_ShareCS(i,:),sm);
    all_dsync_matched(i,:)=smooth(all_dsync_matched(i,:),sm);
        all_sync_CSon_matched4shared(i,:)=smooth(all_sync_CSon_matched4shared(i,:),sm);
    all_sync_CSoff_matched4shared(i,:)=smooth(all_sync_CSoff_matched4shared(i,:),sm);
    all_dsync_matched4shared(i,:)=smooth(all_dsync_matched4shared(i,:),sm);
end
all_sync_CSon(find(all_sync_CSon<0))=nan;% smoothing problem
all_sync_CSon_matched(find(all_sync_CSon_matched<0))=nan;% smoothing problem
all_sync_CSoff(find(all_sync_CSoff<0))=nan;% smoothing problem
all_sync_CSoff_matched(find(all_sync_CSoff_matched<0))=nan;% smoothing problem
all_sync_ShareCS(find(all_sync_ShareCS<0))=nan;% smoothing problem
all_sync_CSon_matched4shared(find(all_sync_CSon_matched4shared<0))=nan;% smoothing problem
all_sync_CSoff_matched4shared(find(all_sync_CSoff_matched4shared<0))=nan;% smoothing problem
all_dsync_matched4shared(find(all_dsync_matched4shared<0))=nan;% smoothing problem
%%
sm=1;
axx=[-tbf:nbin/1000:taf]; axx=axx(1:end-1);
subplot(2,3,1)
ix= find(all_gain(:,1)<-.05 & all_gain(:,2)<-.05);
m=smooth(nanmean(all_sync_CSon(ix,:),1),sm); m=m(:)';
s=smooth(nanstd(all_sync_CSon(ix,:),[],1),sm)./sqrt(length(ix)); s=s(:)';
xx=[ax ax(end) ax(end) fliplr(ax) ax(1) ax(1)];
yy=[m-s m(end)-s(end) m(end)+s(end) fliplr(m+s) m(1)+s(1) m(1)-s(1)];
patch(xx,yy,[0.5 0.5 1],'EdgeColor','none');
hold on
m=smooth(nanmean(all_sync_CSoff(ix,:),1),sm); m=m(:)';
s=smooth(nanstd(all_sync_CSoff(ix,:),[],1),sm)./sqrt(length(ix)); s=s(:)';
xx=[ax ax(end) ax(end) fliplr(ax) ax(1) ax(1)];
yy=[m-s m(end)-s(end) m(end)+s(end) fliplr(m+s) m(1)+s(1) m(1)-s(1)];
patch(xx,yy,[0.5 0.5 0.5],'EdgeColor','none');
m=smooth(nanmean(all_sync_ShareCS(ix,:),1),sm); m=m(:)';
s=smooth(nanstd(all_sync_ShareCS(ix,:),[],1),sm)./sqrt(length(ix)); s=s(:)';
xx=[ax ax(end) ax(end) fliplr(ax) ax(1) ax(1)];
yy=[m-s m(end)-s(end) m(end)+s(end) fliplr(m+s) m(1)+s(1) m(1)-s(1)];
patch(xx,yy,[0.5 1 0.5],'EdgeColor','none');
plot(ax,smooth(nanmean(all_sync_CSon(ix,:),1),sm),'Color','b','LineWidth',1)
plot(ax,smooth(nanmean(all_sync_ShareCS(ix,:),1),sm),'Color','g','LineWidth',1)
plot(ax,smooth(nanmean(all_sync_CSoff(ix,:),1),sm),'Color','k','LineWidth',1)
a=axis;
plot([0 0],[a(3) a(4)],'--k','LineWidth',2)
[~,p]= ttest(all_sync_CSoff(ix,:),all_sync_CSon(ix,:));
scatter(ax(find(p<.05)),smooth(nanmean(all_sync_CSoff(ix,find(p<.05)),1),sm),'*k')
yyaxis right
histogram(all_CS_times,10,'Normalization','probability')
hold off
set(gca,'TickDir','out'); box off
yyaxis left
ylabel('SI'); xlabel('Time (s)')
title('Pausers')

subplot(2,3,2)
ix= find(all_gain(:,1)>.05 & all_gain(:,2)>.05);
m=smooth(nanmean(all_sync_CSon(ix,:),1),sm); m=m(:)';
s=smooth(nanstd(all_sync_CSon(ix,:),[],1),sm)./sqrt(length(ix)); s=s(:)';
xx=[ax ax(end) ax(end) fliplr(ax) ax(1) ax(1)];
yy=[m-s m(end)-s(end) m(end)+s(end) fliplr(m+s) m(1)+s(1) m(1)-s(1)];
patch(xx,yy,[0.5 0.5 1],'EdgeColor','none');
hold on
m=smooth(nanmean(all_sync_CSoff(ix,:),1),sm); m=m(:)';
s=smooth(nanstd(all_sync_CSoff(ix,:),[],1),sm)./sqrt(length(ix)); s=s(:)';
xx=[ax ax(end) ax(end) fliplr(ax) ax(1) ax(1)];
yy=[m-s m(end)-s(end) m(end)+s(end) fliplr(m+s) m(1)+s(1) m(1)-s(1)];
patch(xx,yy,[0.5 0.5 0.5],'EdgeColor','none');
m=smooth(nanmean(all_sync_ShareCS(ix,:),1),sm); m=m(:)';
s=smooth(nanstd(all_sync_ShareCS(ix,:),[],1),sm)./sqrt(length(ix)); s=s(:)';
xx=[ax ax(end) ax(end) fliplr(ax) ax(1) ax(1)];
yy=[m-s m(end)-s(end) m(end)+s(end) fliplr(m+s) m(1)+s(1) m(1)-s(1)];
patch(xx,yy,[0.5 1 0.5],'EdgeColor','none');
plot(ax,smooth(nanmean(all_sync_CSon(ix,:),1),sm),'Color','b','LineWidth',1)
plot(ax,smooth(nanmean(all_sync_ShareCS(ix,:),1),sm),'Color','g','LineWidth',1)
plot(ax,smooth(nanmean(all_sync_CSoff(ix,:),1),sm),'Color','k','LineWidth',1)
a=axis;
plot([0 0],[a(3) a(4)],'--k','LineWidth',2)
[~,p]= ttest(all_sync_CSoff(ix,:),all_sync_CSon(ix,:));
scatter(ax(find(p<.05)),smooth(nanmean(all_sync_CSoff(ix,find(p<.05)),1),sm),'*k')
yyaxis right
histogram(all_CS_times,10,'Normalization','probability')
hold off
set(gca,'TickDir','out'); box off
ylabel('SI'); xlabel('Time (s)')
title('Bursters')

subplot(2,3,3)
ix= find((all_gain(:,1)>.05 & all_gain(:,2)<-.05) | (all_gain(:,1)<-.05 & all_gain(:,2)>.05));
m=smooth(nanmean(all_sync_CSon(ix,:),1),sm); m=m(:)';
s=smooth(nanstd(all_sync_CSon(ix,:),[],1),sm)./sqrt(length(ix)); s=s(:)';
xx=[ax ax(end) ax(end) fliplr(ax) ax(1) ax(1)];
yy=[m-s m(end)-s(end) m(end)+s(end) fliplr(m+s) m(1)+s(1) m(1)-s(1)];
patch(xx,yy,[0.5 0.5 1],'EdgeColor','none');
hold on
m=smooth(nanmean(all_sync_CSoff(ix,:),1),sm); m=m(:)';
s=smooth(nanstd(all_sync_CSoff(ix,:),[],1),sm)./sqrt(length(ix)); s=s(:)';
xx=[ax ax(end) ax(end) fliplr(ax) ax(1) ax(1)];
yy=[m-s m(end)-s(end) m(end)+s(end) fliplr(m+s) m(1)+s(1) m(1)-s(1)];
patch(xx,yy,[0.5 0.5 0.5],'EdgeColor','none');
m=smooth(nanmean(all_sync_ShareCS(ix,:),1),sm); m=m(:)';
s=smooth(nanstd(all_sync_ShareCS(ix,:),[],1),sm)./sqrt(length(ix)); s=s(:)';
xx=[ax ax(end) ax(end) fliplr(ax) ax(1) ax(1)];
yy=[m-s m(end)-s(end) m(end)+s(end) fliplr(m+s) m(1)+s(1) m(1)-s(1)];
patch(xx,yy,[0.5 1 0.5],'EdgeColor','none');
plot(ax,smooth(nanmean(all_sync_CSon(ix,:),1),sm),'Color','b','LineWidth',1)
plot(ax,smooth(nanmean(all_sync_ShareCS(ix,:),1),sm),'Color','g','LineWidth',1)
plot(ax,smooth(nanmean(all_sync_CSoff(ix,:),1),sm),'Color','k','LineWidth',1)
a=axis;
plot([0 0],[a(3) a(4)],'--k','LineWidth',2)
[~,p]= ttest(all_sync_CSoff(ix,:),all_sync_CSon(ix,:));
scatter(ax(find(p<.05)),smooth(nanmean(all_sync_CSoff(ix,find(p<.05)),1),sm),'*k')
yyaxis right
histogram(all_CS_times,10,'Normalization','probability')
hold off
set(gca,'TickDir','out'); box off
yyaxis left
ylabel('SI'); xlabel('Time (s)')
title('Bursters-Pausers')


subplot(2,3,4)
ix= 1:size(all_gain,1);
m=smooth(nanmean(all_sync_CSon(ix,:),1),sm); m=m(:)';
s=smooth(nanstd(all_sync_CSon(ix,:),[],1),sm)./sqrt(length(ix)); s=s(:)';
xx=[ax ax(end) ax(end) fliplr(ax) ax(1) ax(1)];
yy=[m-s m(end)-s(end) m(end)+s(end) fliplr(m+s) m(1)+s(1) m(1)-s(1)];
patch(xx,yy,[0.5 0.5 1],'EdgeColor','none');
hold on
m=smooth(nanmean(all_sync_CSoff(ix,:),1),sm); m=m(:)';
s=smooth(nanstd(all_sync_CSoff(ix,:),[],1),sm)./sqrt(length(ix)); s=s(:)';
xx=[ax ax(end) ax(end) fliplr(ax) ax(1) ax(1)];
yy=[m-s m(end)-s(end) m(end)+s(end) fliplr(m+s) m(1)+s(1) m(1)-s(1)];
patch(xx,yy,[0.5 0.5 0.5],'EdgeColor','none');
m=smooth(nanmean(all_sync_ShareCS(ix,:),1),sm); m=m(:)';
s=smooth(nanstd(all_sync_ShareCS(ix,:),[],1),sm)./sqrt(length(ix)); s=s(:)';
xx=[ax ax(end) ax(end) fliplr(ax) ax(1) ax(1)];
yy=[m-s m(end)-s(end) m(end)+s(end) fliplr(m+s) m(1)+s(1) m(1)-s(1)];
patch(xx,yy,[0.5 1 0.5],'EdgeColor','none');
plot(ax,smooth(nanmean(all_sync_CSon(ix,:),1),sm),'Color','b','LineWidth',1)
plot(ax,smooth(nanmean(all_sync_ShareCS(ix,:),1),sm),'Color','g','LineWidth',1)
plot(ax,smooth(nanmean(all_sync_CSoff(ix,:),1),sm),'Color','k','LineWidth',1)
a=axis;
plot([0 0],[a(3) a(4)],'--k','LineWidth',2)
[~,p]= ttest(all_sync_CSoff(ix,:),all_sync_CSon(ix,:));
scatter(ax(find(p<.05)),smooth(nanmean(all_sync_CSoff(ix,find(p<.05)),1),sm),'*k')
yyaxis right
histogram(all_CS_times,10,'Normalization','probability')
hold off
set(gca,'TickDir','out'); box off
yyaxis left
ylabel('SI'); xlabel('Time (s)')
title('All')

subplot(2,3,5)
tix= find(ax>0 & ax<.3);
avgs(:,1)= nanmean(all_sync_ShareCS(:,tix),2);
avgs(:,2)= nanmean(all_sync_CSon(:,tix),2);
avgs(:,3)= nanmean(all_sync_CSoff(:,tix),2);
bar([.9 1 1.1],nanmean(avgs,1),'FaceColor','w','EdgeColor','k')
hold on
plot([.9 1 1.1],[avgs(:,1) avgs(:,2) avgs(:,3)],'Color',[.6 .6 .6])
plot([.9 1 1.1],nanmean(avgs,1),'Color','k')
[~,p]= ttest(avgs(:,2),avgs(:,3));
if p<.05
    text(1,10,'*','FontName','Times')
end
ix= find(all_gain(:,1)<-.05 & all_gain(:,2)<-.05);
bar([1.9 2 2.1],nanmean(avgs(ix,:),1),'FaceColor',[.9 .9 1],'EdgeColor','k')
plot([1.9 2 2.1],avgs(ix,:),'Color','b')
plot([1.9 2 2.1],nanmean(avgs(ix,:),1),'Color','k')
[~,p]= ttest(avgs(ix,2),avgs(ix,3));
if p<.05
    text(2,10,'*','FontName','Times')
end
ix= find(all_gain(:,1)>.05 & all_gain(:,2)>.05);
bar([2.9 3 3.1],nanmean(avgs(ix,:),1),'FaceColor',[1 .9 .9],'EdgeColor','k')
plot([2.9 3 3.1],avgs(ix,:),'Color','r')
plot([2.9 3 3.1],nanmean(avgs(ix,:),1),'Color','k')
[~,p]= ttest(avgs(ix,2),avgs(ix,3));
if p<.05
    text(3,10,'*','FontName','Times')
end
ix= find((all_gain(:,1)>.05 & all_gain(:,2)<-.05) | (all_gain(:,1)<-.05 & all_gain(:,2)>.05));
bar([3.9 4 4.1],nanmean(avgs(ix,:),1),'FaceColor',[1 .9 1],'EdgeColor','k')
plot([3.9 4 4.1],avgs(ix,:),'Color','m')
plot([3.9 4 4.1],nanmean(avgs(ix,:),1),'Color','k')
[~,p]= ttest(avgs(ix,2),avgs(ix,3));
if p<.05
    text(4,10,'*','FontName','Times')
end
hold off
ylabel('Avg. SI (0-0.3 s)')
set(gca,'TickDir','out'); box off
set(gca,'XTick',[1:4],'XTickLabel',{'All' 'Pausers' 'Burster' 'Mixed'})

subplot(2,3,6)
ix= find((all_gain(:,1)<-.05 & all_gain(:,2)<-.05));
m=smooth(nanmean(all_sync_CSoff_matched(ix,:)-all_sync_CSon_matched(ix,:),1),sm); m=m(:)';
s=smooth(nanstd(all_sync_CSoff_matched(ix,:)-all_sync_CSon_matched(ix,:),[],1),sm)./sqrt(length(ix)); s=s(:)';
m=smooth(nanmean(all_dsync_matched(ix,:),1),sm); m=m(:)';
s=smooth(nanstd(all_dsync_matched(ix,:),[],1),sm)./sqrt(length(ix)); s=s(:)';
xx=[ax ax(end) ax(end) fliplr(ax) ax(1) ax(1)];
yy=[m-s m(end)-s(end) m(end)+s(end) fliplr(m+s) m(1)+s(1) m(1)-s(1)];
patch(xx,yy,[0.5 0.5 1],'EdgeColor','none');
hold on
plot(ax,m,'Color','b')
tmp=all_dsync_matched(ix,:);
[~,p]= ttest(tmp);
scatter(ax(find(p<.05)),m(find(p<.05)),'*b')
ix= find((all_gain(:,1)>.05 & all_gain(:,2)>.05));
m=smooth(nanmean(all_dsync_matched(ix,:),1),sm); m=m(:)';
s=smooth(nanstd(all_dsync_matched(ix,:),[],1),sm)./sqrt(length(ix)); s=s(:)';
xx=[ax ax(end) ax(end) fliplr(ax) ax(1) ax(1)];
yy=[m-s m(end)-s(end) m(end)+s(end) fliplr(m+s) m(1)+s(1) m(1)-s(1)];
patch(xx,yy,[1 0.5 .5],'EdgeColor','none');
tmp=all_dsync_matched(ix,:);
plot(ax,m,'Color','r')
[~,p]= ttest(tmp);
scatter(ax(find(p<.05)),m(find(p<.05)),'*r')
ix= find((all_gain(:,1)<-.05 & all_gain(:,2)>.05) | (all_gain(:,1)>.05 & all_gain(:,2)<-.05));
m=smooth(nanmean(all_dsync_matched(ix,:),1),sm); m=m(:)';
s=smooth(nanstd(all_dsync_matched(ix,:),[],1),sm)./sqrt(length(ix)); s=s(:)';
xx=[ax ax(end) ax(end) fliplr(ax) ax(1) ax(1)];
yy=[m-s m(end)-s(end) m(end)+s(end) fliplr(m+s) m(1)+s(1) m(1)-s(1)];
patch(xx,yy,[1 0.5 1],'EdgeColor','none');
plot(ax,m,'Color','m')
tmp=all_dsync_matched(ix,:);
[~,p]= ttest(tmp);
scatter(ax(find(p<.05)),m(find(p<.05)),'*m')
a=axis;
plot([0 0],[a(3) a(4)],'--k','LineWidth',2)
plot([a(1) a(2)],[0 0],'--k','LineWidth',2)
hold off
ylabel('\DeltaSI (off-on)'); xlabel('Time (s)')
title('Matched trials pairs')
set(gca,'TickDir','out'); box off
