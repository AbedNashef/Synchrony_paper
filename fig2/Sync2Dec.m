% for Fig 2
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
prob=1/4;
%%
nax= [-0.5:1/1000:0.5];nax=nax(1:end-1);
t1= find(nax<-0.25); t2= find(nax>0 & nax<.25);
bhv_ax=[-tbf:1/freq:taf];

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
        %         endpoint_time=ReachS(i).out(ix,1);
        endpoint_time=ReachS(i).out(end,1);
        if isempty(endpoint_time)
            ix=find(ReachS(i).filt_kin(:,2)>move_thresh,1,'first');
            endpoint_time=ReachS(i).filt_kin(ix,1);
        end
        ix= find(ReachS(i).filt_kin(:,1)>=endpoint_time-tbf & ReachS(i).filt_kin(:,1)<=endpoint_time+taf);
        dt=median(diff(ReachS(i).filt_kin(ix,1)));
        tq=-tbf:1/freq:taf;
        t=ReachS(i).filt_kin(ix,1)-endpoint_time;
        xV=interp1(t,ReachS(i).filt_kin(ix,6),tq);
        yV=interp1(t,ReachS(i).filt_kin(ix,7),tq);
        V=interp1(t,ReachS(i).filt_kin(ix,5),tq);
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
                Chs{c_this}= cellData(cc).Channels;
                all_gain(c_this)=cellData(cc).gain;
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
                tix= find(ax>-.2 & ax<.0);
                Pr_S1= nanmean(cell1(:,tix),1);
                Pr_S2= nanmean(cell2(:,tix),1);
                Pr_S1S2= nanmean(cell1(:,tix).*cell2(:,tix),1);
                sumSI= nansum((cell1(:,tix).*cell2(:,tix))./repmat(Pr_S1.*Pr_S2,size(time_min,1),1),2);
                
                [~,ix] = sort(sumSI);
                low_ix=ix(1:round(length(ix)*prob));
                high_ix=ix(end-round(length(ix)*prob):end);
                dec_lowSI(1,c_all)=nanmean(xDec(low_ix));
                dec_lowSI(2,c_all)=nanmean(yDec(low_ix));
                dec_highSI(1,c_all)=nanmean(xDec(high_ix));
                dec_highSI(2,c_all)=nanmean(yDec(high_ix));
                
                [~,ix]= sort(abs(xDec));
                Pr_S1= nanmean(cell1(ix(1:round(length(ix)*prob)),:),1);
                Pr_S2= nanmean(cell2(ix(1:round(length(ix)*prob)),:),1);
                Pr_S1S2= nanmean(cell1(ix(1:round(length(ix)*prob)),:).*cell2(ix(1:round(length(ix)*prob)),:),1);
                all_sync_bottom(c_all,:)= Pr_S1S2./(Pr_S1.*Pr_S2);
                Pr_S1= nanmean(cell1(ix(end-round(length(ix)*prob):end),:),1);
                Pr_S2= nanmean(cell2(ix(end-round(length(ix)*prob):end),:),1);
                Pr_S1S2= nanmean(cell1(ix(end-round(length(ix)*prob):end),:).*cell2(ix(end-round(length(ix)*prob):end),:),1);
                all_sync_top(c_all,:)= Pr_S1S2./(Pr_S1.*Pr_S2);
                
                Pr_S1= nanmean(cell1,1);
                Pr_S2= nanmean(cell2,1);
                Pr_S1S2= nanmean(cell1.*cell2,1);
                all_sync(c_all,:)= Pr_S1S2./(Pr_S1.*Pr_S2);
                
                gain(c_all,1)=all_gain(i);
                gain(c_all,2)=all_gain(ii);
                date_id(c_all)= j;
                Dist4cells(c_all)= abs(cellData(i).depth-cellData(ii).depth);
                all_mice_id{c_all}=mouse_id;
            end
        end
    end
    clear FR FR2gain this_S nX nY nZ xVelocity yVelocity Speed xDec yDec
    clear time_min time_max minV maxV Chs
end
%%
orig_all_sync=all_sync;
orig_all_sync_bottom=all_sync_bottom;
orig_all_sync_top=all_sync_top;
%%
all_sync=orig_all_sync;
all_sync_bottom=orig_all_sync_bottom;
all_sync_top=orig_all_sync_top;
%%
alpha=.05;
sm=31;
axx=[-tbf:nbin/1000:taf]; axx=axx(1:end-1);
subplot(2,3,1)
ix= find(gain(:,1)<-.05 & gain(:,2)<-.05);
m=smooth(nanmean(all_sync_bottom(ix,:),1),sm); m=m(:)';
s=smooth(nanstd(all_sync_bottom(ix,:),1),sm)./sqrt(length(ix)); s=s(:)';
xx=[ax ax(end) ax(end) fliplr(ax) ax(1) ax(1)];
yy=[m-s m(end)-s(end) m(end)+s(end) fliplr(m+s) m(1)+s(1) m(1)-s(1)];
patch(xx,yy,[0.5 0.5 1],'EdgeColor','none');
hold on
m=smooth(nanmean(all_sync_top(ix,:),1),sm); m=m(:)';
s=smooth(nanstd(all_sync_top(ix,:),1),sm)./sqrt(length(ix)); s=s(:)';
xx=[ax ax(end) ax(end) fliplr(ax) ax(1) ax(1)];
yy=[m-s m(end)-s(end) m(end)+s(end) fliplr(m+s) m(1)+s(1) m(1)-s(1)];
patch(xx,yy,[1 0.5 0.5],'EdgeColor','none');
plot(ax,smooth(nanmean(all_sync_bottom(ix,:),1),sm),'Color','b','LineWidth',1)
plot(ax,smooth(nanmean(all_sync_top(ix,:),1),sm),'Color','r','LineWidth',1)
a=axis;
plot([0 0],[a(3) a(4)],'--k','LineWidth',2)
[~,p]= ttest(all_sync_top(ix,:),all_sync_bottom(ix,:));
m=smooth(nanmean(all_sync_top(ix,:),1),sm);
scatter(ax(find(p<alpha)),m(find(p<alpha)),40,'sr','fill')
hold off
set(gca,'TickDir','out'); box off
ylabel('SI'); xlabel('Time (s)')
title('Pausers')
set(gca,'TickDir','out'); box off

subplot(2,3,2)
ix= find(gain(:,1)>.05 & gain(:,2)>.05);
m=smooth(nanmean(all_sync_bottom(ix,:),1),sm); m=m(:)';
s=smooth(nanstd(all_sync_bottom(ix,:),1),sm)./sqrt(length(ix)); s=s(:)';
xx=[ax ax(end) ax(end) fliplr(ax) ax(1) ax(1)];
yy=[m-s m(end)-s(end) m(end)+s(end) fliplr(m+s) m(1)+s(1) m(1)-s(1)];
patch(xx,yy,[0.5 0.5 1],'EdgeColor','none');
hold on
m=smooth(nanmean(all_sync_top(ix,:),1),sm); m=m(:)';
s=smooth(nanstd(all_sync_top(ix,:),1),sm)./sqrt(length(ix)); s=s(:)';
xx=[ax ax(end) ax(end) fliplr(ax) ax(1) ax(1)];
yy=[m-s m(end)-s(end) m(end)+s(end) fliplr(m+s) m(1)+s(1) m(1)-s(1)];
patch(xx,yy,[1 0.5 0.5],'EdgeColor','none');
plot(ax,smooth(nanmean(all_sync_bottom(ix,:),1),sm),'Color','b','LineWidth',1)
plot(ax,smooth(nanmean(all_sync_top(ix,:),1),sm),'Color','r','LineWidth',1)
a=axis;
plot([0 0],[a(3) a(4)],'--k','LineWidth',2)
[~,p]= ttest(all_sync_top(ix,:),all_sync_bottom(ix,:));
m=smooth(nanmean(all_sync_top(ix,:),1),sm);
scatter(ax(find(p<alpha)),m(find(p<alpha)),40,'sr','fill')
hold off
set(gca,'TickDir','out'); box off
ylabel('SI'); xlabel('Time (s)')
title('Bursters')
set(gca,'TickDir','out'); box off

subplot(2,3,3)
ix= find((gain(:,1)>.05 & gain(:,2)<-.05) | (gain(:,1)<-.05 & gain(:,2)>.05));
m=smooth(nanmean(all_sync_bottom(ix,:),1),sm); m=m(:)';
s=smooth(nanstd(all_sync_bottom(ix,:),1),sm)./sqrt(length(ix)); s=s(:)';
xx=[ax ax(end) ax(end) fliplr(ax) ax(1) ax(1)];
yy=[m-s m(end)-s(end) m(end)+s(end) fliplr(m+s) m(1)+s(1) m(1)-s(1)];
patch(xx,yy,[0.5 0.5 1],'EdgeColor','none');
hold on
m=smooth(nanmean(all_sync_top(ix,:),1),sm); m=m(:)';
s=smooth(nanstd(all_sync_top(ix,:),1),sm)./sqrt(length(ix)); s=s(:)';
xx=[ax ax(end) ax(end) fliplr(ax) ax(1) ax(1)];
yy=[m-s m(end)-s(end) m(end)+s(end) fliplr(m+s) m(1)+s(1) m(1)-s(1)];
patch(xx,yy,[1 0.5 0.5],'EdgeColor','none');
plot(ax,smooth(nanmean(all_sync_bottom(ix,:),1),sm),'Color','b','LineWidth',1)
plot(ax,smooth(nanmean(all_sync_top(ix,:),1),sm),'Color','r','LineWidth',1)
a=axis;
plot([0 0],[a(3) a(4)],'--k','LineWidth',2)
[~,p]= ttest(all_sync_top(ix,:),all_sync_bottom(ix,:));
m=smooth(nanmean(all_sync_top(ix,:),1),sm);
scatter(ax(find(p<alpha)),m(find(p<alpha)),40,'sr','fill')
hold off
set(gca,'TickDir','out'); box off
ylabel('SI'); xlabel('Time (s)')
title('Mixed')
set(gca,'TickDir','out'); box off
%% Endpoint to SI
% figure
dec_highSI=abs(dec_highSI);
dec_lowSI=abs(dec_lowSI);
um=unique(all_mice_id);
clear Xdec_HIGH_SI Ydec_HIGH_SI XdecLOW_SI Ydec_Low_SI
for i=1:length(um)
    ix= find(strcmpi(all_mice_id,um{i}));
    udays=unique(date_id(ix));
    for ii=1:length(udays)
        ix=find(strcmpi(all_mice_id(:),um{i}) & date_id(:)==udays(ii) & gain(:,1)<-.05 & gain(:,2)<-.05);
        Xdec_HIGH_SI(i,ii)=mean(squeeze(dec_highSI(1,ix)));
        Ydec_HIGH_SI(i,ii)=mean(squeeze(dec_highSI(2,ix)));
        XdecLOW_SI(i,ii)=mean(squeeze(dec_lowSI(1,ix)));
        Ydec_Low_SI(i,ii)=mean(squeeze(dec_lowSI(2,ix)));
    end
end
ix0=find(Xdec_HIGH_SI==0);
Xdec_HIGH_SI(ix0)=nan;
Ydec_HIGH_SI(ix0)=nan;
XdecLOW_SI(ix0)=nan;
Ydec_Low_SI(ix0)=nan;

pause_x_highSI=Xdec_HIGH_SI;
pause_x_lowSI=XdecLOW_SI;
pause_y_highSI=Ydec_HIGH_SI;
pause_y_lowSI=Ydec_Low_SI;


subplot(2,3,4)
plot([.9 1.1],[Xdec_HIGH_SI(:) XdecLOW_SI(:)],'-','Color',[.5 .5 .5])
hold on
plot([1.9 2.1],[Ydec_HIGH_SI(:) Ydec_Low_SI(:)],'-','Color',[.5 .5 .5])
scatter([.9 ],nanmean([Xdec_HIGH_SI(:)],1),50,'sr','fill')
scatter([1.1],nanmean([XdecLOW_SI(:)],1),50,'sb','fill')
scatter([1.9 ],nanmean([Ydec_HIGH_SI(:)],1),50,'sr','fill')
scatter([2.1],nanmean([Ydec_Low_SI(:)],1),50,'sb','fill')
a=axis;
set(gca,'TickDir','out','Xtick',[1 2],'XTickLabel',{'X' 'Y'}); box off
ylabel('max Dec')
[~,p]= ttest(pause_x_highSI(:),pause_x_lowSI(:));
if p<.05
    scatter(1,a(4),50,'*k')
end
[~,p]= ttest(pause_y_highSI(:),pause_y_lowSI(:));
if p<.05
    scatter(1,a(4),50,'*k')
end
hold off

clear Xdec_HIGH_SI Ydec_HIGH_SI XdecLOW_SI Ydec_Low_SI
for i=1:length(um)
    ix= find(strcmpi(all_mice_id,um{i}));
    udays=unique(date_id(ix));
    for ii=1:length(udays)
        ix=find(strcmpi(all_mice_id(:),um{i}) & date_id(:)==udays(ii) & gain(:,1)>.05 & gain(:,2)>.05);
        Xdec_HIGH_SI(i,ii)=mean(squeeze(dec_highSI(1,ix)));
        Ydec_HIGH_SI(i,ii)=mean(squeeze(dec_highSI(2,ix)));
        XdecLOW_SI(i,ii)=mean(squeeze(dec_lowSI(1,ix)));
        Ydec_Low_SI(i,ii)=mean(squeeze(dec_lowSI(2,ix)));
    end
end
ix0=find(Xdec_HIGH_SI==0);
Xdec_HIGH_SI(ix0)=nan;
Ydec_HIGH_SI(ix0)=nan;
XdecLOW_SI(ix0)=nan;
Ydec_Low_SI(ix0)=nan;

Burst_x_highSI=Xdec_HIGH_SI;
Burst_x_lowSI=XdecLOW_SI;
Burst_y_highSI=Ydec_HIGH_SI;
Burst_y_lowSI=Ydec_Low_SI;


subplot(2,3,5)
plot([.9 1.1],[Xdec_HIGH_SI(:) XdecLOW_SI(:)],'-','Color',[.5 .5 .5])
hold on
plot([1.9 2.1],[Ydec_HIGH_SI(:) Ydec_Low_SI(:)],'-','Color',[.5 .5 .5])
scatter([.9 ],nanmean([Xdec_HIGH_SI(:)],1),50,'sr','fill')
scatter([1.1],nanmean([XdecLOW_SI(:)],1),50,'sb','fill')
scatter([1.9 ],nanmean([Ydec_HIGH_SI(:)],1),50,'sr','fill')
scatter([2.1],nanmean([Ydec_Low_SI(:)],1),50,'sb','fill')
a=axis;
set(gca,'TickDir','out','Xtick',[1 2],'XTickLabel',{'X' 'Y'}); box off
ylabel('max Dec')
[~,p]= ttest(Burst_x_lowSI(:),Burst_x_highSI(:));
if p<.05
    scatter(1,a(4),50,'*k')
end
[~,p]= ttest(Burst_y_highSI(:),Burst_y_lowSI(:));
if p<.05
    scatter(1,a(4),50,'*k')
end
hold off

clear Xdec_HIGH_SI Ydec_HIGH_SI XdecLOW_SI Ydec_Low_SI
for i=1:length(um)
    ix= find(strcmpi(all_mice_id,um{i}));
    udays=unique(date_id(ix));
    for ii=1:length(udays)
        ix=find(strcmpi(all_mice_id(:),um{i}) & date_id(:)==udays(ii) & ((gain(:,1)>.05 & gain(:,2)<=.05) | (gain(:,2)>.05 & gain(:,1)<=.05)));
        Xdec_HIGH_SI(i,ii)=mean(squeeze(dec_highSI(1,ix)));
        Ydec_HIGH_SI(i,ii)=mean(squeeze(dec_highSI(2,ix)));
        XdecLOW_SI(i,ii)=mean(squeeze(dec_lowSI(1,ix)));
        Ydec_Low_SI(i,ii)=mean(squeeze(dec_lowSI(2,ix)));
    end
end
ix0=find(Xdec_HIGH_SI==0);
Xdec_HIGH_SI(ix0)=nan;
Ydec_HIGH_SI(ix0)=nan;
XdecLOW_SI(ix0)=nan;
Ydec_Low_SI(ix0)=nan;

mixed_x_highSI=Xdec_HIGH_SI;
mixed_x_lowSI=XdecLOW_SI;
mixed_y_highSI=Ydec_HIGH_SI;
mixed_y_lowSI=Ydec_Low_SI;


subplot(2,3,6)
plot([.9 1.1],[Xdec_HIGH_SI(:) XdecLOW_SI(:)],'-','Color',[.5 .5 .5])
hold on
plot([1.9 2.1],[Ydec_HIGH_SI(:) Ydec_Low_SI(:)],'-','Color',[.5 .5 .5])
scatter([.9 ],nanmean([Xdec_HIGH_SI(:)],1),50,'sr','fill')
scatter([1.1],nanmean([XdecLOW_SI(:)],1),50,'sb','fill')
scatter([1.9 ],nanmean([Ydec_HIGH_SI(:)],1),50,'sr','fill')
scatter([2.1],nanmean([Ydec_Low_SI(:)],1),50,'sb','fill')
a=axis;
set(gca,'TickDir','out','Xtick',[1 2],'XTickLabel',{'X' 'Y'}); box off
ylabel('max Dec')
[~,p]= ttest(mixed_x_highSI(:),mixed_x_lowSI(:));
if p<.05
    scatter(1,a(4),50,'*k')
end
[~,p]= ttest(mixed_y_highSI(:),mixed_y_lowSI(:));
if p<.05
    scatter(1,a(4),50,'*k')
end
hold off
%%
ix= find(gain(:,1)<-.05 & gain(:,2)<-.05);
subplot(2,3,4)
for  i=1:size(dec_highSI,1)
    plot([i-.1 i+.1],[dec_lowSI(i,ix)' dec_highSI(i,ix)'],'k')
    hold on
    plot([i-.1 i+.1],nanmean([dec_lowSI(i,ix)' dec_highSI(i,ix)'],1),'r','LineWidth',2)
    
    [~,p] = ttest(dec_highSI(i,ix),dec_lowSI(i,ix));
    a=axis;
    if p<.05
        scatter(i,a(4),50,'*k')
    end
end
hold off
set(gca,'XTick',[1 2],'XTickLabel',{'X' 'Y'},'TickDir','out')
box off
ylabel('Dec')
title(['n=' num2str(length(ix))])

ix= find(gain(:,1)>.05 & gain(:,2)>.05);
subplot(2,3,5)
for  i=1:size(dec_highSI,1)
    plot([i-.1 i+.1],[dec_lowSI(i,ix)' dec_highSI(i,ix)'],'k')
    hold on
    plot([i-.1 i+.1],nanmean([dec_lowSI(i,ix)' dec_highSI(i,ix)'],1),'r','LineWidth',2)
    
    [~,p] = ttest(dec_highSI(i,ix),dec_lowSI(i,ix));
    a=axis;
    if p<.05
        scatter(i,a(4),50,'*k')
    end
end
hold off
set(gca,'XTick',[1 2],'XTickLabel',{'X' 'Y'},'TickDir','out')
box off
ylabel('Dec')
title(['n=' num2str(length(ix))])

ix= find((gain(:,1)<-.05 & gain(:,2)>.05) | (gain(:,1)>.05 & gain(:,2)<-.05));
subplot(2,3,6)
for  i=1:size(dec_highSI,1)
    plot([i-.1 i+.1],[dec_lowSI(i,ix)' dec_highSI(i,ix)'],'k')
    hold on
    plot([i-.1 i+.1],nanmean([dec_lowSI(i,ix)' dec_highSI(i,ix)'],1),'r','LineWidth',2)
    
    [~,p] = ttest(dec_highSI(i,ix),dec_lowSI(i,ix));
    a=axis;
    if p<.05
        scatter(i,a(4),50,'*k')
    end
end
hold off
set(gca,'XTick',[1 2],'XTickLabel',{'X' 'Y'},'TickDir','out')
box off
ylabel('Dec')
title(['n=' num2str(length(ix))])