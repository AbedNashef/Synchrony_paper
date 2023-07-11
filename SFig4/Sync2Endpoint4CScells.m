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
prob=1/4;
%%
nax= [-0.5:1/1000:0.5];nax=nax(1:end-1);
t1= find(nax<-0.25); t2= find(nax>0 & nax<.25);

% [ res shift_predict psth_pred surprise_mat std_mat NUM_TRIAL pval_mat] = my_JPSTH(cut1, cut2, BIN_SIZE,gauss_filt_std)
%%
for j=1:length(dates)
    load([udir dates(j).name]);
    c_trials=0;
    index=strfind(dates(j).name,'_');
    mouse_id=dates(j).name(index(end-1)+1:index(end)+1);
    if isfield(cellData,'CS_Bin1')
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
                    endpoint_XY(c_trials,:)=[ReachS(i).out(end,2) ReachS(i).out(end,3)];
                    
                    all_gain(c_this)=cellData(cc).gain;
                    CS_cell(c_this)= ~isempty(cellData(cc).CS_Bin1);
                                        Chs{c_this}= cellData(cc).Channels;

                end
            end
        end
        endpoint_XY(:,3)=sqrt(endpoint_XY(:,1).^2 + endpoint_XY(:,2).^2);
        endpoint_XY=endpoint_XY-repmat(nanmean(endpoint_XY,1),size(endpoint_XY,1),1);
        for i=1:size(FR,1)-1
            for ii=i+1:size(FR,1)
                ch1=Chs{i};
                ch2=Chs{ii};
                if ~isempty(ch1) & ~isempty(ch2)
                    overlap= length(intersect(ch1,ch2))/length(unique([ch1,ch2]));
                else
                    overlap=0;
                end
                if CS_cell(i) & CS_cell(ii) & overlap>=0
                    c_all=c_all+1;
                    cell1= squeeze(FR(i,:,:));
                    cell2= squeeze(FR(ii,:,:));
                    
                    tix= find(ax>-.5 & ax<0.75);
                    Pr_S1= nanmean(cell1(:,tix),1);
                    Pr_S2= nanmean(cell2(:,tix),1);
                    Pr_S1S2= nanmean(cell1(:,tix).*cell2(:,tix),1);
                    sumSI= nansum((cell1(:,tix).*cell2(:,tix))./repmat(Pr_S1.*Pr_S2,size(all_gain,1),1),2);
                    
                    [r_real(c_all),p_real(c_all)]= corr(sumSI,endpoint_XY(:,1));
                    
                    [~,ix] = sort(sumSI);
                    low_ix=ix(1:round(length(ix)*prob));
                    endpoint_high2lowSI(1,c_all,1)=nanmean(endpoint_XY(low_ix,1));
                    endpoint_high2lowSI(2,c_all,1)=nanmean(endpoint_XY(low_ix,2));
                    endpoint_high2lowSI(3,c_all,1)=nanmean(endpoint_XY(low_ix,3));
                    high_ix=ix(end-round(length(ix)*prob):end);
                    endpoint_high2lowSI(1,c_all,2)=nanmean(endpoint_XY(high_ix,1));
                    endpoint_high2lowSI(2,c_all,2)=nanmean(endpoint_XY(high_ix,2));
                    endpoint_high2lowSI(3,c_all,2)=nanmean(endpoint_XY(high_ix,3));
                                        
                    
                    sum_SI_all(c_all,:)=[nanmean(sumSI(low_ix)) nanmean(sumSI(high_ix))];
                    
                    Pr_S1= nanmean(cell1,1);
                    Pr_S2= nanmean(cell2,1);
                    [~,ix1]=sort(rand(1,size(cell1,1)));
                    [~,ix2]=sort(rand(1,size(cell2,1)));
                    Pr_S1S2= nanmean(cell1(ix1,:).*cell2(ix2,:),1);
                    
                    [r_rand(c_all),p_rand(c_all)]= corr(nansum((cell1(ix1,:).*cell2(ix2,:))./repmat(Pr_S1.*Pr_S2,size(endpoint_XY,1),1),2),endpoint_XY(:,1));
                    
                    avg_end(c_all,:)=nanmean(endpoint_XY,1);
                    
                    [~,endix]=sort(endpoint_XY(:,1));
                    Pr_S1= nanmean(cell1(endix(1:round(length(endix)*prob)),:),1);
                    Pr_S2= nanmean(cell2(endix(1:round(length(endix)*prob)),:),1);
                    Pr_S1S2= nanmean(cell1(endix(1:round(length(endix)*prob)),:).*cell2(endix(1:round(length(endix)*prob)),:),1);
                    all_sync_bottom(c_all,:)= Pr_S1S2./(Pr_S1.*Pr_S2);
                    Pr_S1= nanmean(cell1(endix(end-round(length(endix)*prob):end),:),1);
                    Pr_S2= nanmean(cell2(endix(end-round(length(endix)*prob):end),:),1);
                    Pr_S1S2= nanmean(cell1(endix(end-round(length(endix)*prob):end),:).*cell2(endix(end-round(length(endix)*prob):end),:),1);
                    all_sync_top(c_all,:)= Pr_S1S2./(Pr_S1.*Pr_S2);
                    avg_endpoint_change_all(1,c_all,1)=nanmean(endpoint_XY(endix(1:round(length(endix)*prob)),1));
                    avg_endpoint_change_all(2,c_all,1)=nanmean(endpoint_XY(endix(end-round(length(endix)*prob)),1));
                    avg_endpoint_change_all(1,c_all,2)=nanmean(endpoint_XY(endix(1:round(length(endix)*prob)),2));
                    avg_endpoint_change_all(2,c_all,2)=nanmean(endpoint_XY(endix(end-round(length(endix)*prob)),2));
                    
                    gain(c_all,1)=all_gain(i);
                    gain(c_all,2)=all_gain(ii);
                    date_id(c_all)= j;
                    Dist4cells(c_all)= abs(cellData(i).depth-cellData(ii).depth);
                    all_mice_id{c_all}=mouse_id;
                end
            end
        end
    end
    clear FR FR2gain this_S nX nY nZ endpoint_XY CS_cell Chs
end
%%
orig_all_sync_bottom=all_sync_bottom; % = undershoot
orig_all_sync_top=all_sync_top; % = overshoot
%%
all_sync_bottom=orig_all_sync_bottom;
all_sync_top=orig_all_sync_top;
%%
% sm=5;
% gwin=gausswin(sm);
% gwin=gwin./sum(gwin);
% for i=1:size(all_sync_bottom,1)
%     all_sync_bottom(i,:)= smooth(all_sync_bottom(i,:),sm);
%     all_sync_top(i,:)= smooth(all_sync_top(i,:),sm);
% end
%%
alpha=.05;
gain_thresh=.05;
sm=21;
axx=[-tbf:nbin/1000:taf]; axx=axx(1:end-1);
subplot(2,2,1)
ix= find(gain(:,1)<-gain_thresh & gain(:,2)<-gain_thresh);
m=smooth(nanmean(all_sync_bottom(ix,:),1),sm); m=m(:)';
s=smooth(nanstd(all_sync_bottom(ix,:),[],1),sm)./sqrt(length(ix)); s=s(:)';
xx=[ax ax(end) ax(end) fliplr(ax) ax(1) ax(1)];
yy=[m-s m(end)-s(end) m(end)+s(end) fliplr(m+s) m(1)+s(1) m(1)-s(1)];
patch(xx,yy,[0.5 0.5 1],'EdgeColor','none');
hold on
m=smooth(nanmean(all_sync_top(ix,:),1),sm); m=m(:)';
s=smooth(nanstd(all_sync_top(ix,:),[],1),sm)./sqrt(length(ix)); s=s(:)';
xx=[ax ax(end) ax(end) fliplr(ax) ax(1) ax(1)];
yy=[m-s m(end)-s(end) m(end)+s(end) fliplr(m+s) m(1)+s(1) m(1)-s(1)];
patch(xx,yy,[1 0.5 0.5],'EdgeColor','none');
plot(ax,smooth(nanmean(all_sync_bottom(ix,:),1),sm),'Color','b','LineWidth',1)
plot(ax,smooth(nanmean(all_sync_top(ix,:),1),sm),'Color','r','LineWidth',1)
a=axis;
plot([0 0],[a(3) a(4)],'--k','LineWidth',2)
[~,p]= ttest(all_sync_top(ix,:),all_sync_bottom(ix,:));
scatter(ax(find(p<alpha)),ones(1,length(find(p<alpha))).*a(4),20,'*k')
hold off
set(gca,'TickDir','out'); box off
ylabel('SI'); xlabel('Time (s)')
title('Pausers')
set(gca,'TickDir','out'); box off

subplot(2,2,2)
ix= find(gain(:,1)>gain_thresh & gain(:,2)>gain_thresh);
m=smooth(nanmean(all_sync_bottom(ix,:),1),sm); m=m(:)';
s=smooth(nanstd(all_sync_bottom(ix,:),[],1),sm)./sqrt(length(ix)); s=s(:)';
xx=[ax ax(end) ax(end) fliplr(ax) ax(1) ax(1)];
yy=[m-s m(end)-s(end) m(end)+s(end) fliplr(m+s) m(1)+s(1) m(1)-s(1)];
patch(xx,yy,[0.5 0.5 1],'EdgeColor','none');
hold on
m=smooth(nanmean(all_sync_top(ix,:),1),sm); m=m(:)';
s=smooth(nanstd(all_sync_top(ix,:),[],1),sm)./sqrt(length(ix)); s=s(:)';
xx=[ax ax(end) ax(end) fliplr(ax) ax(1) ax(1)];
yy=[m-s m(end)-s(end) m(end)+s(end) fliplr(m+s) m(1)+s(1) m(1)-s(1)];
patch(xx,yy,[1 0.5 0.5],'EdgeColor','none');
plot(ax,smooth(nanmean(all_sync_bottom(ix,:),1),sm),'Color','b','LineWidth',1)
plot(ax,smooth(nanmean(all_sync_top(ix,:),1),sm),'Color','r','LineWidth',1)
a=axis;
plot([0 0],[a(3) a(4)],'--k','LineWidth',2)
[~,p]= ttest(all_sync_top(ix,:),all_sync_bottom(ix,:));
scatter(ax(find(p<alpha)),ones(1,length(find(p<alpha))).*a(4),20,'*k')
hold off
set(gca,'TickDir','out'); box off
ylabel('SI'); xlabel('Time (s)')
title('Bursters')
set(gca,'TickDir','out'); box off

subplot(2,2,3)
ix= find((gain(:,1)>gain_thresh & gain(:,2)<-gain_thresh) | (gain(:,1)<-gain_thresh & gain(:,2)>gain_thresh));
m=smooth(nanmean(all_sync_bottom(ix,:),1),sm); m=m(:)';
s=smooth(nanstd(all_sync_bottom(ix,:),[],1),sm)./sqrt(length(ix)); s=s(:)';
xx=[ax ax(end) ax(end) fliplr(ax) ax(1) ax(1)];
yy=[m-s m(end)-s(end) m(end)+s(end) fliplr(m+s) m(1)+s(1) m(1)-s(1)];
patch(xx,yy,[0.5 0.5 1],'EdgeColor','none');
hold on
m=smooth(nanmean(all_sync_top(ix,:),1),sm); m=m(:)';
s=smooth(nanstd(all_sync_top(ix,:),[],1),sm)./sqrt(length(ix)); s=s(:)';
xx=[ax ax(end) ax(end) fliplr(ax) ax(1) ax(1)];
yy=[m-s m(end)-s(end) m(end)+s(end) fliplr(m+s) m(1)+s(1) m(1)-s(1)];
patch(xx,yy,[1 0.5 0.5],'EdgeColor','none');
plot(ax,smooth(nanmean(all_sync_bottom(ix,:),1),sm),'Color','b','LineWidth',1)
plot(ax,smooth(nanmean(all_sync_top(ix,:),1),sm),'Color','r','LineWidth',1)
a=axis;
plot([0 0],[a(3) a(4)],'--k','LineWidth',2)
[~,p]= ttest(all_sync_top(ix,:),all_sync_bottom(ix,:));
scatter(ax(find(p<alpha)),ones(1,length(find(p<alpha))).*a(4),20,'*k')
hold off
set(gca,'TickDir','out'); box off
ylabel('SI'); xlabel('Time (s)')
title('Mixed')
set(gca,'TickDir','out'); box off

m=[];
subplot(2,2,4)
ix= find(gain(:,1)<-gain_thresh & gain(:,2)<-gain_thresh);
histogram(r_real(ix),20,'Normalization','probability','FaceColor','b','EdgeColor','none')
m(1)=nanmean(r_real(ix));
hold on
ix= find(gain(:,1)>gain_thresh & gain(:,2)>gain_thresh);
histogram(r_real(ix),20,'Normalization','probability','FaceColor','r','EdgeColor','none')
m(2)=nanmean(r_real(ix));
ix= find((gain(:,1)>gain_thresh & gain(:,2)<-gain_thresh) | (gain(:,1)<-gain_thresh & gain(:,2)>gain_thresh));
histogram(r_real(ix),20,'Normalization','probability','FaceColor','m','EdgeColor','none')
m(3)=nanmean(r_real(ix));
a=axis;
plot([m(1) m(1)],[a(4) a(4)-(a(4)-a(3))/10],'b','LineWidth',2)
plot([m(2) m(2)],[a(4) a(4)-(a(4)-a(3))/10],'r','LineWidth',2)
plot([m(3) m(3)],[a(4) a(4)-(a(4)-a(3))/10],'m','LineWidth',2)
hold off
xlabel('Corr (SI-Endpoint)');
ylabel('%');
set(gca,'TickDir','out'); box off
%% Endpoint to SI
figure
um=unique(all_mice_id);
clear XENDPOINT_HIGH_SI YENDPOINT_HIGH_SI XENDPOINT_LOW_SI YENDPOINT_LOW_SI
for i=1:length(um)
    ix= find(strcmpi(all_mice_id,um{i}));
    udays=unique(date_id(ix));
    for ii=1:length(udays)
        ix=find(strcmpi(all_mice_id(:),um{i}) & date_id(:)==udays(ii) & gain(:,1)<-gain_thresh & gain(:,2)<-gain_thresh);
        XENDPOINT_HIGH_SI(i,ii)=mean(squeeze(endpoint_high2lowSI(1,ix,2)));
        YENDPOINT_HIGH_SI(i,ii)=mean(squeeze(endpoint_high2lowSI(2,ix,2)));
        XYENDPOINT_HIGH_SI(i,ii)=mean(squeeze(endpoint_high2lowSI(3,ix,2)));
        XENDPOINT_LOW_SI(i,ii)=mean(squeeze(endpoint_high2lowSI(1,ix,1)));
        YENDPOINT_LOW_SI(i,ii)=mean(squeeze(endpoint_high2lowSI(2,ix,1)));
        XYENDPOINT_LOW_SI(i,ii)=mean(squeeze(endpoint_high2lowSI(3,ix,1)));
    end
end
ix=find(gain(:,1)<-gain_thresh & gain(:,2)<-gain_thresh);
XENDPOINT_HIGH_SI=squeeze(endpoint_high2lowSI(1,ix,2));
XENDPOINT_LOW_SI=squeeze(endpoint_high2lowSI(1,ix,1));
YENDPOINT_HIGH_SI=squeeze(endpoint_high2lowSI(2,ix,2));
YENDPOINT_LOW_SI=squeeze(endpoint_high2lowSI(2,ix,1));
XYENDPOINT_HIGH_SI=squeeze(endpoint_high2lowSI(3,ix,2));
XYENDPOINT_LOW_SI=squeeze(endpoint_high2lowSI(3,ix,1));
% ix0=find(XENDPOINT_HIGH_SI==0);
% XENDPOINT_HIGH_SI(ix0)=nan;
% YENDPOINT_HIGH_SI(ix0)=nan;
% XENDPOINT_LOW_SI(ix0)=nan;
% YENDPOINT_LOW_SI(ix0)=nan;

pause_x_highSI=XENDPOINT_HIGH_SI;
pause_x_lowSI=XENDPOINT_LOW_SI;
pause_y_highSI=YENDPOINT_HIGH_SI;
pause_y_lowSO=YENDPOINT_LOW_SI;

subplot(2,4,1)
scatter(XENDPOINT_LOW_SI(:),YENDPOINT_LOW_SI(:),20,'o','MarkerFaceColor','b','MarkerFaceAlpha',.5,'MarkerEdgeColor','none')
hold on
scatter(XENDPOINT_HIGH_SI(:),YENDPOINT_HIGH_SI(:),20,'o','MarkerFaceColor','r','MarkerFaceAlpha',.5,'MarkerEdgeColor','none')
scatter(nanmean(nanmean(XENDPOINT_HIGH_SI,2)),nanmean(nanmean(YENDPOINT_HIGH_SI,2)),70,'o','MarkerFaceColor','r','MarkerFaceAlpha',1,'MarkerEdgeColor','none')
scatter(nanmean(nanmean(XENDPOINT_LOW_SI,2)),nanmean(nanmean(YENDPOINT_LOW_SI,2)),70,'o','MarkerFaceColor','b','MarkerFaceAlpha',1,'MarkerEdgeColor','none')
a=axis;
plot([0 0],a(3:4),'--k')
plot(a(1:2),[0 0],'--k')
hold off
set(gca,'TickDir','out'); box off
[~,p]= ttest([YENDPOINT_HIGH_SI(:); XENDPOINT_HIGH_SI(:)],[YENDPOINT_LOW_SI(:); XENDPOINT_LOW_SI(:)]);
title(p);
xlabel('\DeltaX (cm)');
ylabel('\DeltaY (cm)');
legend({'Low SI' 'High SI'})

subplot(2,4,5)
plot([.9 1.1],[XENDPOINT_HIGH_SI(:) XENDPOINT_LOW_SI(:)],'-','Color',[.5 .5 .5])
hold on
plot([1.9 2.1],[YENDPOINT_HIGH_SI(:) YENDPOINT_LOW_SI(:)],'-','Color',[.5 .5 .5])
plot([2.9 3.1],[XYENDPOINT_HIGH_SI(:) XYENDPOINT_LOW_SI(:)],'-','Color',[.5 .5 .5])
scatter([.9 1.1],nanmean([XENDPOINT_HIGH_SI(:) XENDPOINT_LOW_SI(:)],1),50,'sk','fill')
scatter([1.9 2.1],nanmean([YENDPOINT_HIGH_SI(:) YENDPOINT_LOW_SI(:)],1),50,'sk','fill')
scatter([2.9 3.1],nanmean([XYENDPOINT_HIGH_SI(:) XYENDPOINT_LOW_SI(:)],1),50,'sk','fill')
a=axis;
plot(a(1:2),[0 0],'--k')
hold off
set(gca,'TickDir','out','Xtick',[1 2],'XTickLabel',{'X' 'Y'}); box off
ylabel('\DeltaEndpoint')

clear XENDPOINT_HIGH_SI YENDPOINT_HIGH_SI XENDPOINT_LOW_SI YENDPOINT_LOW_SI
for i=1:length(um)
    ix= find(strcmpi(all_mice_id,um{i}));
    udays=unique(date_id(ix));
    for ii=1:length(udays)
        ix=find(strcmpi(all_mice_id(:),um{i}) & date_id(:)==udays(ii) & gain(:,1)>gain_thresh & gain(:,2)>gain_thresh);
        XENDPOINT_HIGH_SI(i,ii)=mean(squeeze(endpoint_high2lowSI(1,ix,2)));
        YENDPOINT_HIGH_SI(i,ii)=mean(squeeze(endpoint_high2lowSI(2,ix,2)));
        XYENDPOINT_HIGH_SI(i,ii)=mean(squeeze(endpoint_high2lowSI(3,ix,2)));
        XENDPOINT_LOW_SI(i,ii)=mean(squeeze(endpoint_high2lowSI(1,ix,1)));
        YENDPOINT_LOW_SI(i,ii)=mean(squeeze(endpoint_high2lowSI(2,ix,1)));
        XYENDPOINT_LOW_SI(i,ii)=mean(squeeze(endpoint_high2lowSI(3,ix,1)));
    end
end
ix=find(gain(:,1)>gain_thresh & gain(:,2)>gain_thresh);
XENDPOINT_HIGH_SI=squeeze(endpoint_high2lowSI(1,ix,2));
XENDPOINT_LOW_SI=squeeze(endpoint_high2lowSI(1,ix,1));
YENDPOINT_HIGH_SI=squeeze(endpoint_high2lowSI(2,ix,2));
YENDPOINT_LOW_SI=squeeze(endpoint_high2lowSI(2,ix,1));
XYENDPOINT_HIGH_SI=squeeze(endpoint_high2lowSI(3,ix,2));
XYENDPOINT_LOW_SI=squeeze(endpoint_high2lowSI(3,ix,1));
% 
% ix0=find(XENDPOINT_HIGH_SI==0);
% XENDPOINT_HIGH_SI(ix0)=nan;
% YENDPOINT_HIGH_SI(ix0)=nan;
% XENDPOINT_LOW_SI(ix0)=nan;
% YENDPOINT_LOW_SI(ix0)=nan;

burst_x_highSI=XENDPOINT_HIGH_SI;
burst_x_lowSI=XENDPOINT_LOW_SI;
burst_y_highSI=YENDPOINT_HIGH_SI;
burst_y_lowSO=YENDPOINT_LOW_SI;

subplot(2,4,2)
scatter(XENDPOINT_LOW_SI(:),YENDPOINT_LOW_SI(:),20,'o','MarkerFaceColor','b','MarkerFaceAlpha',.5,'MarkerEdgeColor','none')
hold on
scatter(XENDPOINT_HIGH_SI(:),YENDPOINT_HIGH_SI(:),20,'o','MarkerFaceColor','r','MarkerFaceAlpha',.5,'MarkerEdgeColor','none')
scatter(nanmean(nanmean(XENDPOINT_HIGH_SI,2)),nanmean(nanmean(YENDPOINT_HIGH_SI,2)),70,'o','MarkerFaceColor','r','MarkerFaceAlpha',1,'MarkerEdgeColor','none')
scatter(nanmean(nanmean(XENDPOINT_LOW_SI,2)),nanmean(nanmean(YENDPOINT_LOW_SI,2)),70,'o','MarkerFaceColor','b','MarkerFaceAlpha',1,'MarkerEdgeColor','none')
a=axis;
plot([0 0],a(3:4),'--k')
plot(a(1:2),[0 0],'--k')
hold off
set(gca,'TickDir','out'); box off
[~,p]= ttest([YENDPOINT_HIGH_SI(:); XENDPOINT_HIGH_SI(:)],[YENDPOINT_LOW_SI(:); XENDPOINT_LOW_SI(:)]);
title(p);
xlabel('\DeltaX (cm)');
ylabel('\DeltaY (cm)');

subplot(2,4,6)
plot([.9 1.1],[XENDPOINT_HIGH_SI(:) XENDPOINT_LOW_SI(:)],'-','Color',[.5 .5 .5])
hold on
plot([1.9 2.1],[YENDPOINT_HIGH_SI(:) YENDPOINT_LOW_SI(:)],'-','Color',[.5 .5 .5])
plot([2.9 3.1],[XYENDPOINT_HIGH_SI(:) XYENDPOINT_LOW_SI(:)],'-','Color',[.5 .5 .5])
scatter([.9 1.1],nanmean([XENDPOINT_HIGH_SI(:) XENDPOINT_LOW_SI(:)],1),50,'sk','fill')
scatter([1.9 2.1],nanmean([YENDPOINT_HIGH_SI(:) YENDPOINT_LOW_SI(:)],1),50,'sk','fill')
scatter([2.9 3.1],nanmean([XYENDPOINT_HIGH_SI(:) XYENDPOINT_LOW_SI(:)],1),50,'sk','fill')
a=axis;
plot(a(1:2),[0 0],'--k')
hold off
set(gca,'TickDir','out','Xtick',[1 2],'XTickLabel',{'X' 'Y'}); box off
ylabel('\DeltaEndpoint')


clear XENDPOINT_HIGH_SI YENDPOINT_HIGH_SI XENDPOINT_LOW_SI YENDPOINT_LOW_SI
for i=1:length(um)
    ix= find(strcmpi(all_mice_id,um{i}));
    udays=unique(date_id(ix));
    for ii=1:length(udays)
        ix=find(strcmpi(all_mice_id(:),um{i}) & date_id(:)==udays(ii) & ((gain(:,1)>gain_thresh & gain(:,2)<-gain_thresh) | (gain(:,2)>gain_thresh & gain(:,1)<-gain_thresh)));
        XENDPOINT_HIGH_SI(i,ii)=mean(squeeze(endpoint_high2lowSI(1,ix,2)));
        YENDPOINT_HIGH_SI(i,ii)=mean(squeeze(endpoint_high2lowSI(2,ix,2)));
        XYENDPOINT_HIGH_SI(i,ii)=mean(squeeze(endpoint_high2lowSI(3,ix,2)));
        XENDPOINT_LOW_SI(i,ii)=mean(squeeze(endpoint_high2lowSI(1,ix,1)));
        YENDPOINT_LOW_SI(i,ii)=mean(squeeze(endpoint_high2lowSI(2,ix,1)));
        XYENDPOINT_LOW_SI(i,ii)=mean(squeeze(endpoint_high2lowSI(3,ix,1)));
    end
end
ix=find(((gain(:,1)>gain_thresh & gain(:,2)<-gain_thresh) | (gain(:,2)>gain_thresh & gain(:,1)<-gain_thresh)));
XENDPOINT_HIGH_SI=squeeze(endpoint_high2lowSI(1,ix,2));
XENDPOINT_LOW_SI=squeeze(endpoint_high2lowSI(1,ix,1));
YENDPOINT_HIGH_SI=squeeze(endpoint_high2lowSI(2,ix,2));
YENDPOINT_LOW_SI=squeeze(endpoint_high2lowSI(2,ix,1));
XYENDPOINT_HIGH_SI=squeeze(endpoint_high2lowSI(3,ix,2));
XYENDPOINT_LOW_SI=squeeze(endpoint_high2lowSI(3,ix,1));
% 
% ix0=find(XENDPOINT_HIGH_SI==0);
% XENDPOINT_HIGH_SI(ix0)=nan;
% YENDPOINT_HIGH_SI(ix0)=nan;
% XENDPOINT_LOW_SI(ix0)=nan;
% YENDPOINT_LOW_SI(ix0)=nan;

burst_x_highSI=XENDPOINT_HIGH_SI;
burst_x_lowSI=XENDPOINT_LOW_SI;
burst_y_highSI=YENDPOINT_HIGH_SI;
burst_y_lowSO=YENDPOINT_LOW_SI;

subplot(2,4,3)
scatter(XENDPOINT_LOW_SI(:),YENDPOINT_LOW_SI(:),20,'o','MarkerFaceColor','b','MarkerFaceAlpha',.5,'MarkerEdgeColor','none')
hold on
scatter(XENDPOINT_HIGH_SI(:),YENDPOINT_HIGH_SI(:),20,'o','MarkerFaceColor','r','MarkerFaceAlpha',.5,'MarkerEdgeColor','none')
scatter(nanmean(nanmean(XENDPOINT_HIGH_SI,2)),nanmean(nanmean(YENDPOINT_HIGH_SI,2)),70,'o','MarkerFaceColor','r','MarkerFaceAlpha',1,'MarkerEdgeColor','none')
scatter(nanmean(nanmean(XENDPOINT_LOW_SI,2)),nanmean(nanmean(YENDPOINT_LOW_SI,2)),70,'o','MarkerFaceColor','b','MarkerFaceAlpha',1,'MarkerEdgeColor','none')
a=axis;
plot([0 0],a(3:4),'--k')
plot(a(1:2),[0 0],'--k')
hold off
set(gca,'TickDir','out'); box off
[~,p]= ttest([YENDPOINT_HIGH_SI(:); XENDPOINT_HIGH_SI(:)],[YENDPOINT_LOW_SI(:); XENDPOINT_LOW_SI(:)]);
title(p);
xlabel('\DeltaX (cm)');
ylabel('\DeltaY (cm)');

subplot(2,4,7)
plot([.9 1.1],[XENDPOINT_HIGH_SI(:) XENDPOINT_LOW_SI(:)],'-','Color',[.5 .5 .5])
hold on
plot([1.9 2.1],[YENDPOINT_HIGH_SI(:) YENDPOINT_LOW_SI(:)],'-','Color',[.5 .5 .5])
plot([2.9 3.1],[XYENDPOINT_HIGH_SI(:) XYENDPOINT_LOW_SI(:)],'-','Color',[.5 .5 .5])
scatter([.9 1.1],nanmean([XENDPOINT_HIGH_SI(:) XENDPOINT_LOW_SI(:)],1),50,'sk','fill')
scatter([1.9 2.1],nanmean([YENDPOINT_HIGH_SI(:) YENDPOINT_LOW_SI(:)],1),50,'sk','fill')
scatter([2.9 3.1],nanmean([XYENDPOINT_HIGH_SI(:) XYENDPOINT_LOW_SI(:)],1),50,'sk','fill')
a=axis;
plot(a(1:2),[0 0],'--k')
hold off
set(gca,'TickDir','out','Xtick',[1 2],'XTickLabel',{'X' 'Y'}); box off
ylabel('\DeltaEndpoint')

%%
subplot(2,4,4)
ix= find(gain(:,1)<-gain_thresh & gain(:,2)<-gain_thresh);
plot([.9 1.1],sum_SI_all(ix,:),'Color',[.7 .7 1],'LineWidth',.5)
hold on
plot([.9 1.1],nanmean(sum_SI_all(ix,:),1),'Color','b','LineWidth',2)

ix= find(gain(:,1)>gain_thresh & gain(:,2)>gain_thresh);
plot([1.9 2.1],sum_SI_all(ix,:),'Color',[1 .7 .7],'LineWidth',.5)
plot([1.9 2.1],nanmean(sum_SI_all(ix,:),1),'Color','r','LineWidth',2)

ix= find((gain(:,1)>gain_thresh & gain(:,2)<-gain_thresh) | (gain(:,1)<-gain_thresh & gain(:,2)>gain_thresh));
plot([2.9 3.1],sum_SI_all(ix,:),'Color',[1 .7 1],'LineWidth',.5)
plot([2.9 3.1],nanmean(sum_SI_all(ix,:),1),'Color','m','LineWidth',2)
hold off
set(gca,'TickDir','out','XTick',1:3,'XTickLabel',{'Pausers' 'Bursters' 'Mixed'}); box off;
ylabel('\SigmaSI')

subplot 248
undershoot= squeeze(avg_endpoint_change_all(1,:,:));
overshoot= squeeze(avg_endpoint_change_all(2,:,:));
scatter(undershoot(:,1),undershoot(:,2),'.b')
hold on
scatter(overshoot(:,1),overshoot(:,2),'.r')
scatter(nanmean(undershoot(:,1)),nanmean(undershoot(:,2)),50,'ob','fill')
scatter(nanmean(overshoot(:,1)),nanmean(overshoot(:,2)),50,'or','fill')
hold off
xlabel('\DeltaX')
ylabel('\DeltaY')
set(gca,'TickDir','out'); box off
legend({'Undershooth' 'Overshoot'})