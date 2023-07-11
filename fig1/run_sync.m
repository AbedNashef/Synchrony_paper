% script to calculate and plot the SI for the first figure.
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
bhv_ax=-tbf:1/freq:taf;
% [ res shift_predict psth_pred surprise_mat std_mat NUM_TRIAL pval_mat] = my_JPSTH(cut1, cut2, BIN_SIZE,gauss_filt_std)
%%
for j=1:length(dates)
    load([udir dates(j).name]);
    index=strfind(dates(j).name,'_');
    mouse_id=dates(j).name(index(end-1)+1:index(end)+1);
    c_trials=0;
    for i=1:length(ReachS)
        times = ReachS(i).filt_kin(:,1);
        ix= find(ReachS(i).out(:,2)>move_thresh,1,'first');
        endpoint_time=ReachS(i).out(ix,1); % find time of threshold crossing
%         endpoint_time=ReachS(i).out(end,1);
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
        thisx=interp1(t,ReachS(i).filt_kin(ix,2),tq);
        times= ReachS(i).filt_kin(:,1);
        
        vStimMode = isfield(ReachS(i),'stim');
        if vStimMode, vStim=ReachS(i).stim; else, vStim=0; end
        if isempty(vStim), vStim=0; end
        vExclude=ReachS(i).exclude;
        c_this=0;
        if ~vExclude & ~vStim
            c_trials=c_trials+1;
            all_x(c_trials,:)= thisx;
            for cc=1:length(cellData)
                c_this=c_this+1;
                this_trc=cellData(cc).Bin1;
                index= find(this_trc(:,1)>=endpoint_time-tbf & this_trc(:,1)<=endpoint_time+taf);
                FR(c_this,c_trials,:)= this_trc(index,2);
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
                
                % Synchrony index calcualtion
                Pr_S1= nanmean(cell1,1); %%%%%%%%%%%%%%%%%%%%
                Pr_S2= nanmean(cell2,1); %%%%%%%%%%%%%%%%%%%%
                Pr_S1S2= nanmean(cell1.*cell2,1); %%%%%%%%%%%
                
                all_sync(c_all,:)= Pr_S1S2./(Pr_S1.*Pr_S2);%%%%%%%%%%%%%%%%%%%%
                c=cov(cell1,cell2);
                all_cov(c_all,1)=c(2,1);
                
                Pr_S1= nanmean(cell1,1);
                Pr_S2= nanmean(cell2,1);
                [~,ix1]=sort(rand(1,size(cell1,1)));
                [~,ix2]=sort(rand(1,size(cell2,1)));
                Pr_S1S2= nanmean(cell1(ix1,:).*cell2(ix2,:),1);
                c=cov(cell1(ix1,:),cell2(ix2,:));
                all_cov(c_all,2)=c(2,1);
                all_sync_random(c_all,:)= Pr_S1S2./(Pr_S1.*Pr_S2);
                
                gain(c_all,1)=all_gain(i);
                gain(c_all,2)=all_gain(ii);
                date_id(c_all)= j;
                all_mouse_id{c_all}= mouse_id;
                X_all_days(c_all,:)= nanmean(all_x,1);
                Dist4cells(c_all)= abs(cellData(i).depth-cellData(ii).depth);
            end
        end
    end
    clear FR FR2gain this_S nX nY nZ all_gain all_x Chs
end
bax= [-tbf:1/freq:taf];
%% keep original data
orig_all_sync=all_sync;
orig_all_sync_random=all_sync_random;

%%
sm=21;
axx=[-tbf:nbin/1000:taf]; axx=axx(1:end-1);
subplot(2,2,1)
ix= find(gain(:,1)<-.05 & gain(:,2)<-.05);
m=smooth(nanmean(all_sync(ix,:),1),sm); m=m(:)';
s=smooth(nanstd(all_sync(ix,:),[],1),sm)./sqrt(length(ix)); s=s(:)';
xx=[ax ax(end) ax(end) fliplr(ax) ax(1) ax(1)];
yy=[m-s m(end)-s(end) m(end)+s(end) fliplr(m+s) m(1)+s(1) m(1)-s(1)];
patch(xx,yy,[0.5 0.5 1],'EdgeColor','none');
hold on
m=smooth(nanmean(all_sync_random(ix,:),1),sm); m=m(:)';
s=smooth(nanstd(all_sync_random(ix,:),[],1),sm)./sqrt(length(ix)); s=s(:)';
xx=[ax ax(end) ax(end) fliplr(ax) ax(1) ax(1)];
yy=[m-s m(end)-s(end) m(end)+s(end) fliplr(m+s) m(1)+s(1) m(1)-s(1)];
patch(xx,yy,[0.5 0.5 0.5],'EdgeColor','none');
plot(ax,smooth(nanmean(all_sync(ix,:),1),sm),'Color','b','LineWidth',1)
plot(ax,smooth(nanmean(all_sync_random(ix,:),1),sm),'Color','k','LineWidth',1)
a=axis;
plot([0 0],[a(3) a(4)],'--k','LineWidth',2)
hold off
set(gca,'TickDir','out'); box off
ylabel('SI'); xlabel('Time (s)')
title('Pausers')

subplot(2,2,2)
ix= find(gain(:,1)>.05 & gain(:,2)>.05);
m=smooth(nanmean(all_sync(ix,:),1),sm); m=m(:)';
s=smooth(nanstd(all_sync(ix,:),[],1),sm)./sqrt(length(ix)); s=s(:)';
xx=[ax ax(end) ax(end) fliplr(ax) ax(1) ax(1)];
yy=[m-s m(end)-s(end) m(end)+s(end) fliplr(m+s) m(1)+s(1) m(1)-s(1)];
patch(xx,yy,[1 0.5 0.5],'EdgeColor','none');
hold on
m=smooth(nanmean(all_sync_random(ix,:),1),sm); m=m(:)';
s=smooth(nanstd(all_sync_random(ix,:),[],1),sm)./sqrt(length(ix)); s=s(:)';
xx=[ax ax(end) ax(end) fliplr(ax) ax(1) ax(1)];
yy=[m-s m(end)-s(end) m(end)+s(end) fliplr(m+s) m(1)+s(1) m(1)-s(1)];
patch(xx,yy,[0.5 0.5 0.5],'EdgeColor','none');
plot(ax,smooth(nanmean(all_sync(ix,:),1),sm),'Color','r','LineWidth',1)
plot(ax,smooth(nanmean(all_sync_random(ix,:),1),sm),'Color','k','LineWidth',1)
a=axis;
plot([0 0],[a(3) a(4)],'--k','LineWidth',2)
hold off
set(gca,'TickDir','out'); box off
ylabel('SI'); xlabel('Time (s)')
title('Bursters')

subplot(2,2,3)
ix= find((gain(:,1)>.05 & gain(:,2)<-.05) | (gain(:,1)<-.05 & gain(:,2)>.05));
m=smooth(nanmean(all_sync(ix,:),1),sm); m=m(:)';
s=smooth(nanstd(all_sync(ix,:),[],1),sm)./sqrt(length(ix)); s=s(:)';
xx=[ax ax(end) ax(end) fliplr(ax) ax(1) ax(1)];
yy=[m-s m(end)-s(end) m(end)+s(end) fliplr(m+s) m(1)+s(1) m(1)-s(1)];
patch(xx,yy,[0.5 1 0.5],'EdgeColor','none');
hold on
m=smooth(nanmean(all_sync_random(ix,:),1),sm); m=m(:)';
s=smooth(nanstd(all_sync_random(ix,:),[],1),sm)./sqrt(length(ix)); s=s(:)';
xx=[ax ax(end) ax(end) fliplr(ax) ax(1) ax(1)];
yy=[m-s m(end)-s(end) m(end)+s(end) fliplr(m+s) m(1)+s(1) m(1)-s(1)];
patch(xx,yy,[0.5 0.5 0.5],'EdgeColor','none');
plot(ax,smooth(nanmean(all_sync(ix,:),1),sm),'Color','g','LineWidth',1)
plot(ax,smooth(nanmean(all_sync_random(ix,:),1),sm),'Color','k','LineWidth',1)
a=axis;
plot([0 0],[a(3) a(4)],'--k','LineWidth',2)
hold off
set(gca,'TickDir','out'); box off
ylabel('SI'); xlabel('Time (s)')
title('Bursters-Pausers')


subplot(2,2,4)
ix= find((gain(:,1)<-.05 & gain(:,2)<-.05));
m=smooth(nanmean(all_sync(ix,:),1),sm); m=m(:)';
s=smooth(nanstd(all_sync(ix,:),[],1),sm)./sqrt(length(ix)); s=s(:)';
xx=[ax ax(end) ax(end) fliplr(ax) ax(1) ax(1)];
yy=[m-s m(end)-s(end) m(end)+s(end) fliplr(m+s) m(1)+s(1) m(1)-s(1)];
patch(xx,yy,[0.5 0.5 1],'EdgeColor','none');
hold on
plot(ax,smooth(nanmean(all_sync(ix,:),1),sm),'Color','b','LineWidth',1)
ix= find((gain(:,1)>.05 & gain(:,2)>.05));
m=smooth(nanmean(all_sync(ix,:),1),sm); m=m(:)';
s=smooth(nanstd(all_sync(ix,:),[],1),sm)./sqrt(length(ix)); s=s(:)';
xx=[ax ax(end) ax(end) fliplr(ax) ax(1) ax(1)];
yy=[m-s m(end)-s(end) m(end)+s(end) fliplr(m+s) m(1)+s(1) m(1)-s(1)];
patch(xx,yy,[1 0.5 0.5],'EdgeColor','none');
plot(ax,smooth(nanmean(all_sync(ix,:),1),sm),'Color','r','LineWidth',1)

ix= find((gain(:,1)>.05 & gain(:,2)<-.05) | (gain(:,1)<-.05 & gain(:,2)>.05));
m=smooth(nanmean(all_sync(ix,:),1),sm); m=m(:)';
s=smooth(nanstd(all_sync(ix,:),[],1),sm)./sqrt(length(ix)); s=s(:)';
xx=[ax ax(end) ax(end) fliplr(ax) ax(1) ax(1)];
yy=[m-s m(end)-s(end) m(end)+s(end) fliplr(m+s) m(1)+s(1) m(1)-s(1)];
patch(xx,yy,[0.5 1 0.5],'EdgeColor','none');
plot(ax,smooth(nanmean(all_sync(ix,:),1),sm),'Color','g','LineWidth',1)
a=axis;
plot([0 0],[a(3) a(4)],'--k','LineWidth',2)
hold off
set(gca,'TickDir','out'); box off
ylabel('SI'); xlabel('Time (s)')
