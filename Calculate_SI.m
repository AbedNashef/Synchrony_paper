% By Abdulraheem Nashef July 11, 2023
% General function to calculate the synchrony index between PCs
clear
udir= '../new_data/';
dates = dir([udir 'D*']);
bin2use=1; % the synchrony will be calculated in 1 ms time bin
addpath ../helper_functions/
tbf=0.5; taf=.75; % in seconds
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
iterations=50; % for shuffled data
%% time axes
nax= [-0.5:1/1000:0.5];nax=nax(1:end-1);
t1= find(nax<-0.25); t2= find(nax>0 & nax<.25);
bhv_ax=-tbf:1/freq:taf;
%% SI calculation
for j=1:length(dates)
    load([udir dates(j).name]); % Load data
    index=strfind(dates(j).name,'_');
    mouse_id=dates(j).name(index(end-1)+1:index(end)+1);
    c_trials=0;
    for i=1:length(ReachS)
        times = ReachS(i).filt_kin(:,1);
        ix= find(ReachS(i).out(:,2)>move_thresh,1,'first');
        endpoint_time=ReachS(i).out(ix,1); % find time of threshold crossing
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
        thisx=interp1(t,ReachS(i).filt_kin(ix,2),tq);
        thisv=interp1(t,ReachS(i).filt_kin(ix,6),tq);
        times= ReachS(i).filt_kin(:,1);
        
        vStimMode = isfield(ReachS(i),'stim');
        if vStimMode, vStim=ReachS(i).stim; else, vStim=0; end
        if isempty(vStim), vStim=0; end
        vExclude=ReachS(i).exclude;
        c_this=0;
        if ~vExclude & ~vStim
            c_trials=c_trials+1;
            all_x(c_trials,:)= thisx;
            all_v(c_trials,:)= thisv;
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
                
                all_cov(c_all,2)=0;
                for iter=1:iterations
                    Pr_S1= nanmean(cell1,1);
                    Pr_S2= nanmean(cell2,1);
                    [~,ix1]=sort(rand(1,size(cell1,1)));
                    [~,ix2]=sort(rand(1,size(cell2,1)));
                    Pr_S1S2= nanmean(cell1(ix1,:).*cell2(ix2,:),1);
                    c=cov(cell1(ix1,:),cell2(ix2,:));
                    all_cov(c_all,2)=all_cov(c_all,2)+c(2,1);
                    tmp(iter,:)= Pr_S1S2./(Pr_S1.*Pr_S2);
                end
                all_sync_random(c_all,:)=nanmean(tmp,1);
                gain(c_all,1)=all_gain(i);
                gain(c_all,2)=all_gain(ii);
                date_id(c_all)= j;
                all_mouse_id{c_all}= mouse_id;
                X_all_days(c_all,:)= nanmean(all_x,1);
                V_all_days(c_all,:)= nanmean(all_v,1);
                Dist4cells(c_all)= abs(cellData(i).depth-cellData(ii).depth);
            end
        end
    end
    clear FR FR2gain this_S nX nY nZ all_gain all_x Chs all_v
end
bax= [-tbf:1/freq:taf];
%% keep original data
orig_all_sync=all_sync;
orig_all_sync_random=all_sync_random;

%% smooth data
sm=11;
all_sync=orig_all_sync;
all_sync_random=orig_all_sync_random;
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
%%
figure
d=20;
axx=[-tbf:nbin/1000:taf]; axx=axx(1:end-1);
subplot(2,2,1)
ix= find(gain(:,1)<-.05 & gain(:,2)<-.05);
scatter(axx,nanmean(all_sync(ix,:)),20,'MarkerFaceColor',[0.7 0.7 1],'MarkerEdgeColor','w');
hold on
plot(axx(1:d:end),nanmean(all_sync(ix,1:d:end),1),'b','LineWidth',2)
a=axis;
plot([0 0],[a(3) a(4)],'--k','LineWidth',2)
plot([a(1) a(2)],[1 1],'--k','LineWidth',2)
hold off
set(gca,'TickDir','out'); box off
ylabel('SI'); xlabel('Time (s)')
title('Pausers')

subplot(2,2,2)
ix= find(gain(:,1)>.05 & gain(:,2)>.05);
scatter(axx,nanmean(all_sync(ix,:)),20,'MarkerFaceColor',[1 0.7 .7],'MarkerEdgeColor','w');
hold on
plot(axx(1:d:end),nanmean(all_sync(ix,1:d:end),1),'r','LineWidth',2)
a=axis;
plot([0 0],[a(3) a(4)],'--k','LineWidth',2)
plot([a(1) a(2)],[1 1],'--k','LineWidth',2)
hold off
set(gca,'TickDir','out'); box off
ylabel('SI'); xlabel('Time (s)')
title('Bursters')

subplot(2,2,3)
ix= find((gain(:,1)>.05 & gain(:,2)<-.05) | (gain(:,1)<-.05 & gain(:,2)>.05));
scatter(axx,nanmean(all_sync(ix,:)),20,'MarkerFaceColor',[0.7 1 .7],'MarkerEdgeColor','w');
hold on
plot(axx(1:d:end),nanmean(all_sync(ix,1:d:end),1),'g','LineWidth',2)
a=axis;
plot([0 0],[a(3) a(4)],'--k','LineWidth',2)
plot([a(1) a(2)],[1 1],'--k','LineWidth',2)
hold off
set(gca,'TickDir','out'); box off
ylabel('SI'); xlabel('Time (s)')
title('Bursters-Pausers')

subplot(2,2,4)
ix= find((gain(:,1)<-.05 & gain(:,2)<-.05));
m=smooth(nanmean(all_sync(ix,:),1),sm); m=m(:)';
s=smooth(nanstd(all_sync(ix,:),1),sm)./sqrt(length(ix)); s=s(:)';
xx=[ax ax(end) ax(end) fliplr(ax) ax(1) ax(1)];
yy=[m-s m(end)-s(end) m(end)+s(end) fliplr(m+s) m(1)+s(1) m(1)-s(1)];
patch(xx,yy,[0.5 0.5 1],'EdgeColor','none');
hold on
plot(ax,smooth(nanmean(all_sync(ix,:),1),sm),'Color','b','LineWidth',1)
ix= find((gain(:,1)>.05 & gain(:,2)>.05));
m=smooth(nanmean(all_sync(ix,:),1),sm); m=m(:)';
s=smooth(nanstd(all_sync(ix,:),1),sm)./sqrt(length(ix)); s=s(:)';
xx=[ax ax(end) ax(end) fliplr(ax) ax(1) ax(1)];
yy=[m-s m(end)-s(end) m(end)+s(end) fliplr(m+s) m(1)+s(1) m(1)-s(1)];
patch(xx,yy,[1 0.5 0.5],'EdgeColor','none');
plot(ax,smooth(nanmean(all_sync(ix,:),1),sm),'Color','r','LineWidth',1)

ix= find((gain(:,1)>.05 & gain(:,2)<-.05) | (gain(:,1)<-.05 & gain(:,2)>.05));
m=smooth(nanmean(all_sync(ix,:),1),sm); m=m(:)';
s=smooth(nanstd(all_sync(ix,:),1),sm)./sqrt(length(ix)); s=s(:)';
xx=[ax ax(end) ax(end) fliplr(ax) ax(1) ax(1)];
yy=[m-s m(end)-s(end) m(end)+s(end) fliplr(m+s) m(1)+s(1) m(1)-s(1)];
patch(xx,yy,[0.5 1 0.5],'EdgeColor','none');
plot(ax,smooth(nanmean(all_sync(ix,:),1),sm),'Color','g','LineWidth',1)
a=axis;
plot([0 0],[a(3) a(4)],'--k','LineWidth',2)
hold off
set(gca,'TickDir','out'); box off
ylabel('SI'); xlabel('Time (s)')
%%
%%
sm=11;
gwin=gausswin(sm)';
gwin=gwin./sum(gwin);
%
% for i=1:size(all_sync,1)
%     all_sync(i,:)= smooth(all_sync(i,:),sm);%gwin,'same');
% end
%%
sm=5;
gwin=gausswin(sm)';
gwin=gwin./sum(gwin);
n=5;
pauser_dist= Dist4cells(find(gain(:,1)<-.05 & gain(:,2)<-.05));
burster_dist= Dist4cells(find(gain(:,1)>.05 & gain(:,2)>.05));
mix_dist= Dist4cells(find((gain(:,1)>.05 & gain(:,2)<-.05) | (gain(:,1)<-.05 & gain(:,2)>.05)));
pausers=all_sync(find(gain(:,1)<-.05 & gain(:,2)<-.05),:);
bursters=all_sync(find(gain(:,1)>.05 & gain(:,2)>.05),:);
mixers=all_sync(find((gain(:,1)>.05 & gain(:,2)<-.05) | (gain(:,1)<-.05 & gain(:,2)>.05)),:);

[~,dix_pause]= sort(pauser_dist);
[~,dix_burst]= sort(burster_dist);
[~,dix_mix]= sort(mix_dist);

sax= round(linspace(1,length(dix_pause),n));
tmp=linspace(0,1,n);
col=[flipud(tmp(:)) ones(n,1).*.1 tmp(:)];
subplot 221
for i=1:n-1
    plot(ax,smooth(nanmean(pausers(dix_pause(sax(i):sax(i+1)),:),1),sm),'color',col(i,:),'LineWidth',2)
    hold on
    DD(i)= nanmean(pauser_dist(dix_pause(sax(i):sax(i+1))));
    avg_max(i)= nanmean(max(pausers(dix_pause(sax(i):sax(i+1)),:),[],2));
    smax(i)= nanmean(max(pausers(dix_pause(sax(i):sax(i+1)),:),[],2))./sqrt(sax(i+1)-sax(i)+1);
end
xlim([-tbf+0.01 taf-0.01])
a=axis;
plot([0 0],a(3:4),'--k');
plot(a(1:2),[1 1],'--k');
hold off
xlabel('Time (s)'); ylabel('SI');
title('Pausers')
set(gca,'TickDir','out'); box off

sax= round(linspace(1,length(dix_pause),11));
for i=1:10
    DD(i)= nanmean(pauser_dist(dix_pause(sax(i):sax(i+1))));
    avg_max(i)= nanmean(nansum(pausers(dix_pause(sax(i):sax(i+1)),find(ax>-.1 & ax<.2)),2));
    smax(i)= nanmean(nansum(pausers(dix_pause(sax(i):sax(i+1)),find(ax>-.1 & ax<.2)),2))./sqrt(sax(i+1)-sax(i)+1);
end

subplot 247
errorbar(DD,avg_max,smax,'-s','MarkerFaceColor','b','MarkerEdgeColor','b','Color','b','CapSize',0,'LineWidth',2);

[RR,PP]= corr(pauser_dist(:),pausers);
subplot(6,4,16)
plot(ax,RR,'b','LineWidth',2)
hold on
scatter(ax(find(PP<.05/sm)),RR(find(PP<.05/sm)),30,'o','MarkerFaceColor','b','MarkerEdgeColor','none')
plot([ax(1) ax(end)],[0 0],'--k');
a=axis;
plot([0 0],a(3:4),'--k');
hold off
set(gca,'TickDir','out'); box off
xlabel('Time (s)'); ylabel('Corr. Coeff;');
xlim([-tbf+0.01 taf-0.01])
title('Dist. to SI corr');

sax= round(linspace(1,length(dix_burst),n));
subplot 222
for i=1:n-1
    plot(ax,smooth(nanmean(bursters(dix_burst(sax(i):sax(i+1)),:),1),sm),'color',col(i,:),'LineWidth',2)
    hold on
end
xlim([-tbf+0.01 taf-0.01])
a=axis;
plot([0 0],a(3:4),'--k');
plot(a(1:2),[1 1],'--k');
hold off
xlabel('Time (s)'); ylabel('SI');
title('Bursters')
set(gca,'TickDir','out'); box off

sax= round(linspace(1,length(dix_burst),11));
for i=1:10
    DD(i)= nanmean(burster_dist(dix_burst(sax(i):sax(i+1))));
    avg_max(i)= nanmean(nansum(bursters(dix_burst(sax(i):sax(i+1)),find(ax>-.1 & ax<.2)),2));
    smax(i)= nanmean(nansum(bursters(dix_burst(sax(i):sax(i+1)),find(ax>-.1 & ax<.2)),2))./sqrt(sax(i+1)-sax(i)+1);
end

subplot 247
hold on
errorbar(DD,avg_max,smax,'-s','MarkerFaceColor','r','MarkerEdgeColor','r','Color','r','CapSize',0,'LineWidth',2);

[RR,PP]= corr(burster_dist(:),bursters);
subplot(6,4,20)
plot(ax,RR,'r','LineWidth',2)
hold on
scatter(ax(find(PP<.05/sm)),RR(find(PP<.05/sm)),30,'o','MarkerFaceColor','r','MarkerEdgeColor','none')
plot([ax(1) ax(end)],[0 0],'--k');
a=axis;
plot([0 0],a(3:4),'--k');
hold off
set(gca,'TickDir','out'); box off
xlabel('Time (s)');
xlim([-tbf+0.01 taf-0.01])

sax= round(linspace(1,length(dix_mix),n));
subplot 223
for i=1:n-1
    plot(ax,smooth(nanmean(mixers(dix_mix(sax(i):sax(i+1)),:),1),sm),'color',col(i,:),'LineWidth',2)
    hold on
end
xlim([-tbf+0.01 taf-0.01])
a=axis;
plot([0 0],a(3:4),'--k');
plot(a(1:2),[1 1],'--k');
hold off
xlabel('Time (s)'); ylabel('SI');
title('Mixed')
set(gca,'TickDir','out'); box off
legend;
h=get(gca,'Legend');
h.String=h.String(1:n-1);
h.String{1}='Closest'; h.String{end}='Furthest';
for i=2:length(h.String)-1
    h.String{i}='';
end

sax= round(linspace(1,length(dix_mix),11));
for i=1:10
    DD(i)= nanmean(mix_dist(dix_mix(sax(i):sax(i+1))));
    avg_max(i)= nanmean(nansum(mixers(dix_mix(sax(i):sax(i+1)),find(ax>-.1 & ax<.2)),2));
    smax(i)= nanmean(nansum(mixers(dix_mix(sax(i):sax(i+1)),find(ax>-.1 & ax<.2)),2))./sqrt(sax(i+1)-sax(i)+1);
end

subplot 247
hold on
errorbar(DD,avg_max,smax,'-s','MarkerFaceColor','g','MarkerEdgeColor','g','Color','g','CapSize',0,'LineWidth',2);
hold off
xlabel('Avg. Dist');
ylabel('\SigmaSI (-0.1 - 0.2 s)');
set(gca,'TickDir','out'); box off

[RR,PP]= corr(mix_dist(:),mixers);
subplot(6,4,24)
plot(ax,RR,'g','LineWidth',2)
hold on
scatter(ax(find(PP<.05/sm)),RR(find(PP<.05/sm)),30,'o','MarkerFaceColor','g','MarkerEdgeColor','none')
plot([ax(1) ax(end)],[0 0],'--k');
a=axis;
plot([0 0],a(3:4),'--k');
hold off
set(gca,'TickDir','out'); box off
xlabel('Time (s)');
xlim([-tbf+0.01 taf-0.01])
%% per mouse
sm=31;
% all_sync= orig_all_sync;
% all_sync_random= orig_all_sync_random;

um= unique(all_mouse_id);
for i=1:length(um)
    ix= find(strcmpi(all_mouse_id',um{i}) & gain(:,1)<-.05 & gain(:,2)<-.05);
    subplot(3,length(um)+1,i)
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
    title(['Pausers; Mouse #' num2str(i) '; n=' num2str(length(ix))])
    all_pausers(i,:)=smooth(nanmean(all_sync(ix,:),1),sm);
    all_pausers_random(i,:)=smooth(nanmean(all_sync_random(ix,:),1),sm);
    
    ix= find(strcmpi(all_mouse_id',um{i}) & gain(:,1)>.05 & gain(:,2)>.05);
    subplot(3,length(um)+1,i+length(um)+1)
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
    title(['Bursters; Mouse #' num2str(i) '; n=' num2str(length(ix))])
    all_Bursters(i,:)=smooth(nanmean(all_sync(ix,:),1),sm);
    all_bursters_random(i,:)=smooth(nanmean(all_sync_random(ix,:),1),sm);
    
    ix= find(strcmpi(all_mouse_id',um{i}) & ((gain(:,1)<-.05 & gain(:,2)>.05) | (gain(:,1)>.05 & gain(:,2)<-.05)));
    subplot(3,length(um)+1,i+2*(length(um)+1))
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
    title(['Mixed; Mouse #' num2str(i) '; n=' num2str(length(ix))])
    all_mixed(i,:)=smooth(nanmean(all_sync(ix,:),1),sm);
    all_mixed_random(i,:)=smooth(nanmean(all_sync_random(ix,:),1),sm);
end
%%
tix= find(ax>-.3 & ax<.3);
figure
sm=21;
subplot(2,3,1)
m=smooth(nanmean(all_pausers,1),sm); m=m(:)';
s=smooth(nanstd(all_pausers,[],1),sm)./sqrt(size(all_pausers,1)); s=s(:)';
xx=[ax ax(end) ax(end) fliplr(ax) ax(1) ax(1)];
yy=[m-s m(end)-s(end) m(end)+s(end) fliplr(m+s) m(1)+s(1) m(1)-s(1)];
patch(xx,yy,[0.5 0.5 1],'EdgeColor','none');
hold on
plot(ax,smooth(nanmean(all_pausers,1),sm),'Color','b','LineWidth',1)
a=axis;
plot(a(1:2),[1 1],'--k')
plot([0 0],a(3:4),'--k')
hold off
xlabel('Time (s)'); ylabel('SI');
subplot(2,3,4)
plot([1 2], [nanmedian(all_pausers_random(:,tix),2)./nanmedian(all_pausers_random(:,tix),2) nanmedian(all_pausers(:,tix),2)./nanmedian(all_pausers_random(:,tix),2)])

subplot(2,3,2)
m=smooth(nanmean(all_Bursters,1),sm); m=m(:)';
s=smooth(nanstd(all_Bursters,[],1),sm)./sqrt(size(all_Bursters,1)); s=s(:)';
xx=[ax ax(end) ax(end) fliplr(ax) ax(1) ax(1)];
yy=[m-s m(end)-s(end) m(end)+s(end) fliplr(m+s) m(1)+s(1) m(1)-s(1)];
patch(xx,yy,[1 0.5 0.5],'EdgeColor','none');
hold on
plot(ax,smooth(nanmean(all_Bursters,1),sm),'Color','r','LineWidth',1)
a=axis;
plot(a(1:2),[1 1],'--k')
plot([0 0],a(3:4),'--k')
hold off
xlabel('Time (s)'); ylabel('SI');
subplot(2,3,5)
plot([1 2], [nanmean(all_bursters_random(:,tix),2)./nanmean(all_bursters_random(:,tix),2) nanmean(all_Bursters(:,tix),2)./nanmean(all_bursters_random(:,tix),2)])

subplot(2,3,3)
m=smooth(nanmean(all_mixed,1),sm); m=m(:)';
s=smooth(nanstd(all_mixed,[],1),sm)./sqrt(size(all_mixed,1)); s=s(:)';
xx=[ax ax(end) ax(end) fliplr(ax) ax(1) ax(1)];
yy=[m-s m(end)-s(end) m(end)+s(end) fliplr(m+s) m(1)+s(1) m(1)-s(1)];
patch(xx,yy,[1 0.5 1],'EdgeColor','none');
hold on
plot(ax,smooth(nanmean(all_mixed,1),sm),'Color','m','LineWidth',1)
a=axis;
plot(a(1:2),[1 1],'--k')
plot([0 0],a(3:4),'--k')
hold off
xlabel('Time (s)'); ylabel('SI');
subplot(2,3,6)
plot([1 2], [nanmean(all_mixed_random(:,tix),2)./nanmean(all_mixed_random(:,tix),2) nanmean(all_mixed(:,tix),2)./nanmean(all_mixed_random(:,tix),2)])
%%
N=10;
ix= find(gain(:,1)<-.05 & gain(:,2)<-.05);
pausers_cov=all_cov(ix,1);
[~,cov_ix]= sort(pausers_cov);
D= round(1:length(ix)/N:length(ix));
D=[D length(ix)];
for i=1:length(D)-1
    INDEX=ix(cov_ix(D(i):D(i+1)));
    all_2Cov(i,:)= nanmean(all_sync(INDEX,:),1);
    avg_cov= nanmean(pausers_cov(cov_ix(D(i):D(i+1))))
end
%%
ix= find(gain(:,1)<-.05 & gain(:,2)<-.05);
pausers_SI= all_sync(ix,:);
for i=1:size(pausers_SI,1)
    pausers_SI(i,:)= smooth(pausers_SI(i,:),11);
end
%%
[coeff,score,latent] = pca(pausers_SI');
subplot(2,4,1)
norm_si= (pausers_SI-repmat(nanmean(pausers_SI(:,find(ax<-.3)),2),1,length(ax)));
norm_si=norm_si./repmat(max(norm_si,[],2),1,length(ax));
[~,ix]= sort(max(norm_si,[],2));
imagesc(ax,1:size(pausers_SI,1),norm_si(ix,:))
a=axis;
hold on
plot([0 0],a(3:4),'--w')
colormap jet
xlabel('Time to endpoint (s)');
ylabel('Pair #');
set(gca,'TickDir','out'); box off
caxis([-1 1]);

subplot(2,4,2)
plot(1:length(latent),cumsum(latent)./sum(latent),'Color','k','LineWidth',2);
set(gca,'TickDir','out'); box off
xlabel('PC #'); 
ylabel('Var. exp');

k=5;
npc=5;
subplot(2,4,3)
h=fitgmdist(coeff(:,1:npc),k);
idx= cluster(h,coeff(:,1:npc));
for i=1:k
    ix= find(idx==i);
    scatter3(coeff(ix,1),coeff(ix,2),coeff(ix,3),'.')
    hold on
end
hold off
xlabel('Coeff 1'); ylabel('Coeff 2'); zlabel('Coeff 3');

subplot(2,4,4)
for i=1:npc
    plot(ax,score(:,i),'LineWidth',2)
    hold on
end
hold off
set(gca,'TickDir','out'); box off
xlabel('Time (s)');
ylabel('a.u.')

subplot(2,4,5)
for i=1:k
    ix= find(idx==i);
    X=nanmean(pausers_SI(ix,:),1);
    X=X-nanmean(X(find(ax<-.3)));
    X=X./max(X);
    plot(ax,X,'LineWidth',2)
    hold on
    Ls{i}= num2str(length(ix));
end
a=axis;
plot([0 0],a(3:4),'--k')
hold off
legend(Ls);
xlabel('Time to endpoint (s)');
ylabel('Norm. SI');
set(gca,'TickDir','out'); box off
%
subplot(2,4,6)
m=[]; s=[];
for i=1:k
    ix= find(idx==i);
    m= nanmean(coeff(ix,1:npc));
    s= nanstd(coeff(ix,1:npc))./sqrt(length(ix));
    errorbar(1:npc,m,s)
    hold on
end
hold off
xlabel('PC#');
ylabel('average coeff (a.u.)');
set(gca,'TickDir','out'); box off
%% k-means
k=5;
[idx,C]=kmeans(pausers_SI,k);
subplot(2,4,7)
plot(ax,C,'LineWidth',2)
xlabel('Time (s)');
ylabel('a.u.');
title('Centroids; kmeans');

clear L
subplot(2,4,8)
for i=1:k
    ix= find(idx==i);
    X=nanmean(pausers_SI(ix,:),1);
    X=X-nanmean(X(find(ax<-.3)));
    X=X./max(X);
    plot(ax,X,'LineWidth',2)
    L{i}=num2str(length(ix));
    hold on
end
legend(L)
a=axis;
plot([0 0],a(3:4),'--k');
hold off
xlabel('Time (s)');
ylabel('Norm. SI');
title('Clusters; kmeans');