%
udir= '../new_data/';
dates = dir([udir 'D*']);
bin2use=1;
addpath ../helper_functions/
tbf=1; taf=1;
win2corr=.2; d2corr= .020; ax2corr=-tbf:d2corr:taf-win2corr;
ax=-tbf*1000:bin2use:taf*1000; ax=ax(1:end-1); ax=ax./1000;
c_all=0;
c_allv=0;
c_alla=0;
bins2tune=50;
n=50;
class2take={'PC'};
freq= 120;
L= ((taf+tbf)*freq);
move_thresh=.2;
%%
for j=1:length(dates)
    load([udir dates(j).name]);
    all_x=[]; all_y=[]; all_z=[];
    all_xv=[]; all_yv=[]; all_zv=[];
    all_xacc=[]; all_yacc=[]; all_zacc=[]; this_S=[];
    for i=1:length(ReachS)
        xyz= ReachS(i).filt_kin(:,2:4);
        v= ReachS(i).filt_kin(:,5);
        xv= ReachS(i).filt_kin(:,6);
        yv= ReachS(i).filt_kin(:,7);
        zv= ReachS(i).filt_kin(:,8);
        times = ReachS(i).filt_kin(:,1);
        ix= find(ReachS(i).out(:,2)>move_thresh,1,'first');
        endpoint_time=ReachS(i).out(ix,1);
        if isempty(endpoint_time)
            ix=find(ReachS(i).filt_kin(:,2)>move_thresh,1,'first');
            endpoint_time=ReachS(i).filt_kin(ix,1);
        end
        times= ReachS(i).filt_kin(:,1);
        P1= repmat(ReachS(i).filt_kin(1,2:4),size(ReachS(i).filt_kin,1),1);
        
        trial_x= xyz(find(times>=endpoint_time-tbf & times<=endpoint_time+taf),1);
        trial_y= xyz(find(times>=endpoint_time-tbf & times<=endpoint_time+taf),2);
        trial_z= xyz(find(times>=endpoint_time-tbf & times<=endpoint_time+taf),3);
        tt=times(find(times>=endpoint_time-tbf & times<=endpoint_time+taf));
        
        trial_v= v(find(times>=endpoint_time-tbf & times<=endpoint_time+taf));
        trial_xv= xv(find(times>=endpoint_time-tbf & times<=endpoint_time+taf));
        trial_yv= yv(find(times>=endpoint_time-tbf & times<=endpoint_time+taf));
        trial_zv= zv(find(times>=endpoint_time-tbf & times<=endpoint_time+taf));
        
        
        
        vStimMode = isfield(ReachS(i),'stim');
        if vStimMode, vStim=ReachS(i).stim; else, vStim=0; end
        if isempty(vStim), vStim=0; end
        vExclude=~ReachS(i).exclude;
        if length(trial_x)>L*0.6
            t=tt; d=diff(t);
            xtmp=trial_x;
            ytmp=trial_y;
            ztmp=trial_z;
            tq=[t(1):1/freq:t(end)];
            xx=interp1(t,xtmp,tq);
            yy=interp1(t,ytmp,tq);
            zz=interp1(t,ztmp,tq);
            xx=[xx nan(1,L-length(xx))];
            yy=[yy nan(1,L-length(yy))];
            zz=[zz nan(1,L-length(zz))];
            all_x(end+1,:)= xx;
            all_y(end+1,:)= yy;
            all_z(end+1,:)= zz;
            
            xtmp=trial_xv;
            ytmp=trial_yv;
            ztmp=trial_zv;
            xx=interp1(t,xtmp,tq);
            yy=interp1(t,ytmp,tq);
            zz=interp1(t,ztmp,tq);
            
            xx=[xx nan(1,L-length(xx))];
            yy=[yy nan(1,L-length(yy))];
            zz=[zz nan(1,L-length(zz))];
            all_xv(end+1,:)= xx;
            all_yv(end+1,:)= yy;
            all_zv(end+1,:)= zz;
            this_S(end+1)=vStim;
            all_xacc(i,:)= [0 diff(xx)];
            all_yacc(i,:)= [0 diff(yy)];
            all_zacc(i,:)= [0 diff(zz)];
            
            c_this=0;
            for cc=1:length(cellData)
                c_this=c_this+1;
                this_trc=cellData(cc).Bin10smooth;
                index= find(this_trc(:,1)>=endpoint_time-tbf & this_trc(:,1)<=endpoint_time+taf);
                FR(c_this,i,:)= this_trc(index,2);
                spike_times= cellData(cc).spikeTimes;
                ISI = diff(spike_times);
                m=[];
                for isi_count=1:length(ISI)-1
                    m(isi_count)= abs(ISI(isi_count)-ISI(isi_count+1))/(ISI(isi_count)+ISI(isi_count+1));
                end
                cv2(c_this)=nanmean(m);
                cv(c_this)=nanstd(ISI)/nanmean(ISI);
                if isfield(cellData,'CS_Bin1')
                    CSon(c_this)=~isempty(cellData(cc).CS_Bin1);
                else
                    CSon(c_this)=0;
                end
                gain(c_this)=cellData(cc).gain;
            end
        end
    end
    for i=1:size(FR,1)
        c_all=c_all+1;
        all_cells(c_all,:)= squeeze(nanmean(FR(i,:,:),2));
        all_CSperCells(c_all)=CSon(i);
        all_bhv(1,c_all,:)= nanmean(all_x,1);
        all_bhv(2,c_all,:)= nanmean(all_y,1);
        all_bhv(3,c_all,:)= nanmean(all_z,1);
        all_cv(c_all)=cv(i);
        all_cv2(c_all)=cv2(i);
        all_gain(c_all)=gain(i);
    end
    clear FR all_x all_y all_z this_S CSon cv cv2 gain
end
%%
nax= [-tbf:1/100:taf];nax=nax(1:end-1);
bax=[-tbf:1/freq:taf]; bax=bax(1:end-1);
%%
t1= find(nax<-0.5); t2= find(nax>0 & nax<.5);
for i=1:size(all_cells,1)
    fr1= nanmean(all_cells(i,t1)); fr2= nanmean(all_cells(i,t2));
    gain(i)= (fr2-fr1)/(fr1+fr2);
end
[~,ix]= sort(gain);
norm_all_cells= all_cells./repmat(max(all_cells,[],2),1,size(all_cells,2));

bursters= norm_all_cells(find(gain>=.05),:);
pausers= norm_all_cells(find(gain<=-.05),:);
mix= norm_all_cells(find(gain>=-.05 & gain<=.05),:);
[~,tmp]= max(bursters,[],2); [~,ix]= sort(tmp);
bursters=bursters(ix,:);
[~,tmp]= min(pausers,[],2); [~,ix]= sort(tmp);
pausers=pausers(ix,:);
[~,tmp]= max(mix,[],2); [~,ix]= sort(tmp);
mix=mix(ix,:);
new_norm_all_cells=[bursters; mix; pausers];
subplot 221
imagesc(nax,1:size(all_cells,1),new_norm_all_cells);
xlabel('Time(s)');
ylabel('Cell #');
hold on
plot([nax(1) nax(end)],[size(bursters,1) size(bursters,1)],'-r','LineWidth',2)
plot([nax(1) nax(end)],[length(gain)-size(pausers,1) length(gain)-size(pausers,1)],'-b','LineWidth',2)
hold off
set(gca,'TickDir','out'); box off

subplot 443
y=nanmean(norm_all_cells(find(gain<=-.05),:),1);
s=nanstd(norm_all_cells(find(gain<=-.05),:),1)./sqrt(length(find(gain<=-.05)));
xx=[nax nax(end) nax(end) fliplr(nax) nax(1) nax(1)];
yy=[y-s y(end)-s(end) y(end)+s(end) fliplr(y+s) y(1)+s(1) y(1)-s(1)];
patch(xx,yy,[0.6 0.6 1],'EdgeColor','none')
hold on
y=nanmean(norm_all_cells(find(gain>=.05),:),1);
s=nanstd(norm_all_cells(find(gain>=.05),:),1)./sqrt(length(find(gain>=.05)));
yy=[y-s y(end)-s(end) y(end)+s(end) fliplr(y+s) y(1)+s(1) y(1)-s(1)];
patch(xx,yy,[1 0.6 0.6],'EdgeColor','none')
y=nanmean(norm_all_cells(find(gain<=.05 & gain>=-.05),:),1);
s=nanstd(norm_all_cells(find(gain<=.05 & gain>=-.05),:),1)./sqrt(length(find(gain<=.05 & gain>=-.05)));
yy=[y-s y(end)-s(end) y(end)+s(end) fliplr(y+s) y(1)+s(1) y(1)-s(1)];
patch(xx,yy,[.6 0.6 0.6],'EdgeColor','none')
plot(nax,nanmean(norm_all_cells(find(gain<=-.05),:),1),'b','LineWidth',2)
plot(nax,nanmean(norm_all_cells(find(gain>=.05),:),1),'r','LineWidth',2)
plot(nax,nanmean(norm_all_cells(find(gain<=.05 & gain>=-.05),:),1),'k','LineWidth',2)
a=axis;
plot([0 0],[a(3) a(4)],'--k');
hold off
ylabel('normalized Fr'); xlabel('Time from MO (s)');
set(gca,'TickDir','out'); box off

subplot(4,4,7)
plot(bax,squeeze(nanmean(all_bhv(1,:,:),2)),'LineWidth',1)
hold on
plot(bax,squeeze(nanmean(all_bhv(2,:,:),2)),'LineWidth',1)
plot(bax,squeeze(nanmean(all_bhv(3,:,:),2)),'LineWidth',1)
hold off
ylabel('Position (cm)')
xlabel('Time (s)')

subplot 448
scatter3(squeeze(nanmean(all_bhv(1,:,:),2)),squeeze(nanmean(all_bhv(3,:,:),2)),squeeze(nanmean(all_bhv(2,:,:),2)),50,bax,'fill')
xlabel('X (cm)')
ylabel('z (cm)')
zlabel('y (cm)')

subplot 444
y=nanmean(all_cells(find(gain<=-.05),:),1);
s=nanstd(all_cells(find(gain<=-.05),:),1)./sqrt(length(find(gain<=-.05)));
xx=[nax nax(end) nax(end) fliplr(nax) nax(1) nax(1)];
yy=[y-s y(end)-s(end) y(end)+s(end) fliplr(y+s) y(1)+s(1) y(1)-s(1)];
patch(xx,yy,[0.6 0.6 1],'EdgeColor','none')
hold on
y=nanmean(all_cells(find(gain>=.05),:),1);
s=nanstd(all_cells(find(gain>=.05),:),1)./sqrt(length(find(gain>=.05)));
yy=[y-s y(end)-s(end) y(end)+s(end) fliplr(y+s) y(1)+s(1) y(1)-s(1)];
patch(xx,yy,[1 0.6 0.6],'EdgeColor','none')
y=nanmean(all_cells(find(gain<.05 & gain>-.05),:),1);
s=nanstd(all_cells(find(gain<.05 & gain>-.05),:),1)./sqrt(length(find(gain<.05 & gain>-.05)));
yy=[y-s y(end)-s(end) y(end)+s(end) fliplr(y+s) y(1)+s(1) y(1)-s(1)];
patch(xx,yy,[.6 0.6 0.6],'EdgeColor','none')
plot(nax,nanmean(all_cells(find(gain<=-.05),:),1),'b','LineWidth',2)
plot(nax,nanmean(all_cells(find(gain>=.05),:),1),'r','LineWidth',2)
plot(nax,nanmean(all_cells(find(gain<.05 & gain>-.05),:),1),'k','LineWidth',2)
a=axis;
plot([0 0],[a(3) a(4)],'--k');
hold off
ylabel('Fr (sp/s)'); xlabel('Time from MO (s)');
set(gca,'TickDir','out'); box off
%%
base_ix= find(nax<-.5 & nax>-1);
move_ix= find(nax>-.2 & nax<.3);
FR_base=nanmean(all_cells(:,base_ix),2);
FR_move=nanmean(all_cells(:,move_ix),2);
pause_ix= find(gain<=-.05);
burst_ix= find(gain>=.05);
mix_ix= find(gain<.05 & gain>-.05);

subplot(2,6,7)
bar(1,nanmean(FR_base(pause_ix,:)),'FaceColor','b','EdgeColor','k')
hold on
scatter(1+rand(1,length(pause_ix))./10,FR_base(pause_ix,:),20,'ob','fill','MarkerEdgeColor','w')
bar(2,nanmean(FR_base(burst_ix,:)),'FaceColor','r','EdgeColor','k')
scatter(2+rand(1,length(burst_ix))./10,FR_base(burst_ix,:),20,'or','fill','MarkerEdgeColor','w')
bar(3,nanmean(FR_base(mix_ix,:)),'FaceColor','k','EdgeColor','k')
scatter(3+rand(1,length(mix_ix))./10,FR_base(mix_ix,:),20,'ok','fill','MarkerEdgeColor','w')
hold off
set(gca,'TickDir','out'); box off
ylabel('Baseline FR (sp/s)')
set(gca,'XTick',1:3,'XTickLabel',{'Pausers' 'Bursters' 'Mixed'})

subplot(2,6,8)
bar(1,nanmean(FR_move(pause_ix,:)),'FaceColor','b','EdgeColor','k')
hold on
scatter(1+rand(1,length(pause_ix))./10,FR_move(pause_ix,:),20,'ob','fill','MarkerEdgeColor','w')
bar(2,nanmean(FR_move(burst_ix,:)),'FaceColor','r','EdgeColor','k')
scatter(2+rand(1,length(burst_ix))./10,FR_move(burst_ix,:),20,'or','fill','MarkerEdgeColor','w')
bar(3,nanmean(FR_move(mix_ix,:)),'FaceColor','k','EdgeColor','k')
scatter(3+rand(1,length(mix_ix))./10,FR_move(mix_ix,:),20,'ok','fill','MarkerEdgeColor','w')
hold off
set(gca,'TickDir','out'); box off
ylabel('Move FR (sp/s)')
set(gca,'XTick',1:3,'XTickLabel',{'Pausers' 'Bursters' 'Mixed'})

subplot(2,6,9)
bar(1,nanmean(all_cv(pause_ix)),'FaceColor','b','EdgeColor','k')
hold on
scatter(1+rand(1,length(pause_ix))./10,all_cv(pause_ix),20,'ob','fill','MarkerEdgeColor','w')
bar(2,nanmean(all_cv(burst_ix)),'FaceColor','r','EdgeColor','k')
scatter(2+rand(1,length(burst_ix))./10,all_cv(burst_ix),20,'or','fill','MarkerEdgeColor','w')
bar(3,nanmean(all_cv(mix_ix)),'FaceColor','k','EdgeColor','k')
scatter(3+rand(1,length(mix_ix))./10,all_cv(mix_ix),20,'ok','fill','MarkerEdgeColor','w')
hold off
set(gca,'TickDir','out'); box off
ylabel('CV')
set(gca,'XTick',1:3,'XTickLabel',{'Pausers' 'Bursters' 'Mixed'})

subplot(2,6,10)
bar(1,nanmean(all_cv2(pause_ix)),'FaceColor','b','EdgeColor','k')
hold on
scatter(1+rand(1,length(pause_ix))./10,all_cv2(pause_ix),20,'ob','fill','MarkerEdgeColor','w')
bar(2,nanmean(all_cv2(burst_ix)),'FaceColor','r','EdgeColor','k')
scatter(2+rand(1,length(burst_ix))./10,all_cv2(burst_ix),20,'or','fill','MarkerEdgeColor','w')
bar(3,nanmean(all_cv2(mix_ix)),'FaceColor','k','EdgeColor','k')
scatter(3+rand(1,length(mix_ix))./10,all_cv2(mix_ix),20,'ok','fill','MarkerEdgeColor','w')
hold off
set(gca,'TickDir','out'); box off
ylabel('CV2')
set(gca,'XTick',1:3,'XTickLabel',{'Pausers' 'Bursters' 'Mixed'})

subplot(2,6,11)
bar(1,nanmean(abs(FR_move(pause_ix,:)-FR_base(pause_ix,:))),'FaceColor','b','EdgeColor','k')
hold on
scatter(1+rand(1,length(pause_ix))./10,abs(FR_move(pause_ix,:)-FR_base(pause_ix,:)),20,'ob','fill','MarkerEdgeColor','w')
bar(2,nanmean(abs(FR_base(burst_ix,:)-FR_move(burst_ix,:))),'FaceColor','r','EdgeColor','k')
scatter(2+rand(1,length(burst_ix))./10,abs(FR_move(burst_ix,:)-FR_base(burst_ix,:)),20,'or','fill','MarkerEdgeColor','w')
bar(3,nanmean(abs(FR_move(mix_ix,:)-FR_base(mix_ix,:))),'FaceColor','k','EdgeColor','k')
scatter(3+rand(1,length(mix_ix))./10,abs(FR_move(mix_ix,:)-FR_base(mix_ix,:)),20,'ok','fill','MarkerEdgeColor','w')
hold off
set(gca,'TickDir','out'); box off
ylabel('Modulation (sp/s)')
set(gca,'XTick',1:3,'XTickLabel',{'Pausers' 'Bursters' 'Mixed'})

subplot(2,6,12)
histogram(all_gain,50)