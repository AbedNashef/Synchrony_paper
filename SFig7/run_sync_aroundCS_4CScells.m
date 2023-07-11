%
clear
udir= '../new_data/';
dates = dir([udir 'D*']);
bin2use=1;
addpath ../helper_functions/
tbf=0.05; taf=0.1;
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
    if isfield(cellData,'CS_Bin1')
        cell2take=[];
        for i=1:length(cellData)
            if ~isempty(cellData(i).CS_Bin1)
                cell2take=[cell2take i];
            end
        end
        if length(cell2take)>1
            for i= 1:length(cell2take)
                ix =find(cellData(cell2take(i)).CS_Bin1(:,2));
                times= cellData(cell2take(i)).CS_Bin1(ix,1);
                
                
                tt= repmat(ix,1,length(ax))+repmat(ax.*1000,length(times),1);
                max_l= size(cellData(cell2take(i)).CS_Bin1,1);
                tt=tt(find(tt(:,1)>0 & tt(:,end)<max_l),:);
                c_this=0;
                for cc=1:length(cell2take)
                    c_this=c_this+1;
                    this_trc=cellData(cell2take(cc)).Bin1(:,2);
                    FR(c_this,:,:)= this_trc(tt);
                    this_gain(c_this)=cellData(cell2take(cc)).gain;
                    Chs{c_this}= cellData(cc).Channels;
                    
                end
                
                for jj=1:size(FR,1)-1
                    for ii=jj+1:size(FR,1)
                        ch1=Chs{jj};
                        ch2=Chs{ii};
                        if ~isempty(ch1) & ~isempty(ch2)
                            overlap= length(intersect(ch1,ch2))/length(unique([ch1,ch2]));
                        else
                            overlap=0;
                        end
                        if overlap==0
                            c_all=c_all+1;
                            cell1= squeeze(FR(jj,:,:));
                            cell2= squeeze(FR(ii,:,:));
                            
                            Pr_S1= nanmean(cell1,1);
                            Pr_S2= nanmean(cell2,1);
                            Pr_S1S2= nanmean(cell1.*cell2,1);
                            
                            all_sync(c_all,:)= Pr_S1S2./(Pr_S1.*Pr_S2);
                            
                            
                            Pr_S1= nanmean(cell1,1);
                            Pr_S2= nanmean(cell2,1);
                            [~,ix1]=sort(rand(1,size(cell1,1)));
                            [~,ix2]=sort(rand(1,size(cell2,1)));
                            Pr_S1S2= nanmean(cell1(ix1,:).*cell2(ix2,:),1);
                            
                            all_sync_random(c_all,:)= Pr_S1S2./(Pr_S1.*Pr_S2);
                            [r_all(c_all,1),p_all(c_all,1)]= corr(Pr_S1(:),all_sync(c_all,:)');
                            [r_all(c_all,2),p_all(c_all,2)]= corr(Pr_S2(:),all_sync(c_all,:)');
                            gain(c_all,:)=this_gain([jj,ii]);
                            
                            all_Cell1(c_all,:)= nanmean(cell1,1);
                            all_Cell2(c_all,:)= nanmean(cell2,1);
                            date_id(c_all)= j;
                            Dist4cells(c_all)= abs(cellData(jj).depth-cellData(ii).depth);
                        end
                    end
                end
                clear FR FR2gain this_S nX nY nZ this_gain
            end
            
        end
    end
end

%%
sm=5;
gwin=gausswin(sm)'; gwin=gwin./sum(gwin);
axx=[-tbf:nbin/1000:taf]; axx=axx(1:end-1);
subplot(2,2,1)
ix= find(gain(:,1)<-.05 & gain(:,2)<-.05);
% m=conv2(nanmean(all_sync(ix,:),1),gwin,'same'); m=m(:)';
% s=conv2(nanstd(all_sync(ix,:),1),gwin,'same')./sqrt(length(ix)); s=s(:)';
sm=21;
m=smooth(nanmean(all_sync(ix,:),1),sm); m=m(:)';
s=smooth(nanstd(all_sync(ix,:),[],1),sm)./sqrt(length(ix)); s=s(:)';
xx=[ax ax(end) ax(end) fliplr(ax) ax(1) ax(1)];
yy=[m-s m(end)-s(end) m(end)+s(end) fliplr(m+s) m(1)+s(1) m(1)-s(1)];
patch(xx,yy,[0.5 0.5 1],'EdgeColor','none');
hold on
% m=conv2(nanmean(all_sync_random(ix,:),1),gwin,'same'); m=m(:)';
% s=conv2(nanstd(all_sync_random(ix,:),1),gwin,'same')./sqrt(length(ix)); s=s(:)';
m=smooth(nanmean(all_sync_random(ix,:),1),sm); m=m(:)';
s=smooth(nanstd(all_sync_random(ix,:),[],1),sm)./sqrt(length(ix)); s=s(:)';
xx=[ax ax(end) ax(end) fliplr(ax) ax(1) ax(1)];
yy=[m-s m(end)-s(end) m(end)+s(end) fliplr(m+s) m(1)+s(1) m(1)-s(1)];
patch(xx,yy,[0.5 0.5 0.5],'EdgeColor','none');
plot(ax,conv2(nanmean(all_sync(ix,:),1),gwin,'same'),'Color','b','LineWidth',1)
plot(ax,conv2(nanmean(all_sync_random(ix,:),1),gwin,'same'),'Color','k','LineWidth',1)
a=axis;
plot([0 0],[a(3) a(4)],'--k','LineWidth',2)
hold off
set(gca,'TickDir','out'); box off
ylabel('SI'); xlabel('Time (s)')
title(['Pausers; n=' num2str(length(ix))])

subplot(2,2,2)
ix= find(gain(:,1)>.05 & gain(:,2)>.05);
% m=conv2(nanmean(all_sync(ix,:),1),gwin,'same'); m=m(:)';
% s=conv2(nanstd(all_sync(ix,:),1),gwin,'same')./sqrt(length(ix)); s=s(:)';
m=smooth(nanmean(all_sync(ix,:),1),sm); m=m(:)';
s=smooth(nanstd(all_sync(ix,:),[],1),sm)./sqrt(length(ix)); s=s(:)';
xx=[ax ax(end) ax(end) fliplr(ax) ax(1) ax(1)];
yy=[m-s m(end)-s(end) m(end)+s(end) fliplr(m+s) m(1)+s(1) m(1)-s(1)];
patch(xx,yy,[1 0.5 0.5],'EdgeColor','none');
hold on
% m=conv2(nanmean(all_sync_random(ix,:),1),gwin,'same'); m=m(:)';
% s=conv2(nanstd(all_sync_random(ix,:),1),gwin,'same')./sqrt(length(ix)); s=s(:)';
m=smooth(nanmean(all_sync_random(ix,:),1),sm); m=m(:)';
s=smooth(nanstd(all_sync_random(ix,:),[],1),sm)./sqrt(length(ix)); s=s(:)';
xx=[ax ax(end) ax(end) fliplr(ax) ax(1) ax(1)];
yy=[m-s m(end)-s(end) m(end)+s(end) fliplr(m+s) m(1)+s(1) m(1)-s(1)];
patch(xx,yy,[0.5 0.5 0.5],'EdgeColor','none');
plot(ax,conv2(nanmean(all_sync(ix,:),1),gwin,'same'),'Color','r','LineWidth',1)
plot(ax,conv2(nanmean(all_sync_random(ix,:),1),gwin,'same'),'Color','k','LineWidth',1)
a=axis;
plot([0 0],[a(3) a(4)],'--k','LineWidth',2)
hold off
set(gca,'TickDir','out'); box off
ylabel('SI'); xlabel('Time (s)')
title(['Bursters; n=' num2str(length(ix))])
subplot(2,2,3)
ix= find((gain(:,1)>.05 & gain(:,2)<-.05) | (gain(:,1)<-.05 & gain(:,2)>.05));
% m=conv2(nanmean(all_sync(ix,:),1),gwin,'same'); m=m(:)';
% s=conv2(nanstd(all_sync(ix,:),1),gwin,'same')./sqrt(length(ix)); s=s(:)';
m=smooth(nanmean(all_sync(ix,:),1),sm); m=m(:)';
s=smooth(nanstd(all_sync(ix,:),[],1),sm)./sqrt(length(ix)); s=s(:)';
xx=[ax ax(end) ax(end) fliplr(ax) ax(1) ax(1)];
yy=[m-s m(end)-s(end) m(end)+s(end) fliplr(m+s) m(1)+s(1) m(1)-s(1)];
patch(xx,yy,[0.5 1 0.5],'EdgeColor','none');
hold on
% m=conv2(nanmean(all_sync_random(ix,:),1),gwin,'same'); m=m(:)';
% s=conv2(nanstd(all_sync_random(ix,:),1),gwin,'same')./sqrt(length(ix)); s=s(:)';
m=smooth(nanmean(all_sync_random(ix,:),1),sm); m=m(:)';
s=smooth(nanstd(all_sync_random(ix,:),[],1),sm)./sqrt(length(ix)); s=s(:)';
xx=[ax ax(end) ax(end) fliplr(ax) ax(1) ax(1)];
yy=[m-s m(end)-s(end) m(end)+s(end) fliplr(m+s) m(1)+s(1) m(1)-s(1)];
patch(xx,yy,[0.5 0.5 0.5],'EdgeColor','none');
plot(ax,conv2(nanmean(all_sync(ix,:),1),gwin,'same'),'Color','g','LineWidth',1)
plot(ax,conv2(nanmean(all_sync_random(ix,:),1),gwin,'same'),'Color','k','LineWidth',1)
a=axis;
plot([0 0],[a(3) a(4)],'--k','LineWidth',2)
hold off
set(gca,'TickDir','out'); box off
ylabel('SI'); xlabel('Time (s)')
title(['Bursters-Pausers; n=' num2str(length(ix))])


subplot(2,2,4)
ix= find((gain(:,1)<-.05 & gain(:,2)<-.05));
m=conv2(nanmean(all_sync(ix,:),1),gwin,'same'); m=m(:)';
s=conv2(nanstd(all_sync(ix,:),1),gwin,'same')./sqrt(length(ix)); s=s(:)';
xx=[ax ax(end) ax(end) fliplr(ax) ax(1) ax(1)];
yy=[m-s m(end)-s(end) m(end)+s(end) fliplr(m+s) m(1)+s(1) m(1)-s(1)];
patch(xx,yy,[0.5 0.5 1],'EdgeColor','none');
hold on
plot(ax,conv2(nanmean(all_sync(ix,:),1),gwin,'same'),'Color','b','LineWidth',1)
ix= find((gain(:,1)>.05 & gain(:,2)>.05));
m=conv2(nanmean(all_sync(ix,:),1),gwin,'same'); m=m(:)';
s=conv2(nanstd(all_sync(ix,:),1),gwin,'same')./sqrt(length(ix)); s=s(:)';
xx=[ax ax(end) ax(end) fliplr(ax) ax(1) ax(1)];
yy=[m-s m(end)-s(end) m(end)+s(end) fliplr(m+s) m(1)+s(1) m(1)-s(1)];
patch(xx,yy,[1 0.5 0.5],'EdgeColor','none');
plot(ax,conv2(nanmean(all_sync(ix,:),1),gwin,'same'),'Color','r','LineWidth',1)

ix= find((gain(:,1)>.05 & gain(:,2)<-.05) | (gain(:,1)<-.05 & gain(:,2)>.05));
m=conv2(nanmean(all_sync(ix,:),1),gwin,'same'); m=m(:)';
s=conv2(nanstd(all_sync(ix,:),1),gwin,'same')./sqrt(length(ix)); s=s(:)';
xx=[ax ax(end) ax(end) fliplr(ax) ax(1) ax(1)];
yy=[m-s m(end)-s(end) m(end)+s(end) fliplr(m+s) m(1)+s(1) m(1)-s(1)];
patch(xx,yy,[0.5 1 0.5],'EdgeColor','none');
plot(ax,conv2(nanmean(all_sync(ix,:),1),gwin,'same'),'Color','g','LineWidth',1)
a=axis;
plot([0 0],[a(3) a(4)],'--k','LineWidth',2)
hold off
set(gca,'TickDir','out'); box off
ylabel('SI'); xlabel('Time (s)')

%%
base_ix= find(ax<-.01);
sm=5;
gwin=gausswin(sm)'; gwin=gwin./sum(gwin);
axx=[-tbf:nbin/1000:taf]; axx=axx(1:end-1);
subplot(2,2,1)
ix= find(gain(:,1)<-.05 & gain(:,2)<-.05);
% m=conv2(nanmean(all_sync(ix,:),1),gwin,'same'); m=m(:)';
% s=conv2(nanstd(all_sync(ix,:),1),gwin,'same')./sqrt(length(ix)); s=s(:)';
sm=21;
m=smooth(nanmean(all_sync(ix,:),1),sm)-nanmean(nanmean(all_sync(ix,base_ix))); m=m(:)';
s=smooth(nanstd(all_sync(ix,:),[],1),sm)./sqrt(length(ix)); s=s(:)';
xx=[ax ax(end) ax(end) fliplr(ax) ax(1) ax(1)];
yy=[m-s m(end)-s(end) m(end)+s(end) fliplr(m+s) m(1)+s(1) m(1)-s(1)];
% patch(xx,yy,[0.5 0.5 1],'EdgeColor','none');
hold on
% m=conv2(nanmean(all_sync_random(ix,:),1),gwin,'same'); m=m(:)';
% s=conv2(nanstd(all_sync_random(ix,:),1),gwin,'same')./sqrt(length(ix)); s=s(:)';
m=smooth(nanmean(all_sync_random(ix,:),1),sm)-nanmean(nanmean(all_sync_random(ix,base_ix))); m=m(:)';
s=smooth(nanstd(all_sync_random(ix,:),[],1),sm)./sqrt(length(ix)); s=s(:)';
xx=[ax ax(end) ax(end) fliplr(ax) ax(1) ax(1)];
yy=[m-s m(end)-s(end) m(end)+s(end) fliplr(m+s) m(1)+s(1) m(1)-s(1)];
% patch(xx,yy,[0.5 0.5 0.5],'EdgeColor','none');
plot(ax,smooth(nanmean(all_sync(ix,:),1),sm)-nanmean(nanmean(all_sync(ix,base_ix))),'Color','b','LineWidth',1)
plot(ax,smooth(nanmean(all_sync_random(ix,:),1),sm)-nanmean(nanmean(all_sync_random(ix,base_ix))),'Color','k','LineWidth',1)
a=axis;
plot([0 0],[a(3) a(4)],'--k','LineWidth',2)
plot(a(1:2),[0 0],'--k','LineWidth',2)
hold off
set(gca,'TickDir','out'); box off
ylabel('SI'); xlabel('Time (s)')
title(['Pausers; n=' num2str(length(ix))])

subplot(2,2,2)
ix= find(gain(:,1)>.05 & gain(:,2)>.05);
% m=conv2(nanmean(all_sync(ix,:),1),gwin,'same'); m=m(:)';
% s=conv2(nanstd(all_sync(ix,:),1),gwin,'same')./sqrt(length(ix)); s=s(:)';
m=smooth(nanmean(all_sync(ix,:),1),sm)-nanmean(nanmean(all_sync(ix,base_ix))); m=m(:)';
s=smooth(nanstd(all_sync(ix,:),[],1),sm)./sqrt(length(ix)); s=s(:)';
xx=[ax ax(end) ax(end) fliplr(ax) ax(1) ax(1)];
yy=[m-s m(end)-s(end) m(end)+s(end) fliplr(m+s) m(1)+s(1) m(1)-s(1)];
% patch(xx,yy,[1 0.5 0.5],'EdgeColor','none');
hold on
% m=conv2(nanmean(all_sync_random(ix,:),1),gwin,'same'); m=m(:)';
% s=conv2(nanstd(all_sync_random(ix,:),1),gwin,'same')./sqrt(length(ix)); s=s(:)';
m=smooth(nanmean(all_sync_random(ix,:),1),sm)-nanmean(nanmean(all_sync_random(ix,base_ix))); m=m(:)';
s=smooth(nanstd(all_sync_random(ix,:),[],1),sm)./sqrt(length(ix)); s=s(:)';
xx=[ax ax(end) ax(end) fliplr(ax) ax(1) ax(1)];
yy=[m-s m(end)-s(end) m(end)+s(end) fliplr(m+s) m(1)+s(1) m(1)-s(1)];
% patch(xx,yy,[0.5 0.5 0.5],'EdgeColor','none');
plot(ax,smooth(nanmean(all_sync(ix,:),1),sm)-nanmean(nanmean(all_sync(ix,base_ix))),'Color','r','LineWidth',1)
plot(ax,smooth(nanmean(all_sync_random(ix,:),1),sm)-nanmean(nanmean(all_sync_random(ix,base_ix))),'Color','k','LineWidth',1)
a=axis;
plot([0 0],[a(3) a(4)],'--k','LineWidth',2)
plot(a(1:2),[0 0],'--k','LineWidth',2)
hold off
set(gca,'TickDir','out'); box off
ylabel('SI'); xlabel('Time (s)')
title(['Bursters; n=' num2str(length(ix))])

subplot(2,2,3)
ix= find((gain(:,1)>.05 & gain(:,2)<-.05) | (gain(:,1)<-.05 & gain(:,2)>.05));
% m=conv2(nanmean(all_sync(ix,:),1),gwin,'same'); m=m(:)';
% s=conv2(nanstd(all_sync(ix,:),1),gwin,'same')./sqrt(length(ix)); s=s(:)';
m=smooth(nanmean(all_sync(ix,:),1),sm)-nanmean(nanmean(all_sync(ix,base_ix))); m=m(:)';
s=smooth(nanstd(all_sync(ix,:),[],1),sm)./sqrt(length(ix)); s=s(:)';
xx=[ax ax(end) ax(end) fliplr(ax) ax(1) ax(1)];
yy=[m-s m(end)-s(end) m(end)+s(end) fliplr(m+s) m(1)+s(1) m(1)-s(1)];
% patch(xx,yy,[0.5 1 0.5],'EdgeColor','none');
hold on
% m=conv2(nanmean(all_sync_random(ix,:),1),gwin,'same'); m=m(:)';
% s=conv2(nanstd(all_sync_random(ix,:),1),gwin,'same')./sqrt(length(ix)); s=s(:)';
m=smooth(nanmean(all_sync_random(ix,:),1),sm)-nanmean(nanmean(all_sync_random(ix,base_ix))); m=m(:)';
s=smooth(nanstd(all_sync_random(ix,:),[],1),sm)./sqrt(length(ix)); s=s(:)';
xx=[ax ax(end) ax(end) fliplr(ax) ax(1) ax(1)];
yy=[m-s m(end)-s(end) m(end)+s(end) fliplr(m+s) m(1)+s(1) m(1)-s(1)];
% patch(xx,yy,[0.5 0.5 0.5],'EdgeColor','none');
plot(ax,smooth(nanmean(all_sync(ix,:),1),sm)-nanmean(nanmean(all_sync(ix,base_ix))),'Color','g','LineWidth',1)
plot(ax,smooth(nanmean(all_sync_random(ix,:),1),sm)-nanmean(nanmean(all_sync_random(ix,base_ix))),'Color','k','LineWidth',1)
a=axis;
plot([0 0],[a(3) a(4)],'--k','LineWidth',2)
plot(a(1:2),[0 0],'--k','LineWidth',2)
hold off
set(gca,'TickDir','out'); box off
ylabel('SI'); xlabel('Time (s)')
title(['Bursters-Pausers; n=' num2str(length(ix))])

subplot(2,6,10)
tix= find(ax>-.05 & ax<-.01);
base_all=nanmean(all_sync(:,tix),2);
tix= find(ax>.005 & ax<.045);
after_all=nanmean(all_sync(:,tix),2);
ix= find(gain(:,1)<-.05 & gain(:,2)<-.05);
M(1,:)= [nanmean(base_all(ix)) nanmean(after_all(ix))];
S(1,:)= [nanstd(base_all(ix)) nanstd(after_all(ix))]./sqrt(length(ix));
histogram(base_all(ix)-after_all(ix),'BinWidth',.05,'Normalization','count','FaceColor',[.5 .5 1])
hold on
a=axis;
plot([0 0],a(3:4),'--k')
plot([nanmean(base_all(ix)-after_all(ix)) nanmean(base_all(ix)-after_all(ix))],a(3:4),'-b')
hold off
[~,p(1)]= ttest(base_all(ix),after_all(ix));
title(p(1))
set(gca,'TickDir','out'); box off
xlabel('Base-AfterCS');
ylabel('#')

subplot(2,6,11)
ix= find(gain(:,1)>.05 & gain(:,2)>.05);
M(2,:)= [nanmean(base_all(ix)) nanmean(after_all(ix))];
S(2,:)= [nanstd(base_all(ix)) nanstd(after_all(ix))]./sqrt(length(ix));
histogram(base_all(ix)-after_all(ix),'BinWidth',.05,'Normalization','count','FaceColor',[1 .5 .5])
hold on
a=axis;
plot([0 0],a(3:4),'--k')
plot([nanmean(base_all(ix)-after_all(ix)) nanmean(base_all(ix)-after_all(ix))],a(3:4),'-r')
hold off
[~,p(2)]= ttest(base_all(ix),after_all(ix));
title(p(2))
set(gca,'TickDir','out'); box off
xlabel('Base-AfterCS');
ylabel('#')

subplot(2,6,12)
ix= find((gain(:,1)>.05 & gain(:,2)<-.05) | (gain(:,1)<-.05 & gain(:,2)>.05));
M(2,:)= [nanmean(base_all(ix)) nanmean(after_all(ix))];
S(2,:)= [nanstd(base_all(ix)) nanstd(after_all(ix))]./sqrt(length(ix));
histogram(base_all(ix)-after_all(ix),'BinWidth',.05,'Normalization','count','FaceColor',[1 .5 1])
hold on
a=axis;
plot([0 0],a(3:4),'--k')
plot([nanmean(base_all(ix)-after_all(ix)) nanmean(base_all(ix)-after_all(ix))],a(3:4),'-m')
hold off
[~,p(3)]= ttest(base_all(ix),after_all(ix));
title(p(3))
set(gca,'TickDir','out'); box off
xlabel('Base-AfterCS');
ylabel('#')
%%
figure
d=20;
axx=[-tbf:nbin/1000:taf]; axx=axx(1:end-1);
subplot(2,2,1)
ix= find(gain(:,1)<-.05 & gain(:,2)<-.05);
scatter(axx,nanmean(all_sync(ix,:),1),20,'MarkerFaceColor',[0.7 0.7 1],'MarkerEdgeColor','w');
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
scatter(axx,nanmean(all_sync(ix,:),1),20,'MarkerFaceColor',[1 0.7 .7],'MarkerEdgeColor','w');
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