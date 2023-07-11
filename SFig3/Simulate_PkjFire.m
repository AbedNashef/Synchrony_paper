% create simulated data
addpath ../helper_functions/
clear
iterations=1000;
t=2000;
taxis= linspace(-1,1,t);
bin=10;
ntr=60;
all_cells_mat=zeros(iterations,ntr,t);
for iter=1:iterations
    MO=1000;
    downgrad=60;
    upgrad=500;
    freq=randi([50 110]);
    mo_frq=randi([0 freq]);
    avg_isi= 1000/freq;
    mavg_isi= 1000/mo_frq;
    downward= linspace(freq,mo_frq,downgrad/bin);
    upward= linspace(mo_frq,freq,upgrad/bin);
    TT=0:bin:t;
    Prob=[ones(1,(t-MO-downgrad/2)/bin).*freq downward upward];
    Prob=smooth([Prob ones(1,t/bin+1-length(Prob)).*freq],51);
    for i=1:length(TT)
        x(:,i)=poissrnd(1000/Prob(i),[ntr,1]);
    end
    ix= 1:size(x,2);%round(linspace(1,size(x,2),freq));
    x=x(:,ix);
    index= cumsum(x,2);
    index(find(index>t))=nan;
    for i=1:ntr
        ix= find(~isnan(index(i,:)));
        ixx=index(i,ix);
        ixx=ixx(find(ixx>0));
        all_cells_mat(iter,i,ixx)=1;
    end
end
%%
orig_mat=all_cells_mat;
orig_ax=taxis;
%%
all_cells_mat=orig_mat;
taxis=orig_ax;
clear means_trc;
%%
ix= find(taxis>-.9 & taxis<.9);
taxis=taxis(ix);
all_cells_mat=all_cells_mat(:,:,ix);
for i=1:size(all_cells_mat,1)
    means_trc(i,:)= squeeze(nanmean(all_cells_mat(i,:,:),2));
end
subplot(2,2,1)
plot(taxis,nanmean(means_trc,1),'b')
a=axis;
hold on
plot([0 0],a(3:4),'--k')
hold off
set(gca,'TickDir','out'); box off
xlabel('Time')
ylabel('Pkj. Firing porbability');
%%
npercell=5;
c_all=0;
for i=1:size(all_cells_mat,1)
    cell1=squeeze(all_cells_mat(i,:,:));
    [~,ix]= sort(rand(1,iter)); ix=ix(find(ix~=i));
    index2use=ix(1:npercell);
    for ii=1:npercell
        c_all=c_all+1;
        
        cell2=squeeze(all_cells_mat(index2use(ii),:,:));
        
        Pr_S1= nanmean(cell1,1);
        Pr_S2= nanmean(cell2,1);
        Pr_S1S2= nanmean(cell1.*cell2,1);
        all_sync(c_all,:)= Pr_S1S2./(Pr_S1.*Pr_S2);
        
       
        [raw, shift_predict, pred, surprise, std, ~, ~] = my_JPSTH(cell1', cell2', 10,0.5);
        all_jpsth(c_all,:,:)=  (raw-pred)./std;
    end
end
%%
orig_sync= all_sync;
%%
sm=21;
for i=1:size(all_sync,1)
    all_sync(i,:)=smooth(all_sync(i,:),sm);
end
%%
sm=1;
m=smooth(nanmean(all_sync,1),sm); m=m(:)'; 
s=smooth(nanstd(all_sync,1),sm)./sqrt(size(all_sync,1)); s=s(:)';
xx=[taxis taxis(end) taxis(end) fliplr(taxis) taxis(1) taxis(1)];
yy=[m-s m(end)-s(end) m(end)+s(end) fliplr(m+s) m(1)+s(1) m(1)-s(1)];
subplot 222
patch(xx,yy,[0.5 0.5 1],'EdgeColor','none');
hold on
plot(taxis,nanmean(all_sync,1),'Color','b','LineWidth',2)
a=axis;
plot([0 0],a(3:4),'--k');
plot(a(1:2),[1 1],'--k');
hold off
set(gca,'TickDir','out'); box off
xlabel('Time');
ylabel('SI')
%%
ax= taxis(1):.010:taxis(end);
subplot 223
pcolor(ax,ax,squeeze(nanmean(all_jpsth,1)));
hold on
plot([ax(1) ax(end)],[0 0],'Color','w','LineWidth',2)
plot([0 0],[ax(1) ax(end)],'Color','w','LineWidth',2)
plot([ax(1) ax(end)],[ax(1) ax(end)],'Color','w','LineWidth',2)
hold off
shading flat
colormap jet
set(gca,'TickDir','out'); box off
xlabel('Time'); ylabel('Time');

for i=1:size(all_jpsth,1)
    d(i,:)=diag(squeeze(all_jpsth(i,:,:)));
end
subplot(2,2,4)
m=smooth(nanmean(d,1),sm); m=m(:)'; 
s=smooth(nanstd(d,1),sm)./sqrt(size(d,1)); s=s(:)';
xx=[ax ax(end) ax(end) fliplr(ax) ax(1) ax(1)];
yy=[m-s m(end)-s(end) m(end)+s(end) fliplr(m+s) m(1)+s(1) m(1)-s(1)];
subplot 224
patch(xx,yy,[0.5 0.5 1],'EdgeColor','none');
hold on
plot(ax,nanmean(d,1),'Color','b','LineWidth',2)
a=axis;
plot([0 0],a(3:4),'--k');
plot(a(1:2),[0 0],'--k');
hold off
set(gca,'TickDir','out'); box off
xlabel('Time'); ylabel('Avg. CoFiring (sp^2/s^2)');
%%
all_sync=orig_sync;
%% Synchronize some data
Sync=[.25 .5 .75 1];
diff_synchronized(1,:,:)=all_sync;
diff_jPSTH(1,:,:,:)=all_jpsth;
index= find(taxis>0 & taxis<.5);
for ss=1:length(Sync)
    npercell=5;
    c_all=0;
    for i=1:size(all_cells_mat,1)
        cell1=squeeze(all_cells_mat(i,:,:));
        [~,ix]= sort(rand(1,iter)); ix=ix(find(ix~=i));
        index2use=ix(1:npercell);
        for ii=1:npercell
            c_all=c_all+1;
            
            Temp= cell1(:,index);
            six= find(Temp); l= length(six);
            [~,new_ix]= sort(rand(1,length(six)));
            six=six(new_ix);
            Temp(six(1:round((1-Sync(ss))*l)))=0;
            cell2=squeeze(all_cells_mat(index2use(ii),:,:));
            cell2(:,index)=cell2(:,index)+Temp; 
            cell2(find(cell2>1))=1;
            
            Pr_S1= nanmean(cell1,1);
            Pr_S2= nanmean(cell2,1);
            Pr_S1S2= nanmean(cell1.*cell2,1);
            diff_synchronized(1+ss,c_all,:)= Pr_S1S2./(Pr_S1.*Pr_S2);
            
            
            [raw, shift_predict, pred, surprise, std, ~, ~] = my_JPSTH(cell1', cell2', 10,0.5);
            diff_jPSTH(1+ss,c_all,:,:)=  (raw-pred)./std;
        end
    end
end

%%
orig_all_sync=all_sync;
%%
C=colormap(jet);
Cols=C(round(linspace(1,size(C,1),length(Sync)+1)),:);
figure
sm=1;
for i=1:size(diff_synchronized,1)
    all_sync=squeeze(diff_synchronized(i,:,:));
    m=smooth(nanmean(all_sync,1),sm); m=m(:)';
    s=smooth(nanstd(all_sync,1),sm)./sqrt(size(all_sync,1)); s=s(:)';
    xx=[taxis taxis(end) taxis(end) fliplr(taxis) taxis(1) taxis(1)];
    yy=[m-s m(end)-s(end) m(end)+s(end) fliplr(m+s) m(1)+s(1) m(1)-s(1)];
    subplot 221
%     patch(xx,yy,[0.5 0.5 1],'EdgeColor','none');
    hold on
    plot(taxis,nanmean(all_sync,1),'LineWidth',2,'Color',Cols(i,:))
    a=axis;
    plot([0 0],a(3:4),'--k');
    plot(a(1:2),[1 1],'--k');
    set(gca,'TickDir','out'); box off
    xlabel('Time');
    ylabel('SI')
    
    subplot 222
    hold on
    bar(i,nanmean(nanmean(all_sync(:,index))),'FaceColor',Cols(i,:))
    errorbar(i,nanmean(nanmean(all_sync(:,index))),nanstd(nanmean(all_sync(:,index),2))./sqrt(size(all_sync,1)),'Color','k')
end
subplot 221; hold off
subplot 222; hold off
ylabel('Avg. SI (0-500ms)')
set(gca,'TickDir','out'); box off
set(gca,'XTick',1:length(Sync)+1,'XTickLabel',[0 Sync])
%%
orig_jpsth=all_jpsth;
%%
ax= taxis(1):.010:taxis(end);
for i=1:size(diff_jPSTH,1)
    subplot(2,size(diff_jPSTH,1),i+size(diff_jPSTH,1))
    all_jpsth=squeeze(diff_jPSTH(i,:,:,:));
    pcolor(ax,ax,squeeze(nanmean(all_jpsth,1)));
    hold on
    plot([ax(1) ax(end)],[0 0],'Color','w','LineWidth',2)
    plot([0 0],[ax(1) ax(end)],'Color','w','LineWidth',2)
%     plot([ax(1) ax(end)],[ax(1) ax(end)],'Color','w','LineWidth',2)
    hold off
    shading flat
    colormap jet
    set(gca,'TickDir','out'); box off
    xlabel('Time'); ylabel('Time');
end