clear
addpath ../helper_functions/
n=40; %Purkinje cells
fr=40;%sp/s
isi=1000/fr;
t=1000;
iterations=1e3;
tau=2.5; 
lambda=1/tau;
time4decay=15;
decay_trace=exp(-lambda.*[1:time4decay]);
x=poissrnd(isi,n,fr);
T=cumsum(x,2);
MAT=zeros(n,max(max(T)));
ISI=[];
for i=1:n
    MAT(i,T(i,:))=1;
    ISI=[ISI diff(T(i,:))];
end
MAT=MAT(:,1:1000);
subplot(4,3,[1 4])
plot_raster(MAT,1:size(MAT,2))
set(gca,'TickDir','out'); box off
title('Simulated PC spikes');

subplot(4,3,7)
tmp = sum(MAT,1);
tmp(find(tmp>0))=1;
ix= find(tmp>0);
plot([ix; ix],[0 1],'k','LineWidth',.5)
set(gca,'TickDir','out'); box off
non_sync_fr=tmp;
title('Converging spikes');

Overall_spike_times=ix;

Nuclear_V= rand(1,t)./100;
for i=1:length(Overall_spike_times)
    if Overall_spike_times(i)+time4decay-1<t
        tmp= Nuclear_V(Overall_spike_times(i):Overall_spike_times(i)+time4decay-1);
        tmp= tmp+decay_trace;
%         tmp=tmp./max(tmp);
        Nuclear_V(Overall_spike_times(i):Overall_spike_times(i)+time4decay-1)=tmp;
    end
end
subplot(4,3,10)
plot(1:t,Nuclear_V,'k','LineWidth',2)
set(gca,'TickDir','out'); box off
% ylim([0 1.2])
ylabel('simulated Vmem')
xlabel('Time');

Non_sync_V=Nuclear_V;

%%
sync_fraction=[.05:.05:.6]; %50% synchrony
ISI=[];
for ss=1:length(sync_fraction)
    x=poissrnd(isi,n,fr);
    T=cumsum(x,2);
    MAT=zeros(n,max(max(T)));
    for i=1:n
        MAT(i,T(i,:))=1;
        ISI(i,ss)=nanmean(diff(T(i,:)));
    end
    MAT=MAT(:,1:1000);
    trace2sync= MAT(1,:);
    index= find(trace2sync);
    [~,rix] = sort(rand(1,length(index)));
    index2add=index(rix(1:round(length(rix)*sync_fraction(ss))));
    for i=2:size(MAT,1)
        trace2sync= MAT(i,:);
        index= find(trace2sync);
        [~,rix] = sort(rand(1,length(index)));
        index=index(rix(1:round(length(rix)*sync_fraction(ss))));
        MAT(i,index)=0;
        MAT(i,index2add)=1;
    end
    
    
    tmp = sum(MAT,1);
    tmp(find(tmp>0))=1;
    ix= find(tmp>0);
    sync_fr(ss,:)=tmp;
    
    Overall_spike_times=ix;
    
    Nuclear_V= rand(1,t)./100;
    for i=1:length(Overall_spike_times)
        if Overall_spike_times(i)+time4decay-1<t
            tmp= Nuclear_V(Overall_spike_times(i):Overall_spike_times(i)+time4decay-1);
            tmp= tmp+decay_trace;
%             tmp=tmp./max(tmp);
            Nuclear_V(Overall_spike_times(i):Overall_spike_times(i)+time4decay-1)=tmp;
        end
    end
    
    sync_V(ss,:)= Nuclear_V;
    sync_Mat(ss,:,:)=MAT;
end

index= find(round(sync_fraction.*100)==30);

subplot(4,3,[2 5])
plot_raster(squeeze(sync_Mat(index,:,:)),1:size(MAT,2))
set(gca,'TickDir','out'); box off
title('Simulated PC spikes with 30% sync');

subplot(4,3,8)
tmp =sum(squeeze(sync_Mat(index,:,:)),1);
tmp(find(tmp>0))=1;
ix= find(tmp>0);
plot([ix; ix],[0 1],'k','LineWidth',.5)
set(gca,'TickDir','out'); box off
title('Converging spikes');

subplot(4,3,11)
plot(1:t,sync_V(index,:),'k','LineWidth',2)
set(gca,'TickDir','out'); box off
% ylim([0 1.2])
ylabel('simulated Vmem')
xlabel('Time');

%% iterations
all_fr=[20:2:120];
MAT=[];
Nuclear_V=[];
all_sync_fr=[];
ISI=[];
sync_fraction=[0:.01:.6]; %50% synchrony
iterations=1000;
for ff=1:length(all_fr)
    for ss=1:length(sync_fraction)
        for iter=1:iterations
            isi=1000/all_fr(ff);
            x=poissrnd(isi,n,all_fr(ff));
            x(x==0)=1;
            T=cumsum(x,2);
            tmp=zeros(n,max(max(T)));
            for i=1:n
                tmp(i,T(i,:))=1;
            end
            XX=tmp(:,1:1000);
            x=XX(:);
            [~,ix]=sort(rand(1,length(x))); x=x(ix);
            XX=reshape(x,n,[]);
%             XX=desynch_mat(XX);
            MAT(iter,:,:)=XX;
            trace2sync= MAT(iter,1,:);
            index= find(trace2sync);
            [~,rix] = sort(rand(1,length(index)));
            index2add=index(rix(1:round(length(rix)*sync_fraction(ss))));
            for i=2:size(MAT,2)
                trace2sync= MAT(iter,i,:);
                index= find(trace2sync);
                [~,rix] = sort(rand(1,length(index)));
                index=index(rix(1:round(length(rix)*sync_fraction(ss))));
                MAT(iter,i,index)=0;
                MAT(iter,i,index2add)=1;
            end
            
            for i=1:size(MAT,2)
                tmp= diff(find(squeeze(MAT(iter,i,:))));
            end
%             ISI(iter,ff,ss)=nanmean(tmp);
            tmp = sum(squeeze(MAT(iter,:,:)),1);
            tmp(tmp>0)=1;
            ix= find(tmp>0);
            Overall_spike_times=ix;
            all_sync_fr(ff,ss,iter,:)=tmp;
            ISI(iter,ff,ss)=nanmean(diff(find(tmp)));
            Nuclear_V(ff,ss,iter,:)= rand(1,t)./100;
            for i=1:length(Overall_spike_times)
                if Overall_spike_times(i)+time4decay-1<t
                    tmp= squeeze(Nuclear_V(ff,ss,iter,Overall_spike_times(i):Overall_spike_times(i)+time4decay-1));
                    tmp= tmp(:)+decay_trace(:);
                    %         tmp=tmp./max(tmp);
                    Nuclear_V(ff,ss,iter,Overall_spike_times(i):Overall_spike_times(i)+time4decay-1)=tmp;
                end
            end
        end
    end
end
%%
addpath ../synchrony/
COLS=flipud(brewermap(length(all_fr),'reds'));
sync_fr2plot=squeeze(nansum(all_sync_fr,4));
subplot(4,3,3)
pcolor(sync_fraction,all_fr,squeeze(nanmean(sync_fr2plot,3)))
shading flat
title('Converging spikes/s')
xlabel('%Synchrony')
ylabel('FR')
set(gca,'TickDir','out'); box off
colormap jet 

subplot(4,3,6)
tmp=[];
for i=1:length(all_fr)
%     errorbar([0:10:100],(((nanmean(sync_fr2plot(i,:,:)-sync_fr2plot(i,1,:),3))./n)./all_fr(i))*100,((((nanstd(sync_fr2plot(i,:,:)-sync_fr2plot(i,1,:),[],3))./n)./all_fr(i))*100)./sqrt(iterations),'Color',COLS(i,:));%[nansum(non_sync_fr); nansum(sync_fr,2)],'-ok');
    tmp(i,:)=(((nanmean(sync_fr2plot(i,:,:)-sync_fr2plot(i,1,:),3))./n)./all_fr(i))*100;
end
pcolor(sync_fraction,all_fr,tmp)

% errorbar([0:10:100],(((nanmean(sync_fr2plot-sync_fr2plot(1,:),2))./n)./fr)*100,((((nanstd(sync_fr2plot-sync_fr2plot(1,:),[],2))./n)./fr)*100)./sqrt(iterations),'Color','k');%[nansum(non_sync_fr); nansum(sync_fr,2)],'-ok');
% plot([0 sync_fraction.*100],(([[nansum(non_sync_fr); nansum(sync_fr,2)]-nansum(non_sync_fr)]./n)./fr).*100,'-ok');
title('\Delta % spikes/s')
xlabel('%Synchrony')
ylabel('FR')
set(gca,'TickDir','out'); box off
hold off
shading flat
cla

subplot(4,3,6)
tmp=[];
for i=1:length(all_fr)
%     errorbar([0:10:100],(((nanmean(sync_fr2plot(i,:,:)-sync_fr2plot(i,1,:),3))./n)./all_fr(i))*100,((((nanstd(sync_fr2plot(i,:,:)-sync_fr2plot(i,1,:),[],3))./n)./all_fr(i))*100)./sqrt(iterations),'Color',COLS(i,:));%[nansum(non_sync_fr); nansum(sync_fr,2)],'-ok');
    tmp(i,:)=squeeze(nanmean(ISI(:,i,:),1))./squeeze(nanmean(ISI(:,i,1),1));
end
pcolor(sync_fraction,all_fr,tmp)
title('avg. ISI (ms)')
ylabel('FR')
xlabel('%Synchrony')
set(gca,'TickDir','out'); box off
shading flat


subplot(4,3,9)
tmp=[];
for i=1:length(all_fr)
%     errorbar([0:10:100],(((nanmean(sync_fr2plot(i,:,:)-sync_fr2plot(i,1,:),3))./n)./all_fr(i))*100,((((nanstd(sync_fr2plot(i,:,:)-sync_fr2plot(i,1,:),[],3))./n)./all_fr(i))*100)./sqrt(iterations),'Color',COLS(i,:));%[nansum(non_sync_fr); nansum(sync_fr,2)],'-ok');
    tmp(i,:)=squeeze(nanstd(ISI(:,i,:),[],1));
end
pcolor(sync_fraction,all_fr,tmp)
title('st.d. ISI (ms)')
ylabel('FR')
xlabel('%Synchrony')
set(gca,'TickDir','out'); box off
shading flat


subplot(4,3,12)
tmp=[];
tmp2=[];
for i=1:length(all_fr)
    tmp(i,:)=nanmean(nansum(Nuclear_V(i,:,:,:),4),3);
    tmp2(i,:)=nanstd(nansum(Nuclear_V(i,:,:,:),4),[],3);
    % errorbar(sync_fraction.*100,nanmean(nanmean(Nuclear_V,2),3),nanstd(nanmean(Nuclear_V,2),[],3)./sqrt(iterations),'Color','k')
    % plot([0:10:100],[nansum(Non_sync_V); nansum(sync_V,2)],'-ok');
end
pcolor(sync_fraction,all_fr,tmp)
set(gca,'TickDir','out'); box off
title('\SigmaIPSC/s')
ylabel('FR')
xlabel('%Synchrony')
shading flat
% % % 
% % % subplot(4,3,12)
% % % addpath ../synchrony/
% % % C=flipud(brewermap(length(sync_fraction),'YlGnBu'));
% % % ax=1:t;
% % % X=[ax ax(end) ax(end) flipud(ax) ax(1) ax(1)];
% % % ix= find(all_fr==40);
% % % for i=1:size(Nuclear_V,2)
% % %     y=smooth(squeeze(nanmean(Nuclear_V(ix,i,:,:),3)),21)';
% % %     s=smooth(squeeze(nanstd(Nuclear_V(ix,i,:,:),[],3)),21)';
% % %     YY=[y-s y(end)-s(end) y(end)+s(end) fliplr(y+s) y(1)+s(1) y(1)-s(1)];
% % %     %     patch(X,YY,C(i,:),'EdgeColor','none')
% % %     plot(ax,y,'Color',C(i,:),'LineWidth',2)
% % % %     plot(1:t,smooth(squeeze(nanmean(Nuclear_V(i,:,:),2)),21),'Color',C(i,:),'LineWidth',2)
% % %     hold on
% % % end
% % % hold off
% % % set(gca,'TickDir','out'); box off
% % % xlabel('Time');
% % % ylabel('Simulated Vmem')