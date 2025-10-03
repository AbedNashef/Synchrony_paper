time0=-100/1000;
for i=1:4
    eval(['load([''E:/Dylans_data/jPSTHs/Cluster_jPETHs/normalized_jpsth_' num2str(i) '-' num2str(i) '.mat''])']);
    [ax,anti_diag] = get_anti_diag(all_jpsth_norm,[-.5:1/1000:.499],time0);
    [peakiness,vpeak]= get_peakiness(anti_diag,ax);
    [widths,vwidth]=get_width(anti_diag);
    subplot(4,4,(i-1)*4+1)
    real_ax=ax;
    ax=[-size(anti_diag,2)/2:size(anti_diag,2)/2]./1000;
    ax=ax(1:end-1);
    errorbar(ax,nanmean(anti_diag,1),nanstd(anti_diag,1)./sqrt(size(anti_diag,1)),'CapSize',0,'Color','k')

    subplot(4,4,(i-1)*4+2)
    histogram(peakiness,'BinWidth',.1,'FaceColor','k')
    title([length(find(vpeak<.05))/length(vpeak)])


    subplot(4,4,(i-1)*4+3)
    histogram(1-vpeak,'BinWidth',.1,'FaceColor','k')
    title([length(find(vpeak<.05))/length(vpeak)])

    subplot(4,4,(i-1)*4+4)
    histogram(widths,'BinWidth',2,'FaceColor','k','EdgeColor','none')
    hold on
    histogram(widths(find(vpeak<.05)),'BinWidth',2,'FaceColor','b','EdgeColor','none')
    hold  off

    all_anti_diag{i}=anti_diag;
    all_peakiness{i}=[peakiness(:) vpeak(:)];
    all_width{i}=[widths(:) vwidth(:)];

end
%%
gwin= gausswin(21);
gwin=gwin./sum(gwin);
colors=[0 0 1; 1 0 0; 0 1 0; 1 0 1];
clf
for i=1:length(all_anti_diag)
    z_anti_diag=[]; smoothed_anti_diag=[];
    anti_diag=all_anti_diag{i};
    peakiness= all_peakiness{i};
    width= all_width{i};
    for ii=1:size(anti_diag,1)
        smoothed_anti_diag(ii,:)=conv2(anti_diag(ii,:),gwin','same');
        z_anti_diag(ii,:)=conv2(anti_diag(ii,:),gwin','same');
        z_anti_diag(ii,:)=(z_anti_diag(ii,:)-nanmean(z_anti_diag(ii,:)))./nanstd(z_anti_diag(ii,:));
    end

    subplot(2,2,1)
    m=nanmean(smoothed_anti_diag,1);
    s=nanstd(smoothed_anti_diag,[],1)./sqrt(size(smoothed_anti_diag,1));
    patch([ax fliplr(ax)],[m-s fliplr(m+s)],colors(i,:),'edgeColor','none')
    hold on
    plot(ax,m,'Color',colors(i,:).*.5,'LineWidth',2)
    set(gca,'TickDir','out'); box off

    subplot(2,2,2)
    m=nanmean(z_anti_diag,1);
    s=nanstd(z_anti_diag,[],1)./sqrt(size(z_anti_diag,1));
    patch([ax fliplr(ax)],[m-s fliplr(m+s)],colors(i,:),'edgeColor','none')
    hold on
    plot(ax,m,'Color',colors(i,:).*.5,'LineWidth',2)
    set(gca,'TickDir','out'); box off

    subplot(2,3,4)
    bar(i,length(find(peakiness(:,2)<.05))./size(peakiness,1),'FaceColor',colors(i,:))
    p=length(find(peakiness(:,2)<.05))./size(peakiness,1);
    q=1-p;
    n=size(peakiness,1);
    se=sqrt((p*q)/n);
    hold on
    line([i i],[p-se p+se],'Color',colors(i,:))
    set(gca,'TickDir','out'); box off
    ylabel('fraction with significant peak')

    subplot(2,3,5)
    barwitherr(nanstd(peakiness(:,1))./sqrt(size(peakiness,1)),nanmean(peakiness(:,1)))
    hold on
    h=get(gca,'Children');
    h(1).XData=i;
    h(2).XData=i;
    h(2).FaceColor=colors(i,:);
    set(gca,'TickDir','out'); box off
    ylabel('Average Prominence')

    subplot(2,3,6)
    peaks=z_anti_diag(:,size(z_anti_diag,2)/2);
    barwitherr(nanstd(peaks)./sqrt(size(peakiness,1)),nanmean(peaks))
    hold on
    h=get(gca,'Children');
    h(1).XData=i;
    h(2).XData=i;
    h(2).FaceColor=colors(i,:);
    set(gca,'TickDir','out'); box off
    ylabel('Averga peak (z-score)')
    all_peaks{i}=peaks;
end