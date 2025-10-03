clearvars -except AllData ax
load ../data/Figure3/FullData_withAutoClusters.mat
c_all=0;
xcorr_lim= [-.3 .3];
xcorr_bin=1;
xcorr_win=[-50 50];
xcorr_ax=xcorr_win(1):xcorr_bin:xcorr_win(2);
tbf=50; taf=50; %50 works
xcorr_ax=-tbf:taf;
iterations=500;
%%
% tix= find(ax>-.3 & ax<.3);
tix=1:length(ax);
for i=1:length(AllData)
    for c1=1:size(AllData(i).Trial(1).realNeu,1)-1
        for c2=c1+1:size(AllData(i).Trial(1).realNeu,1)
            this_ccg=[]; shuffle_ccg_peak=[];
            for ii=1:length(AllData(i).Trial)
                cell1= AllData(i).Trial(ii).realNeu(c1,tix);
                cell2= AllData(i).Trial(ii).realNeu(c2,tix);
                ncell2= AllData(i).Trial(ii).nullNeu(c2,tix);

                ix= find(cell1==1);
                ix= ix(find(ix>tbf & ix<length(cell1)-taf));
                INDEX= repmat(ix(:),1,length(xcorr_ax))+repmat(xcorr_ax,length(ix),1);

                raw=cell2(INDEX);
                null=ncell2(INDEX);
                this_ccg=[this_ccg; raw-null];

                % to calculate significane
                try
                    ix= find(cell1==0);
                    [~,rix]= sort(rand(1,length(ix)));
                    rix=ix(rix(1:length(find(cell1==1))));
                    INDEX= repmat(rix(:),1,length(xcorr_ax))+repmat(xcorr_ax,length(rix),1);
                    INDEX=INDEX(find(INDEX(:,1)>0 & INDEX(:,end)<length(cell1)),:);
                    raw=cell2(INDEX);
                    null=ncell2(INDEX);
                    shuffle_ccg_peak=[shuffle_ccg_peak; raw-null];
                end

            end
            c_all=c_all+1;
            all_ccg(c_all,:)=nanmean(this_ccg,1);
            all_clusters(c_all,:)=[AllData(i).Trial(1).clusters(c1) AllData(i).Trial(1).clusters(c2)];
            all_CS(c_all,:)=[AllData(i).Trial(1).CSv(c1) AllData(i).Trial(1).CSv(c2)];
            [~,all_p_value(c_all)]=ttest(shuffle_ccg_peak(:,find(xcorr_ax==0)),nanmean(this_ccg(:,find(xcorr_ax==0))));
            % all_p_value(c_all)=length(find(shuffle_ccg_peak(:,find(xcorr_ax==0))>all_ccg(c_all,find(xcorr_ax==0))))./size(shuffle_ccg_peak,1);
            [~,all_p_value_another(c_all)]=ttest(all_ccg(c_all,:),all_ccg(c_all,find(xcorr_ax==0)));
        end
    end
end
%%
orig_all_ccg=all_ccg;
orig_all_clusters=all_clusters;
all_ccg=orig_all_ccg;
all_clusters=orig_all_clusters;

%% calculate width
ix0=find(xcorr_ax==0);
V=zeros(size(all_ccg,1),1);
for i=1:size(all_ccg,1)
    tmp=smooth(all_ccg(i,:),5);
    %     tmp=all_ccg(i,:);
    m=tmp(ix0);
    if m<0
        %         all_W(i)=nan;
        tmp=tmp.*-1;
        %         continue
    end
    try
        ix1= ix0-find(flipud(tmp(1:ix0))<m/2,1,'first')+1;
        e1=interp1(tmp(ix1:ix0),xcorr_ax(ix1:ix0),m/2);

        ix2= ix0+find(tmp(ix0+1:end)<m/2,1,'first');
        e2=interp1(tmp(ix0+1:ix2),xcorr_ax(ix0+1:ix2),m/2);

        all_W(i)=e2-e1;
        if m>std(tmp)
            V(i)=1;
        end
    catch me
        all_W(i)=nan;
    end
end

%% calculate prominence
ix0=find(xcorr_ax==0);
V=zeros(size(all_ccg,1),1);
for i=1:size(all_ccg,1)
    tmp=smooth(all_ccg(i,:),5);
    %     tmp=all_ccg(i,:);
    m=tmp(ix0);
    m1=tmp(ix0-2);
    m2=tmp(ix0+2);
    if m<0
        %         all_W(i)=nan;
        tmp=tmp.*-1;
        %         continue
        m=tmp(ix0);
        m1=tmp(ix0-1);
        m2=tmp(ix0+1);
    end
    try
        m0=(abs(m1)+abs(m2))/2;

        all_G(i)=(m-m0)/(m+m0);
        if m>std(tmp)
            V(i)=1;
        end
    catch me
        all_G(i)=nan;
    end
end
%%
hplot=figure;
hhist=figure;
m=all_ccg(:,find(xcorr_ax==0));
ix=find(abs(m)<.2);
all_ccg=all_ccg(ix,:);
all_clusters=all_clusters(ix,:);
all_CS=all_CS(ix,:);
all_p_value=all_p_value(ix);
all_p_value_another=all_p_value_another(ix);
for i=1:4
    for ii=1:4
        ix=find((all_clusters(:,1)==i & all_clusters(:,2)==ii) | (all_clusters(:,1)==ii & all_clusters(:,2)==i)); % add all_CS==1 if you need to test only in CS identified PCs
        figure(hplot)
        subplot(4,4,(i-1)*4+ii)
        errorbar(xcorr_ax,nanmean(all_ccg(ix,:),1),nanstd(all_ccg(ix,:),[],1)./sqrt(length(ix)),'CapSize',0);
        hold on
        a=axis;
        plot([0 0],a(3:4),'--k')
        plot(a(1:2),[0 0],'--k')
        hold off
        xlabel('Time (ms)');
        ylabel('corrected CCG');
        title(['Cluster ' num2str(i) '-cluster ' num2str(ii)])
        xlim([-10 10])

        figure(hhist)
        subplot(4,4,(i-1)*4+ii)
        histogram(all_ccg(ix,find(xcorr_ax==0)),'BinWidth',5/1e4,'FaceColor','k','EdgeColor','none')
        hold on
        a=axis;
        plot([0 0],a(3:4),'--k')
        plot([mean(all_ccg(ix,find(xcorr_ax==0))) mean(all_ccg(ix,find(xcorr_ax==0)))],a(3:4),'b')
        hold off
        xlabel('corrected CCG')
        ylabel('#')
        [~,p]= ttest(all_ccg(ix,find(xcorr_ax==0)));
        title(p)
    end
end
%%
m=max(all_ccg,[],2);
h=isoutlier(m);
figure
for i=1:4
    ii=i;
    ix=find(((all_clusters(:,1)==i & all_clusters(:,2)==ii) | (all_clusters(:,1)==ii & all_clusters(:,2)==i)) & m<.2);
    subplot(3,4,i)
    errorbar(xcorr_ax,nanmean(all_ccg(ix,:),1),nanstd(all_ccg(ix,:),[],1)./sqrt(length(ix)),'CapSize',0);
    hold on
    a=axis;
    plot([0 0],a(3:4),'--k')
    plot(a(1:2),[0 0],'--k')
    hold off
    xlabel('Time (ms)');
    ylabel('corrected CCG');
    title(['Cluster ' num2str(i) '-cluster ' num2str(ii)])
    xlim([-10 10])
    set(gca,'TickDir','out'); box off

    subplot(3,4,4+i)
    histogram(all_ccg(ix,find(xcorr_ax==0)),'BinWidth',5/1e4,'FaceColor','k','EdgeColor','none')
    hold on
    a=axis;
    plot([0 0],a(3:4),'--k')
    plot([mean(all_ccg(ix,find(xcorr_ax==0))) mean(all_ccg(ix,find(xcorr_ax==0)))],a(3:4),'b')
    hold off
    xlabel('corrected CCG')
    ylabel('#')
    [~,p]= ttest(all_ccg(ix,find(xcorr_ax==0)));
    title(p)
    set(gca,'TickDir','out'); box off

    subplot(3,4,8+i)
    histogram(all_W(ix),'BinWidth',.1,'FaceColor','k','EdgeColor','none')
    hold on
    a=axis;
    plot([0 0],a(3:4),'--k')
    plot([nanmean(all_W(ix)) nanmean(all_W(ix))],a(3:4),'b')
    hold off
    xlabel('Width (ms)')
    ylabel('#')
    set(gca,'TickDir','out'); box off
    title(nansum(V(ix))/length(ix))
end
%%
for i=1:4
    p(i)=length(find(all_p_value(:)<.05 & all_clusters(:,1)==i & all_clusters(:,2)==i))/length(find(all_clusters(:,1)==i & all_clusters(:,2)==i));
end
p
%%
ccs='brgm';
D=[];
dax=linspace(min(all_ccg(:,find(xcorr_ax==0))),max(all_ccg(:,find(xcorr_ax==0))),51);
for i=1:4
    ii=i;
    ix=find(((all_clusters(:,1)==i & all_clusters(:,2)==ii) | (all_clusters(:,1)==ii & all_clusters(:,2)==i)) & m<.2);
    x=all_ccg(ix,find(xcorr_ax==0));
    for ii=1:length(dax)-1
        D(i,ii)= length(find(x>=dax(ii) & x<dax(ii+1)))/length(x);
    end
end

subplot(2,3,1)
imagesc(dax(1:end-1),1:4,D)
caxis([-.0 .05])
% set(gca,'YDir','normal')
colormap(flipud(hot))
set(gca,'TickDir','out'); box off
ylabel('Cluster #')
xlabel('cofiring probability at 1-ms')
colorbar
xlim([-.03 .1])

subplot(2,3,2)
for i=1:4
    ii=i;
    ix=find(((all_clusters(:,1)==i & all_clusters(:,2)==ii) | (all_clusters(:,1)==ii & all_clusters(:,2)==i)) & m<.2);
    x=all_ccg(ix,find(xcorr_ax==0));
    qqplot(x);
    hold on
    h=get(gca,'Children');
    h(1).Color=ccs(i);
    h(2).Color=ccs(i);
    h(3).Color=ccs(i);
end
hold off
set(gca,'TickDir','out'); box off


subplot(2,3,3)
for i=1:4
    ii=i;
    ix=find(((all_clusters(:,1)==i & all_clusters(:,2)==ii) | (all_clusters(:,1)==ii & all_clusters(:,2)==i)) & m<.2);
    x=all_ccg(ix,find(xcorr_ax==0));
    ecdf(x);
    hold on
        h=get(gca,'Children');
    h(1).Color=ccs(i);
end
hold off
legend
set(gca,'TickDir','out'); box off

subplot(2,3,4)
for i=1:4
    ii=i;
    ix=find(((all_clusters(:,1)==i & all_clusters(:,2)==ii) | (all_clusters(:,1)==ii & all_clusters(:,2)==i)) & m<.2);
    x=all_ccg(ix,find(xcorr_ax==0));
    probplot('normal',x);
    hold on
        h=get(gca,'Children');
    h(1).Color=ccs(i);
    h(2).Color=ccs(i);
end
hold off
set(gca,'TickDir','out'); box off


subplot(2,3,5)
qax=[0:0.1:1];
for i=1:4
    ii=i;
    ix=find(((all_clusters(:,1)==i & all_clusters(:,2)==ii) | (all_clusters(:,1)==ii & all_clusters(:,2)==i)) & m<.2);
    x=all_ccg(ix,find(xcorr_ax==0));
    quantile(x,qax);
    plot(qax,quantile(x,qax))
    hold on
end
hold off
set(gca,'TickDir','out'); box off

subplot(2,3,6)
for i=1:4
    ii=i;
    ix=find(((all_clusters(:,1)==i & all_clusters(:,2)==ii) | (all_clusters(:,1)==ii & all_clusters(:,2)==i)) & m<.2);
    x=all_ccg(ix,find(xcorr_ax==0));
    s=sign(x);
    absx=abs(x);
    logx=log10(absx);
    newx= logx.*s;
    % scatter(i+randn(1,length(x))./20,x,abs(x)*200,'fill')

    scatter(i+randn(1,length(x))./20,log((1+x)./(1-x))./2,'o','fill')
    hold on
end
hold off
ylabel('Z(ccg')

figure
c=0;
for i=1:4
    for ii=i:4
        c=c+1;
        subplot(4,4,c)
        ix=find(((all_clusters(:,1)==i & all_clusters(:,2)==i) | (all_clusters(:,1)==i & all_clusters(:,2)==i)) & m<.2);
        x=all_ccg(ix,find(xcorr_ax==0));
        ix=find(((all_clusters(:,1)==ii & all_clusters(:,2)==ii) | (all_clusters(:,1)==ii & all_clusters(:,2)==ii)) & m<.2);
        y=all_ccg(ix,find(xcorr_ax==0));
        qqplot(x,y,[0:1:100]);
        hold on
        h=get(gca,'Children');
        h(1).Color=ccs(i);
        h(2).Color=ccs(i);
        h(3).Color=ccs(i);
        title([num2str(i) '- ' num2str(ii)])
    end
end
hold off
set(gca,'TickDir','out'); box off