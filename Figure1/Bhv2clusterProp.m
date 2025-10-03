% this is in response to Reviewer comment:
% there is a rough inverse % correlation across sessions b/w prevalence
% of late suppressing and late facilitating neurons. It would be good
% to test for systematic session-by-session differences in behavior that
% might explain this, to strengthen the suggested idea that it instead
% reflects spatial distribution.
% basically test whether behavior is different across the sessions with the
% different probability of different clusters
clear
all_xv_allTrials=[];
load FullData_withAutoClusters.mat
udir= 'E:/all_data/';
dates = dir([udir 'D*']);
mix=find(bhv_ax>=-1 & bhv_ax<=1);
new_bhv_ax=bhv_ax(mix);
d=(length(mix)-1)/2;
for i=1:length(AllData)
    for c=1:4
        n(i,c)=length(find(AllData(i).Trial(1).clusters==c))/length(AllData(i).Trial(1).clusters);
    end
    for ii=1:length(AllData(i).Trial)
        tmp=AllData(i).Trial(ii).bhv(1,:);
        [~,ix]=max(tmp(mix));
        ix= find(bhv_ax==new_bhv_ax(ix));
        ix=ix-d:ix+d;
        this_x(ii,:)=AllData(i).Trial(ii).bhv(1,ix);
        this_xv(ii,:)=AllData(i).Trial(ii).bhv(5,mix);
    end
    all_x(i,:)=nanmean(this_x,1);
    all_xv(i,:)=nanmean(this_xv,1);
    all_x_std(i,:)=nanstd(this_x,1)./sqrt(size(this_x,1));
    all_xv_std(i,:)=nanstd(this_xv,1)./sqrt(size(this_x,1));


    all_xv_allTrials=[all_xv_allTrials; this_xv];
    index=strfind(dates(i).name,'_');
    mouse_id{i}=dates(i).name(index(end-1)+1:index(end)+1);
    clear this_xv this_x;
end
%%
dax= [0:.2:1];
Cols=slanCM('winter',length(dax)-1);
for i=1:length(dax)-1
    for c=1:4
        index= find(n(:,c)>=dax(i) & n(:,c)<=dax(i+1));
        subplot(2,4,c)
        errorbar(new_bhv_ax,nanmean(all_x(index,:)),nanstd(all_x(index,:))./sqrt(length(index)),'Color',Cols(i,:),'CapSize',0)
        hold on
        xlim([-.5 .5])

        subplot(2,4,c+4)
        errorbar(new_bhv_ax,nanmean(all_xv(index,:)),nanstd(all_xv(index,:))./sqrt(length(index)),'Color',Cols(i,:),'CapSize',0)
        hold on
        xlim([-.5 .5])
    end
end
%%
index= find(n(:,1)>0.5 & n(:,4)<0.5);
subplot(2,2,1)
x=nanmean(all_x(index,:))-nanmean(nanmean(all_x(index,find(new_bhv_ax<-.1))));
errorbar(new_bhv_ax,x,nanstd(all_x(index,:))./sqrt(length(index)),'Color','b','CapSize',0);
index= find(n(:,1)<0.5 & n(:,4)>0.5);
hold on
x=nanmean(all_x(index,:))-nanmean(nanmean(all_x(index,find(new_bhv_ax<-.1))));
errorbar(new_bhv_ax,x,nanstd(all_x(index,:))./sqrt(length(index)),'Color','r','CapSize',0);
% xlim([-.5 .5])
a=axis;
plot([0 0],a(3:4),'--k')
hold off
set(gca,'TickDir','out'); box off
xlabel('time to endpoint (s)')
ylabel('Out position (cm)')

index= find(n(:,1)>0.4 & n(:,4)<0.2);
subplot(2,2,2)
x=nanmean(all_xv(index,:))-nanmean(nanmean(all_xv(index,find(new_bhv_ax<-.1))));
errorbar(new_bhv_ax,x,nanstd(all_xv(index,:))./sqrt(length(index)),'Color','b','CapSize',0);
index= find(n(:,1)>0.2 & n(:,4)<0.4);
hold on
x=nanmean(all_xv(index,:))-nanmean(nanmean(all_xv(index,find(new_bhv_ax<-.1))));
errorbar(new_bhv_ax,x,nanstd(all_xv(index,:))./sqrt(length(index)),'Color','r','CapSize',0);
% xlim([-.5 .5])
a=axis;
plot([0 0],a(3:4),'--k')
hold off
ylabel('Out velocity (cm/s)')
set(gca,'TickDir','out'); box off
xlabel('time to endpoint (s)')

%%
for j=1:length(dates)
    index=strfind(dates(j).name,'_');
    mouse_id{j}=dates(j).name(index(end-1)+1:index(end)+1);
end
[~,ix]=sort(n(:,1));
subplot(2,2,1)
barh(n(ix,:),'stacked')

subplot(2,2,2)
for i=1:length(ix)
    mid=find(strcmpi(mouse_id{ix(i)},unique(mouse_id)));
    scatter(1,i,50,mid,'fill')
    hold on
end
%%