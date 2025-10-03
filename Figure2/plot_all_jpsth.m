function plot_all_jpsth
h_jpsth=figure;
h_diags=figure;
for i=1:3
    eval(['load([''E:\Dylans_data\jPSTHs\Cluster_jPETHs/normalized_by_tr_jpsth_' num2str(i) '-' num2str(i) '.mat''])']);
    plot_jpsth_data(all_jpsth_norm,[-.5:1/1000:.499],i,i,4,4, h_diags,h_jpsth)
    for ii=i+1:4
        eval(['load([''E:\Dylans_data\jPSTHs\Cluster_jPETHs/normalized_by_tr_jpsth_' num2str(i) '-' num2str(ii) '.mat''])']);
        plot_jpsth_data(all_jpsth_norm,[-.5:1/1000:.499],i,ii,4,4, h_diags,h_jpsth)
    end
end
i=4;
eval(['load([''E:\Dylans_data\jPSTHs\Cluster_jPETHs/normalized_by_tr_jpsth_' num2str(i) '-' num2str(i) '.mat''])']);
plot_jpsth_data(all_jpsth_norm,[-.5:1/1000:.499],i,i,4,4, h_diags,h_jpsth)
end

function plot_jpsth_data(mat,ax,i1,i2,r,c,h_diags,h_jpsth)
for i=1:size(mat,1)
    mat(i,:,:)=imgaussfilt(squeeze(mat(i,:,:)),2);
    all_diag(i,:)=diag(squeeze(mat(i,:,:)));
end
figure(h_jpsth)
subplot(r,c,(i1-1)*c+i2)
pcolor(ax,ax,squeeze(nanmean(mat,1)))
shading flat
set(gca,'TickDir','out'); box off
xlabel('Time (s)')
ylabel('Time (s)')

figure(h_diags)
subplot(r,c,(i1-1)*c+i2)
m=nanmean(all_diag,1);
s=nanstd(all_diag,1)./sqrt(size(all_diag,1));
patch([ax fliplr(ax)],[m-s fliplr(m+s)],[.5 .5 .5],'EdgeColor','none');
hold on
plot(ax,m,'k','LineWidth',2)

m=nanmean(all_diag,1)-nanmean(nanmean(all_diag(:,find(ax<-.3))));
patch([ax fliplr(ax)],[m-s fliplr(m+s)],[.5 .5 .5],'EdgeColor','none');
plot(ax,m,'k','LineWidth',2)
hold off

xlabel('Time (s)')
ylabel('R^2');
end