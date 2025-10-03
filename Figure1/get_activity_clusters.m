%% cross-correlation of cells + hierarchical grouping to cluster...
% this is based on Kat's program
function clusters= get_activity_clusters(PsthMaster,ax,times2take,no_clusters)

%xcorr_traces = corrcoef(cell2mat(cellfun(@(x) nanmean(x,1), snips_cat, 'uni', 0))');
PsthMasterClip=PsthMaster(:,times2take); %Clip PSTH so that you're only clustering around the time of reach

%%
%xcorr_traces = corrcoef(PsthMasterClip');
xcorr_traces = corrcoef(PsthMaster');

% set auto-correlations to zero
for ii=1:size(xcorr_traces,1)
    xcorr_traces(ii,ii)=0;
end
Z = linkage(xcorr_traces,'complete','correlation');
clusters = cluster(Z,'Maxclust',no_clusters);

hFig = figure('color','white');
subplot(2,3,1)
% dendrogram_ax = subplot(1, 3, 1);
%dendrogram_ax = axes('Position', [0.05, 0.1, 0.2, 0.8]);

% do cutoff for 'no_clusters' clusters
cutoff = median([Z(end-(no_clusters-1),3) Z(end-(no_clusters-2),3)]);
[H,T,outperm]=dendrogram(Z,0,'Orientation','left','ColorThreshold',cutoff);
set(H,'LineWidth',2)
% hFig=figure(100);
% set(hFig, 'units','normalized','position',[.05 .1 .14 .8])



% figure('color','white');
subplot(2,3,2)
hold on;

%heatmap1_ax = subplot(1, 2, 1);
%heatmap_ax = axes('Position', [0.3, 0.1, 0.6, 0.8]);
% re-sort matrix
for ii=1:size(xcorr_traces,1)
    for jj=1:size(xcorr_traces,1)
        perm_xcorr_traces(ii,jj)=xcorr_traces(outperm(ii),outperm(jj));
    end
end

pcolor(perm_xcorr_traces)

%colormap fire
colormap(brewermap([],('PRGn')))
shading flat
axis tight



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%dFF=cell2mat(cellfun(@(x) nanmean(x,1), snips_cat, 'uni', 0));
%dFF=PsthMasterClip;
dFF=PsthMaster;
% z-score
dFF=(dFF - nanmean(dFF(:,1:200),'all')) ./ std(dFF(:,1:200),[],"all",'omitnan');

% figure('color','white');
subplot(2,3,3)

hold on;
%heatmap_ax2 = axes('Position', [2.5, 0.1, 1.8, 0.8]);
%heatmap2_ax = subplot(1, 3, 2);

sorted_dFF_traces=[];
for ii=1:size(xcorr_traces,1)
    sorted_dFF_traces(ii,:)=dFF(outperm(ii),:);
end

% windowSize=1971;
% time_span = (-windowSize/2:windowSize/2)';
time_span=ax';%find(ax>-.5 & ax<.5)';
%time_span=[-200: 200]';
time_span_mat=repmat(time_span',ii,1);
ROI_mat=repmat(1:ii,length(time_span),1)';

pcolor(time_span_mat,ROI_mat,sorted_dFF_traces)
% pcolor(time_span_mat,ROI_mat,dFF)
%colormap fire
colormap(brewermap([],('*RdPu')))
shading flat
axis tight

line([0 0],get(gca,'ylim'),'color',[0,0,0],'linewidth',0.5,'linestyle','--');
end