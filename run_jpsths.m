% script to calculate and plot panels in Fig. 1
clear
udir= '../new_data/';
dates = dir([udir 'D*']);
bin2use=1;
addpath ../helper_functions/
tbf=.5; taf=.75;
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
thresh=.2;
%%
nax= [-0.5:1/1000:0.5];nax=nax(1:end-1);
t1= find(nax<-0.25); t2= find(nax>0 & nax<.25);

% [ res shift_predict psth_pred surprise_mat std_mat NUM_TRIAL pval_mat] = my_JPSTH(cut1, cut2, BIN_SIZE,gauss_filt_std)
%%
for j=1:length(dates)
    load([udir dates(j).name]);
    c_trials=0;
    for i=1:length(ReachS)
        times = ReachS(i).filt_kin(:,1);
        ix= find(ReachS(i).out(:,2)>thresh,1,'first');
        endpoint_time=ReachS(i).out(ix,1);
        endpoint_time=ReachS(i).out(end,1);
        if isempty(endpoint_time)
            ix=find(ReachS(i).filt_kin(:,2)>thresh,1,'first');
            endpoint_time=ReachS(i).filt_kin(ix,1);
        end
        times= ReachS(i).filt_kin(:,1);
        P1= repmat(ReachS(i).filt_kin(1,2:4),size(ReachS(i).filt_kin,1),1);
        tt=times(find(times>=endpoint_time-tbf & times<=endpoint_time+taf));
        
        vStimMode = isfield(ReachS(i),'stim');
        if vStimMode, vStim=ReachS(i).stim; else, vStim=0; end
        if isempty(vStim), vStim=0; end
        vExclude=ReachS(i).exclude;
        c_this=0;
        if ~vExclude & ~vStim
            c_trials=c_trials+1;
            for cc=1:length(cellData)
                c_this=c_this+1;
                this_trc=cellData(cc).Bin1;
                index= find(this_trc(:,1)>=endpoint_time-tbf & this_trc(:,1)<=endpoint_time+taf);
                FR(c_this,c_trials,:)= this_trc(index,2);
                
                all_gain(c_this)= cellData(cc).gain;
                
                Chs{c_this}= cellData(cc).Channels;
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
            if overlap==0 & ((all_gain(i)<-.05 & all_gain(ii)>.05) | (all_gain(ii)<-.05 & all_gain(i)>.05))
                c_all=c_all+1;
                cell1= squeeze(FR(i,:,:));
                cell2= squeeze(FR(ii,:,:));
                gain(c_all,1)=all_gain(i);
                %             Modulation(c_all,1) =abs(fr2-fr1)*1000;
                gain(c_all,2)=all_gain(ii);
                %             Modulation(c_all,2) =abs(fr2-fr1)*1000;
                if gain(c_all,1)>0 & gain(c_all,2)<0
                    gain(c_all,:)=fliplr(gain(c_all,:));
                    %                 Modulation(c_all,:)=fliplr(Modulation(c_all,:));
                    new_cell1=cell2;
                    new_cell2=cell1;
                    cell1=new_cell1; cell2=new_cell2;
                end
                [raw, shift_predict, pred, surprise, std, ~, ~] = my_JPSTH(cell1', cell2', nbin,0.5); 
                % You can find my_JPSTH in https://zenodo.org/record/6464052
                
                % %             all_jpsth(c_all,:,:)=raw-pred;
                num_trials= size(cell1,1);
                all_jpsth_norm(c_all,:,:)=(raw-pred)./num_trials;%./std;
                date_id(c_all)= j;
                Dist4cells(c_all)= abs(cellData(i).depth-cellData(ii).depth);
            end
        end
    end
    
    clear FR FR2gain this_S nX nY nZ
end
%%
sm=3;
axx=[-tbf:nbin/1000:taf]; axx=axx(1:end-1);
subplot(3,3,1)
ix= find(gain(:,1)<-.05 & gain(:,2)<-.05);
pcolor(axx,axx,imgaussfilt(squeeze(nanmean(all_jpsth_norm(ix,:,:),1)),sm))
shading flat
hold on
plot([axx(1) axx(end)],[0 0],'-w','LineWidth',2)
plot([0 0],[axx(1) axx(end)],'-w','LineWidth',2)
plot([axx(1) axx(end)],[axx(1) axx(end)],'-w','LineWidth',3)
hold off
title(['Pausers; n=' num2str(length(ix))]);

subplot(3,3,2)
ix= find(gain(:,1)>.05 & gain(:,2)>.05);
pcolor(axx,axx,imgaussfilt(squeeze(nanmean(all_jpsth_norm(ix,:,:),1)),sm))
shading flat
hold on
plot([axx(1) axx(end)],[0 0],'-w','LineWidth',2)
plot([0 0],[axx(1) axx(end)],'-w','LineWidth',2)
plot([axx(1) axx(end)],[axx(1) axx(end)],'-w','LineWidth',2)
hold off
title(['Bursters; n=' num2str(length(ix))]);


subplot(3,3,3)
ix= find(gain(:,1)<-.05 & gain(:,2)>.05);
pcolor(axx,axx,imgaussfilt(squeeze(nanmean(all_jpsth_norm(ix,:,:),1)),sm))
shading flat
hold on
plot([axx(1) axx(end)],[0 0],'-w','LineWidth',2)
plot([0 0],[axx(1) axx(end)],'-w','LineWidth',2)
plot([axx(1) axx(end)],[axx(1) axx(end)],'-w','LineWidth',2)
hold off
title(['Pausers-Bursters; n=' num2str(length(ix))]);


diag_win=5; %bins
k2use=[-diag_win:diag_win];
for k=1:size(all_jpsth_norm,1)
    X=imgaussfilt(squeeze(all_jpsth_norm(k,:,:)),sm);
    d=nan(length(k2use),length(axx));
    for i=1:length(k2use)
        if k2use(i)<0
            tmp=diag(X,k2use(i));
            d(i,1:end+k2use(i))=tmp;
        else
            d(i,1+k2use(i):end)=diag(X,k2use(i));
        end
    end
%     all_diag(k,:)=nanmean(d,1);
    all_diag(k,:)=d(find(k2use==0),:);
end

subplot(3,3,4)
ix= find(gain(:,1)<-.05 & gain(:,2)<-.05);
plot(axx,nanmean(all_diag(ix,:),1),'b');
shading flat
hold on
plot(axx,nanmean(all_diag(ix,:),1)+nanstd(all_diag(ix,:),[],1)./sqrt(length(ix)),'b');
plot(axx,nanmean(all_diag(ix,:),1)-nanstd(all_diag(ix,:),[],1)./sqrt(length(ix)),'b');
plot([axx(1) axx(end)],[0 0],'--k','LineWidth',1)
a=axis;
plot([0 0],[a(3) a(4)],'--k','LineWidth',2)
hold off

subplot(3,3,5)
ix= find(gain(:,1)>.05 & gain(:,2)>.05);
plot(axx,nanmean(all_diag(ix,:),1),'b');
shading flat
hold on
plot(axx,nanmean(all_diag(ix,:),1)+nanstd(all_diag(ix,:),[],1)./sqrt(length(ix)),'b');
plot(axx,nanmean(all_diag(ix,:),1)-nanstd(all_diag(ix,:),[],1)./sqrt(length(ix)),'b');
plot([axx(1) axx(end)],[0 0],'--k','LineWidth',1)
a=axis;
plot([0 0],[a(3) a(4)],'--k','LineWidth',2)
hold off


subplot(3,3,6)
ix= find(gain(:,1)<-.05 & gain(:,2)>.05);
plot(axx,nanmean(all_diag(ix,:),1),'b');
shading flat
hold on
plot(axx,nanmean(all_diag(ix,:),1)+nanstd(all_diag(ix,:),[],1)./sqrt(length(ix)),'b');
plot(axx,nanmean(all_diag(ix,:),1)-nanstd(all_diag(ix,:),[],1)./sqrt(length(ix)),'b');
plot([axx(1) axx(end)],[0 0],'--k','LineWidth',1)
a=axis;
plot([0 0],[a(3) a(4)],'--k','LineWidth',2)
hold off

%%
% % % % % k=round(-length(axx)/2:length(axx)/2);
% % % % % for i=1:length(k)
% % % % %     ix= find(gain(:,1)<-.05 & gain(:,2)<-.05);
% % % % %     anti_diag(1,i)= nanmean(diag(imgaussfilt(squeeze(nanmean(all_jpsth_norm(ix,:,:),1)),1),k(i)));
% % % % %     ix= find(gain(:,1)>.05 & gain(:,2)>.05);
% % % % %     anti_diag(2,i)= nanmean(diag(imgaussfilt(squeeze(nanmean(all_jpsth_norm(ix,:,:),1)),1),k(i)));
% % % % %     ix= find(gain(:,1)<-.05 & gain(:,2)>.05);
% % % % %     anti_diag(4,i)= nanmean(diag(imgaussfilt(squeeze(nanmean(all_jpsth_norm(ix,:,:),1)),1),k(i)));
% % % % % end
% % % % % for i=1:size(anti_diag,1)
% % % % %     subplot(3,3,6+i)
% % % % %     plot(k*10,anti_diag(i,:),'k','LineWidth',2)
% % % % %     hold on
% % % % %     a=axis;
% % % % %     plot([0 0],a(3:4),'--k')
% % % % %     plot(a(1:2),[0 0],'--k')
% % % % %     hold off
% % % % % end
