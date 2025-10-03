udir= 'E:/all_data/';
dates = dir([udir 'D*']);
bin2use=1;
addpath ../helper_functions/
addpath shuffle_control_progs\
tbf=.5; taf=.5;
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
load ../LASSO/FullData_withAutoClusters.mat;
%%
nax= [-0.5:1/1000:0.5];nax=nax(1:end-1);
t1= find(nax<-0.25); t2= find(nax>0 & nax<.25);
iterations=500;
gwin=gausswin(21); gwin=gwin./sum(gwin);
% [ res shift_predict psth_pred surprise_mat std_mat NUM_TRIAL pval_mat] = my_JPSTH(cut1, cut2, BIN_SIZE,gauss_filt_std)
%%
for c1=1:4
    c2=c1;
    disp(['Running ' num2str(c1) 'vs. ' num2str(c2)])
    c_all=0;
    for j=1:length(dates)
        load([udir dates(j).name]);
        c_trials=0;
        clusters= AllData(j).Trial(1).clusters;
        ntr=length(AllData(j).Trial);
        if length(clusters)==length(cellData)

            for i=1:length(ReachS)
                times = ReachS(i).filt_kin(:,1);
                ix= find(ReachS(i).out(:,2)>thresh,1,'first');
                [~,ix]=max(gradient(ReachS(i).out(:,6)));
                endpoint_time=ReachS(i).out(1,1);
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
                        if length(index)>length(nax)
                            index=index(1:end-1);
                        end
                        FR(c_this,c_trials,:)= this_trc(index,2);

                        all_gain(c_this)= clusters(cc);

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
                    if overlap==0 & ((all_gain(ii)==c1 & all_gain(i)==c2) | (all_gain(ii)==c2 & all_gain(i)==c1))
                        c_all=c_all+1;
                        cell1= squeeze(FR(i,:,:));
                        cell2= squeeze(FR(ii,:,:));
                        gain(c_all,:)=[all_gain(i) all_gain(ii)];
                        if gain(c_all,1)==c2 & gain(c_all,2)==c1
                            new_cell1=cell2;
                            new_cell2=cell1;
                            cell1=new_cell1; cell2=new_cell2;
                        end
                        [raw, shift_predict, pred, surprise, std, ~, ~] = my_JPSTH(cell1', cell2', nbin,1);
                        num_trials= size(cell1,1);
                        tmp=(raw-pred)./num_trials;%./std;
                        real_diag(c_all,:)= diag(imgaussfilt(tmp,2));
                        [~,tmp]=get_anti_diag(tmp,nax,-.1);
                        real_anti_diag(c_all,:)=tmp;
                        for iter=1:iterations
                            [~,rix]= sort(rand(1,size(cell1,1)));
                            % scell1=montecarlo_jitter_spikes(cell1,75);
                            scell1=cell1(rix,:);
                            [~,rix]= sort(rand(1,size(cell2,1)));
                            % scell2=montecarlo_jitter_spikes(cell2,75);
                            scell2=cell2(rix,:);
                            [raw, shift_predict, pred, surprise, std, ~, ~] = my_JPSTH(scell1', scell2', nbin,1);
                            % You can find my_JPSTH in
                            % %             all_jpsth(c_all,:,:)=raw-pred;
                            tmp=(raw-pred)./num_trials;%./std;
                            tmp_diag(iter,:)= diag(imgaussfilt(tmp,2));
                            [~,tmp_anti_diag(iter,:)]=get_anti_diag(tmp,nax,-.1);
                        end
                        anti_diag_ax=[-size(real_anti_diag,2)/2:size(real_anti_diag,2)/2];
                        MC_shuffle_diag(c_all,:)= nanmean(tmp_diag,1);
                        MC_shuffle_anti_diag(c_all,:)=nanmean(tmp_anti_diag,1);
                        date_id(c_all)= j;
                        p_value_anti_diag(c_all)=length(find(tmp_anti_diag(:,find(anti_diag_ax==0))>=real_anti_diag(c_all,find(anti_diag_ax==0))))/iterations;
                        % Dist4cells(c_all)= abs(cellData(i).depth-cellData(ii).depth);
                    end
                end
            end
        else
            disp([dates(j).name]);
        end

        clear FR FR2gain this_S nX nY nZ
    end
    subplot(4,4,(c1-1)*4+1)
    m=nanmean(real_diag,1);
    s=nanstd(real_diag,1)./sqrt(size(real_diag,1)); s=s.*1.96; % 5% confidence
    patch([nax fliplr(nax)],[m-s fliplr(m+s)],[.5 .5 1],'edgecolor','none')
    hold on
    plot(nax,m,'b','LineWidth',2)
    m=nanmean(MC_shuffle_diag,1);
    s=nanstd(MC_shuffle_diag,1)./sqrt(size(MC_shuffle_diag,1)); s=s.*1.96; % 5% confidence
    patch([nax fliplr(nax)],[m-s fliplr(m+s)],[.5 .5 0.5],'edgecolor','none')
    plot(nax,m,'k','LineWidth',2)
    hold off

    anti_diag_ax=[-size(real_anti_diag,2)/2:size(real_anti_diag,2)/2];
    anti_diag_ax=anti_diag_ax./1000; anti_diag_ax=anti_diag_ax(1:end-1);
    subplot(4,4,(c1-1)*4+2)
    m=nanmean(real_anti_diag,1);
    s=nanstd(real_anti_diag,1)./sqrt(size(real_anti_diag,1)); s=s.*1.96; % 5% confidence
    patch([anti_diag_ax fliplr(anti_diag_ax)],[m-s fliplr(m+s)],[.5 .5 1],'edgecolor','none')
    hold on
    plot(anti_diag_ax,m,'b','LineWidth',2)
    m=nanmean(MC_shuffle_anti_diag,1);
    s=nanstd(MC_shuffle_anti_diag,1)./sqrt(size(MC_shuffle_anti_diag,1)); s=s.*1.96; % 5% confidence
    patch([anti_diag_ax fliplr(anti_diag_ax)],[m-s fliplr(m+s)],[.5 .5 0.5],'edgecolor','none')
    plot(anti_diag_ax,m,'k','LineWidth',2)
    hold off

    subplot(4,4,(c1-1)*4+3)
    histogram(nanmean(real_diag(:,find(nax>-.15 & nax<0)),2),'FaceColor','b')
    hold on
    histogram(nanmean(MC_shuffle_diag(:,find(nax>-.15 & nax<0)),2),'FaceColor','k')
    hold off

    subplot(4,4,(c1-1)*4+4)
    histogram(real_anti_diag(:,find(anti_diag_ax==0)),'FaceColor','b')
    hold on
    histogram(MC_shuffle_diag(:,find(anti_diag_ax==0)),'FaceColor','k')
    hold off

    all_Sig(c1)=length(find(p_value_anti_diag<.05))/length(p_value_anti_diag);
    clear MC_shuffle_diag MC_shuffle_anti_diag real_diag real_anti_diag p_value_anti_diag

end