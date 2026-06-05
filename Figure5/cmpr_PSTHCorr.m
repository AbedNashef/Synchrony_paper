%
clear
udir= 'E:/all_data/';
dates = dir([udir 'D*']);
bin2use=1;
addpath ../helper_functions/
tbf=1; taf=.5;
win2corr=.2; d2corr= .020; ax2corr=-tbf:d2corr:taf-win2corr;
c_all=0;
c_allv=0;
c_alla=0;
n=50;
nbin=1;
class2take={'PC'};
freq= 120;
L= ((taf+tbf)*freq);
move_thresh=0.2;
prob=4;
load ../LASSO/FullData_withAutoClusters.mat;
ax=-tbf*1000:bin2use:taf*1000; ax=ax(1:end-1); ax=ax./1000;
%%
nax= [-tbf:10/1000:taf];nax=nax(1:end-1);
t1= find(nax<-0.25); t2= find(nax>0 & nax<.25);
bhv_ax=[-tbf:1/freq:taf];
sm=3;
gwin=gausswin(sm);
gwin=gwin./sum(gwin);
tixx= find(nax>=-.25 & nax<=.1);
% [ res shift_predict psth_pred surprise_mat std_mat NUM_TRIAL pval_mat] = my_JPSTH(cut1, cut2, BIN_SIZE,gauss_filt_std)
%%
for j=1:length(dates)
    load([udir dates(j).name]);
    c_trials=0;
    index=strfind(dates(j).name,'_');
    mouse_id=dates(j).name(index(end-1)+1:index(end)+1);
    for i=1:length(ReachS)
        times = ReachS(i).filt_kin(:,1);
        ix= find(ReachS(i).out(:,2)>move_thresh,1,'first');
        [~,ix]= max(ReachS(i).out(:,6));

        endpoint_time=ReachS(i).out(ix,1);
        endpoint_time=ReachS(i).out(end,1);
        if isempty(endpoint_time)
            ix=find(ReachS(i).filt_kin(:,2)>move_thresh,1,'first');
            endpoint_time=ReachS(i).filt_kin(ix,1);
        end
        ix= find(ReachS(i).filt_kin(:,1)>=endpoint_time-tbf & ReachS(i).filt_kin(:,1)<=endpoint_time+taf);
        dt=median(diff(ReachS(i).filt_kin(ix,1)));
        tq=-tbf:1/freq:taf;
        t=ReachS(i).filt_kin(ix,1)-endpoint_time;
        xV=interp1(ReachS(i).filt_kin(:,1),ReachS(i).filt_kin(:,6),tq+endpoint_time);
        yV=interp1(ReachS(i).filt_kin(:,1),ReachS(i).filt_kin(:,7),tq+endpoint_time);
        V=interp1(ReachS(i).filt_kin(:,1),ReachS(i).filt_kin(:,5),tq+endpoint_time);
        if ~isempty(find(xV>50))
            ix=find(xV<50);
            V=V(ix);
            yV=yV(ix);
            xV=xV(ix);
            V=interp1(tq(ix)+endpoint_time,V,tq+endpoint_time);
            yV=interp1(tq(ix)+endpoint_time,yV,tq+endpoint_time);
            xV=interp1(tq(ix)+endpoint_time,xV,tq+endpoint_time);
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
                this_trc=cellData(cc).Bin10;
                index= find(this_trc(:,1)>=endpoint_time-tbf & this_trc(:,1)<=endpoint_time+taf);
                FR(c_this,c_trials,:)= this_trc(index,2);

                index= find(this_trc(:,1)>=endpoint_time-0.5 & this_trc(:,1)<=endpoint_time+0.5);
                FR2gain(c_this,c_trials,:)= this_trc(index,2);
                %                 [h,t]=findpeaks(xV);
                %                 maxV(c_trials,1)=h(find(bhv_ax(t)>0,1,'first'));
                %                 mix=t(find(bhv_ax(t)>0,1,'first'));
                %                 time_max(c_trials,1)=bhv_ax(mix);
                %                 h=find(islocalmin(xV));
                %                 minix=h(find(bhv_ax(h)>time_max(c_trials,1),1,'first'));
                %                 minV(c_trials,1)=xV(minix);
                %                 time_min(c_trials,1)=bhv_ax(minix);
                %                 xDec(c_trials)=(maxV(c_trials,1)-minV(c_trials,1))/(time_max(c_trials,1)-time_min(c_trials,1));

                endpoints(c_trials)=ReachS(i).out(end,2);
                [h,t]=findpeaks(xV);
                [~,mix]= max(xV(1:find(bhv_ax==0)));
                maxV(c_trials,1)=xV(mix);
                %                     mix=t(find(bhv_ax(t)>-.10,1,'first'));
                time_max(c_trials,1)=bhv_ax(mix);
                h=find(islocalmin(xV));
                minix=h(find(bhv_ax(h)>time_max(c_trials,1) & xV(h)<0,1,'first'));
                minV(c_trials,1)=xV(minix);
                time_min(c_trials,1)=bhv_ax(minix);
                xDec(c_trials)=(maxV(c_trials,1)-minV(c_trials,1))/(time_max(c_trials,1)-time_min(c_trials,1));

                [h,t]=findpeaks(yV);
                maxV(c_trials,2)=h(find(bhv_ax(t)>0,1,'first'));
                mix=t(find(bhv_ax(t)>0,1,'first'));
                time_max(c_trials,2)=bhv_ax(mix);
                h=find(islocalmin(yV));
                minix=h(find(bhv_ax(h)>time_max(c_trials,2),1,'first'));
                minV(c_trials,2)=yV(minix);
                time_min(c_trials,2)=bhv_ax(minix);
                yDec(c_trials)=(maxV(c_trials,2)-minV(c_trials,2))/(time_max(c_trials,2)-time_min(c_trials,2));
                Chs{c_this}= cellData(cc).Channels;
                all_gain(c_this)=cellData(cc).gain;
            end
        end
    end
    ax=ax(1:10:end);
    all_gain= AllData(j).Trial(1).clusters;
    overall_dec=sqrt(xDec.^2+yDec.^2);
    xDec=-1.*(xDec);
    %     xDec=abs(xDec);
    endpoints=endpoints-nanmean(endpoints);
    for i=1:size(FR,1)-1
        for ii=i+1:size(FR,1)
            ch1=Chs{i};
            ch2=Chs{ii};
            if ~isempty(ch1) & ~isempty(ch2)
                overlap= length(intersect(ch1,ch2))/length(unique([ch1,ch2]));
            else
                overlap=0;
            end
            if overlap==0 & ((all_gain(i)<=2 & all_gain(ii)<=2) | (all_gain(i)<=2 & all_gain(ii)<=2))
                % if overlap==0 & ((all_gain(i)==2 & all_gain(ii)==1) | (all_gain(i)==1 & all_gain(ii)==2))
                c_all=c_all+1;
                cell1= squeeze(FR(i,:,:));
                cell2= squeeze(FR(ii,:,:));
                tix= find(ax>-.15 & ax<.0);

                %                 sim_spikes=nansum((cell1(:,tix).*cell2(:,tix)),2);
                sim_spikes=nanmean((cell1(:,tix).*cell2(:,tix)),2)-nanmean(nanmean((cell1(:,tix).*cell2(:,tix)),2));

                [r_real(c_all),p_real(c_all)]= corr(sim_spikes,xDec(:));
                h=fitlm(maxV(:,1),xDec);
                se= h.Residuals.Raw;

                tix1=find(nax>-.15 & nax<0);
                for nn=1:prob
                    [~,endix]=sort(se);
                    l=floor(length(endix)/prob);
                    Pr_S1= cell1(endix((nn-1)*l+1:nn*l),tix1);
                    Pr_S2= cell2(endix((nn-1)*l+1:nn*l),tix1);
                    tmp1=nanmean(Pr_S1,1);
                    tmp2=nanmean(Pr_S2,1);
                    % tmp1= conv2(tmp1,gwin','same');
                    % tmp2= conv2(tmp2,gwin','same');
                    [all_R_gain(c_all,nn),~]= corr(tmp1',tmp2');
                    [all_xcorr_gain(c_all,nn,:),~]=xcorr(tmp1,tmp2,'normalized');
                    Pr_S1= cell1(endix((nn-1)*l+1:nn*l),tixx);
                    Pr_S2= cell2(endix((nn-1)*l+1:nn*l),tixx);
                    [tr_r_gain(c_all,nn,:,:),~]= corr(Pr_S1,Pr_S2);
                    FR4PCA_gain(c_all,nn,:)=nanmean([nanmean(cell1(endix((nn-1)*l+1:nn*l),:),1); nanmean(cell2(endix((nn-1)*l+1:nn*l),:),1)],1);

                    [~,endix]=sort(endpoints);
                    l=floor(length(endix)/prob);
                    Pr_S1= cell1(endix((nn-1)*l+1:nn*l),tix1);
                    Pr_S2= cell2(endix((nn-1)*l+1:nn*l),tix1);
                    tmp1=nanmean(Pr_S1,1);
                    tmp2=nanmean(Pr_S2,1);
                    % tmp1= conv2(tmp1,gwin','same');
                    % tmp2= conv2(tmp2,gwin','same');
                    [all_R_endpoint(c_all,nn),~]= corr(tmp1',tmp2');
                    [all_xcorr_end(c_all,nn,:),~]=xcorr(tmp1,tmp2,'normalized');
                    Pr_S1= cell1(endix((nn-1)*l+1:nn*l),tixx);
                    Pr_S2= cell2(endix((nn-1)*l+1:nn*l),tixx);
                    [tr_r_End(c_all,nn,:,:),~]= corr(Pr_S1,Pr_S2);
                    FR4PCA_end(c_all,nn,:)=nanmean([nanmean(cell1(endix((nn-1)*l+1:nn*l),:),1); nanmean(cell2(endix((nn-1)*l+1:nn*l),:),1)],1);

                    [~,endix]=sort(maxV(:,1));
                    l=floor(length(endix)/prob);
                    Pr_S1= cell1(endix((nn-1)*l+1:nn*l),tix1);
                    Pr_S2= cell2(endix((nn-1)*l+1:nn*l),tix1);
                    tmp1=nanmean(Pr_S1,1);
                    tmp2=nanmean(Pr_S2,1);
                    % tmp1= conv2(tmp1,gwin','same');
                    % tmp2= conv2(tmp2,gwin','same');
                    [all_R_peakV(c_all,nn),~]= corr(tmp1',tmp2');
                    [all_xcorr_peakv(c_all,nn,:),~]=xcorr(tmp1,tmp2,'normalized');
                    Pr_S1= cell1(endix((nn-1)*l+1:nn*l),tixx);
                    Pr_S2= cell2(endix((nn-1)*l+1:nn*l),tixx);
                    [tr_r_peakV(c_all,nn,:,:),~]= corr(Pr_S1,Pr_S2);
                    FR4PCA_peakV(c_all,nn,:)=nanmean([nanmean(cell1(endix((nn-1)*l+1:nn*l),:),1); nanmean(cell2(endix((nn-1)*l+1:nn*l),:),1)],1);

                    [~,endix]=sort(xDec);
                    l=floor(length(endix)/prob);
                    Pr_S1= cell1(endix((nn-1)*l+1:nn*l),tix1);
                    Pr_S2= cell2(endix((nn-1)*l+1:nn*l),tix1);
                    tmp1=nanmean(Pr_S1,1);
                    tmp2=nanmean(Pr_S2,1);
                    % tmp1= conv2(tmp1,gwin','same');
                    % tmp2= conv2(tmp2,gwin','same');
                    [all_R_xDec(c_all,nn),~]= corr(tmp1',tmp2');
                    [all_xcorr_xdec(c_all,nn,:),~]=xcorr(tmp1,tmp2,'normalized');
                    Pr_S1= cell1(endix((nn-1)*l+1:nn*l),tixx);
                    Pr_S2= cell2(endix((nn-1)*l+1:nn*l),tixx);
                    [tr_r_xdec(c_all,nn,:,:),~]= corr(Pr_S1,Pr_S2);
                    FR4PCA_dec(c_all,nn,:)=nanmean([nanmean(cell1(endix((nn-1)*l+1:nn*l),:),1); nanmean(cell2(endix((nn-1)*l+1:nn*l),:),1)],1);
                end

                N_trials(c_all)=length(endix);
                gain(c_all,:)=[all_gain(i) all_gain(ii)];
                date_id(c_all)= j;
                Dist4cells(c_all)= abs(cellData(i).depth-cellData(ii).depth);
                all_mice_id{c_all}=mouse_id;

            end
        end
    end
    clear FR FR2gain this_S nX nY nZ xVelocity yVelocity Speed xDec yDec endpoints
    clear time_min time_max minV maxV Chs
end
%%
ix1= find(gain(:,1)==1 & gain(:,2)==1);
ix2= find(gain(:,1)==2 & gain(:,2)==2);
ix3= find((gain(:,1)==2 & gain(:,2)==1) | (gain(:,1)==1 & gain(:,2)==2));
%%
for i=1:size(tr_r_xdec,1)
    for ii=1:size(tr_r_xdec,2)
        diag_end(i,ii,:)=diag(squeeze(tr_r_End(i,ii,:,:)));
        diag_xdec(i,ii,:)=diag(squeeze(tr_r_xdec(i,ii,:,:)));
        diag_peakV(i,ii,:)=diag(squeeze(tr_r_peakV(i,ii,:,:)));
        diag_gain(i,ii,:)=diag(squeeze(tr_r_gain(i,ii,:,:)));
    end
end
%%
ci=1.96;
for i=1:3
    this_ix= eval(['ix' num2str(i)]);
    subplot(3,4,(i-1)*4+1)
    plot(nax(tixx),diag(squeeze(nanmean(tr_r_End(this_ix,1,:,:),1))),'m')
    hold on
    plot(nax(tixx),diag(squeeze(nanmean(tr_r_End(this_ix,1,:,:),1)))-diag(squeeze(nanstd(tr_r_End(this_ix,1,:,:),[],1)))./sqrt(length(this_ix))*ci,'--m')
    plot(nax(tixx),diag(squeeze(nanmean(tr_r_End(this_ix,1,:,:),1)))+diag(squeeze(nanstd(tr_r_End(this_ix,1,:,:),[],1)))./sqrt(length(this_ix))*ci,'--m')
    plot(nax(tixx),diag(squeeze(nanmean(tr_r_End(this_ix,4,:,:),1))),'g')
    plot(nax(tixx),diag(squeeze(nanmean(tr_r_End(this_ix,4,:,:),1)))-diag(squeeze(nanstd(tr_r_End(this_ix,4,:,:),[],1)))./sqrt(length(this_ix))*ci,'--g')
    plot(nax(tixx),diag(squeeze(nanmean(tr_r_End(this_ix,4,:,:),1)))+diag(squeeze(nanstd(tr_r_End(this_ix,4,:,:),[],1)))./sqrt(length(this_ix))*ci,'--g')
    a=axis;
    plot(a(1:2),[0 0],'--k')
    plot([0 0],a(3:4),'--k')
    hold off
    set(gca,'TickDir','out'); box off
    xlabel('Time to endpoint (s)');
    ylabel('Corr')
    title(['endpoint (cluster' num2str(i) ')'])


    subplot(3,4,(i-1)*4+2)
    plot(nax(tixx),diag(squeeze(nanmean(tr_r_xdec(this_ix,1,:,:),1))),'g')
    hold on
    plot(nax(tixx),diag(squeeze(nanmean(tr_r_xdec(this_ix,1,:,:),1)))-diag(squeeze(nanstd(tr_r_xdec(this_ix,1,:,:),1)))./sqrt(length(this_ix))*ci,'--g')
    plot(nax(tixx),diag(squeeze(nanmean(tr_r_xdec(this_ix,1,:,:),1)))+diag(squeeze(nanstd(tr_r_xdec(this_ix,1,:,:),1)))./sqrt(length(this_ix))*ci,'--g')
    plot(nax(tixx),diag(squeeze(nanmean(tr_r_xdec(this_ix,4,:,:),1))),'m')
    plot(nax(tixx),diag(squeeze(nanmean(tr_r_xdec(this_ix,4,:,:),1)))-diag(squeeze(nanstd(tr_r_xdec(this_ix,4,:,:),1)))./sqrt(length(this_ix))*ci,'--m')
    plot(nax(tixx),diag(squeeze(nanmean(tr_r_xdec(this_ix,4,:,:),1)))+diag(squeeze(nanstd(tr_r_xdec(this_ix,4,:,:),1)))./sqrt(length(this_ix))*ci,'--m')
    a=axis;
    plot(a(1:2),[0 0],'--k')
    plot([0 0],a(3:4),'--k')
    hold off
    set(gca,'TickDir','out'); box off
    xlabel('Time to endpoint (s)');
    ylabel('Corr')
    title(['decel (cluster' num2str(i) ')'])


    subplot(3,4,(i-1)*4+3)
    plot(nax(tixx),diag(squeeze(nanmean(tr_r_peakV(this_ix,1,:,:),1))),'g')
    hold on
    plot(nax(tixx),diag(squeeze(nanmean(tr_r_peakV(this_ix,1,:,:),1)))-diag(squeeze(nanstd(tr_r_peakV(this_ix,1,:,:),1)))./sqrt(length(this_ix))*ci,'--g')
    plot(nax(tixx),diag(squeeze(nanmean(tr_r_peakV(this_ix,1,:,:),1)))+diag(squeeze(nanstd(tr_r_peakV(this_ix,1,:,:),1)))./sqrt(length(this_ix))*ci,'--g')
    plot(nax(tixx),diag(squeeze(nanmean(tr_r_peakV(this_ix,4,:,:),1))),'m')
    plot(nax(tixx),diag(squeeze(nanmean(tr_r_peakV(this_ix,4,:,:),1)))-diag(squeeze(nanstd(tr_r_peakV(this_ix,4,:,:),1)))./sqrt(length(this_ix))*ci,'--m')
    plot(nax(tixx),diag(squeeze(nanmean(tr_r_peakV(this_ix,4,:,:),1)))+diag(squeeze(nanstd(tr_r_peakV(this_ix,4,:,:),1)))./sqrt(length(this_ix))*ci,'--m')
    a=axis;
    plot(a(1:2),[0 0],'--k')
    plot([0 0],a(3:4),'--k')
    hold off
    set(gca,'TickDir','out'); box off
    xlabel('Time to endpoint (s)');
    ylabel('Corr')
    title(['peak v (cluster' num2str(i) ')'])

    subplot(3,4,(i-1)*4+4)
    plot(nax(tixx),diag(squeeze(nanmean(tr_r_gain(this_ix,1,:,:),1))),'g')
    hold on
    plot(nax(tixx),diag(squeeze(nanmean(tr_r_gain(this_ix,1,:,:),1)))-diag(squeeze(nanstd(tr_r_gain(this_ix,1,:,:),1)))./sqrt(length(this_ix))*ci,'--g')
    plot(nax(tixx),diag(squeeze(nanmean(tr_r_gain(this_ix,1,:,:),1)))+diag(squeeze(nanstd(tr_r_gain(this_ix,1,:,:),1)))./sqrt(length(this_ix))*ci,'--g')
    plot(nax(tixx),diag(squeeze(nanmean(tr_r_gain(this_ix,4,:,:),1))),'m')
    plot(nax(tixx),diag(squeeze(nanmean(tr_r_gain(this_ix,4,:,:),1)))-diag(squeeze(nanstd(tr_r_gain(this_ix,4,:,:),1)))./sqrt(length(this_ix))*ci,'--m')
    plot(nax(tixx),diag(squeeze(nanmean(tr_r_gain(this_ix,4,:,:),1)))+diag(squeeze(nanstd(tr_r_gain(this_ix,4,:,:),1)))./sqrt(length(this_ix))*ci,'--m')
    a=axis;
    plot(a(1:2),[0 0],'--k')
    plot([0 0],a(3:4),'--k')
    hold off
    set(gca,'TickDir','out'); box off
    xlabel('Time to endpoint (s)');
    ylabel('Corr')
    title(['decel error (cluster' num2str(i) ')'])

end
%%
ixx= find(nax(tixx)>=-0.15 & nax(tixx)<=0.);
for i=1:3
    this_ix= eval(['ix' num2str(i)]);
    subplot(3,4,(i-1)*4+1)

    histogram(squeeze(nanmean(diag_end(this_ix,1,ixx),3)),'BinWidth',.1,'FaceColor','m')
    hold on
    histogram(squeeze(nanmean(diag_end(this_ix,4,ixx),3)),'BinWidth',.1,'FaceColor','g')
    a=axis;
    n1=nanmean(squeeze(nanmean(diag_end(this_ix,1,ixx),3)));
    n2=nanmean(squeeze(nanmean(diag_end(this_ix,4,ixx),3)));
    plot([n1 n1],a(3:4),'-g')
    plot([n2 n2],a(3:4),'-m')
    hold off
    [~,p]=ttest(squeeze(nanmean(diag_end(this_ix,1,ixx),3)),squeeze(nanmean(diag_end(this_ix,4,ixx),3)));
    set(gca,'TickDir','out'); box off
    ylabel('#');
    xlabel('Corr')
    title(['endpoint (cluster' num2str(i) '); p=' num2str(p)])


    subplot(3,4,(i-1)*4+2)
    histogram(squeeze(nanmean(diag_xdec(this_ix,1,ixx),3)),'BinWidth',.1,'FaceColor','m')
    hold on
    histogram(squeeze(nanmean(diag_xdec(this_ix,4,ixx),3)),'BinWidth',.1,'FaceColor','g')
    a=axis;
    n1=nanmean(squeeze(nanmean(diag_xdec(this_ix,1,ixx),3)));
    n2=nanmean(squeeze(nanmean(diag_xdec(this_ix,4,ixx),3)));
    plot([n1 n1],a(3:4),'-g')
    plot([n2 n2],a(3:4),'-m')
    hold off
    [~,p]=ttest(squeeze(nanmean(diag_xdec(this_ix,1,ixx),3)),squeeze(nanmean(diag_xdec(this_ix,4,ixx),3)));
    set(gca,'TickDir','out'); box off
    ylabel('#');
    xlabel('Corr')
    title(['xdec (cluster' num2str(i) '); p=' num2str(p)])


    subplot(3,4,(i-1)*4+3)
    histogram(squeeze(nanmean(diag_peakV(this_ix,1,ixx),3)),'BinWidth',.1,'FaceColor','m')
    hold on
    histogram(squeeze(nanmean(diag_peakV(this_ix,4,ixx),3)),'BinWidth',.1,'FaceColor','g')
    a=axis;
    n1=nanmean(squeeze(nanmean(diag_peakV(this_ix,1,ixx),3)));
    n2=nanmean(squeeze(nanmean(diag_peakV(this_ix,4,ixx),3)));
    plot([n1 n1],a(3:4),'-g')
    plot([n2 n2],a(3:4),'-m')
    hold off
    [~,p]=ttest(squeeze(nanmean(diag_peakV(this_ix,1,ixx),3)),squeeze(nanmean(diag_peakV(this_ix,4,ixx),3)));
    set(gca,'TickDir','out'); box off
    ylabel('#');
    xlabel('Corr')
    title(['peak V (cluster' num2str(i) '); p=' num2str(p)])

    subplot(3,4,(i-1)*4+4)
    histogram(squeeze(nanmean(diag_gain(this_ix,1,ixx),3)),'BinWidth',.1,'FaceColor','g')
    hold on
    histogram(squeeze(nanmean(diag_gain(this_ix,4,ixx),3)),'BinWidth',.1,'FaceColor','m')
    a=axis;
    n1=nanmean(squeeze(nanmean(diag_gain(this_ix,1,ixx),3)));
    n2=nanmean(squeeze(nanmean(diag_gain(this_ix,4,ixx),3)));
    plot([n1 n1],a(3:4),'-g')
    plot([n2 n2],a(3:4),'-m')
    hold off
    [~,p]=ttest(squeeze(nanmean(diag_gain(this_ix,1,ixx),3)),squeeze(nanmean(diag_gain(this_ix,4,ixx),3)));
    set(gca,'TickDir','out'); box off
    ylabel('#');
    xlabel('Corr')
    title(['decel error (cluster' num2str(i) '); p=' num2str(p)])
end
%%
tix= find(nax>=-1 & nax<=1);
gwin=gausswin(21);
gwin=gwin./sum(gwin);
clear sFR4PCA_end sFR4PCA_dec sFR4PCA_peakV sFR4PCA_gain
for i=1:size(FR4PCA_dec,1)
    for ii=1:size(FR4PCA_peakV,2)
        tmp=conv2(squeeze(FR4PCA_end(i,ii,:)),gwin,'same');
        tmp=(tmp-mean(tmp(find(nax<-.5))))./std(tmp(find(nax<-.5)));
        sFR4PCA_end(i,ii,:)=tmp;
        tmp=conv2(squeeze(FR4PCA_dec(i,ii,:)),gwin,'same');
        tmp=(tmp-mean(tmp(find(nax<-.5))))./std(tmp(find(nax<-.5)));
        sFR4PCA_dec(i,ii,:)=tmp;%conv2(squeeze(FR4PCA_dec(i,ii,:)),gwin,'same');
        tmp=conv2(squeeze(FR4PCA_peakV(i,ii,:)),gwin,'same');
        tmp=(tmp-mean(tmp(find(nax<-.5))))./std(tmp(find(nax<-.5)));
        sFR4PCA_peakV(i,ii,:)=tmp;%conv2(squeeze(FR4PCA_peakV(i,ii,:)),gwin,'same');
        tmp=conv2(squeeze(FR4PCA_gain(i,ii,:)),gwin,'same');
        tmp=(tmp-mean(tmp(find(nax<-.5))))./std(tmp(find(nax<-.5)));
        sFR4PCA_gain(i,ii,:)=tmp;%conv2(squeeze(FR4PCA_gain(i,ii,:)),gwin,'same');
    end
end

sFR4PCA_end=sFR4PCA_end(:,:,tix);
sFR4PCA_dec=sFR4PCA_dec(:,:,tix);
sFR4PCA_peakV=sFR4PCA_peakV(:,:,tix);
sFR4PCA_gain=sFR4PCA_gain(:,:,tix);
%%
figure
t_index= find(nax(tix)>=-.15 & nax(tix)<=0);
for i=1:3
    this_ix= eval(['ix' num2str(i)]);
    fr4pc= [squeeze(sFR4PCA_end(this_ix,1,:)) squeeze(sFR4PCA_end(this_ix,4,:))];
    [coeff,score,latent] = pca(fr4pc');
    n=length(tix);
    subplot(3,4,(i-1)*4+1)
    plot3(score(1:n,1),score(1:n,2),score(1:n,3),'m')
    hold on
    plot3(score(1+t_index,1),score(1+t_index,2),score(1+t_index,3),'color',[0.5 0 0.5],'LineWidth',2)
    scatter3(score(n,1),score(n,2),score(n,3),40,'^m')
    plot3(score(1+n:end,1),score(1+n:end,2),score(1+n:end,3),'g')
    plot3(score(1+n+t_index,1),score(1+n+t_index,2),score(1+n+t_index,3),'color',[0 .5 0],'LineWidth',2)
    scatter3(score(end,1),score(end,2),score(end,3),40,'^g')
    hold off
    tmp1=[score(1:n,1),score(1:n,2),score(1:n,3)];
    tmp2=[score(1+n:end,1),score(1+n:end,2),score(1+n:end,3)];
    [r,p]=corr(tmp1(:),tmp2(:));
    xlabel('dim1'); ylabel('dim2'); zlabel('dim3');
    title(['end (cluster' num2str(i) ')'])

    fr4pc= [squeeze(sFR4PCA_dec(this_ix,4,:)) squeeze(sFR4PCA_dec(this_ix,1,:))];
    [coeff,score,latent] = pca(fr4pc');
    n=length(tix);
    subplot(3,4,(i-1)*4+2)
    plot3(score(1:n,1),score(1:n,2),score(1:n,3),'m')
    hold on
    plot3(score(1+t_index,1),score(1+t_index,2),score(1+t_index,3),'color',[0.5 0 0.5],'LineWidth',2)
    scatter3(score(n,1),score(n,2),score(n,3),40,'^m')
    plot3(score(1+n:end,1),score(1+n:end,2),score(1+n:end,3),'g')
    plot3(score(1+n+t_index,1),score(1+n+t_index,2),score(1+n+t_index,3),'color',[0 .5 0],'LineWidth',2)
    scatter3(score(end,1),score(end,2),score(end,3),40,'^g')
    hold off
    tmp1=[score(1:n,1),score(1:n,2),score(1:n,3)];
    tmp2=[score(1+n:end,1),score(1+n:end,2),score(1+n:end,3)];
    [r,p]=corr(tmp1(:),tmp2(:));
    xlabel('dim1'); ylabel('dim2'); zlabel('dim3');
    title(['decel  (cluster' num2str(i) ')'])

    fr4pc= [squeeze(sFR4PCA_peakV(this_ix,4,:)) squeeze(sFR4PCA_peakV(this_ix,1,:))];
    [coeff,score,latent] = pca(fr4pc');
    n=length(tix);
    subplot(3,4,(i-1)*4+3)
    plot3(score(1:n,1),score(1:n,2),score(1:n,3),'m')
    hold on
    plot3(score(1+t_index,1),score(1+t_index,2),score(1+t_index,3),'color',[0.5 0 0.5],'LineWidth',2)
    scatter3(score(n,1),score(n,2),score(n,3),40,'^m')
    plot3(score(1+n:end,1),score(1+n:end,2),score(1+n:end,3),'g')
    plot3(score(1+n+t_index,1),score(1+n+t_index,2),score(1+n+t_index,3),'color',[0 .5 0],'LineWidth',2)
    scatter3(score(end,1),score(end,2),score(end,3),40,'^g')
    hold off
    tmp1=[score(1:n,1),score(1:n,2),score(1:n,3)];
    tmp2=[score(1+n:end,1),score(1+n:end,2),score(1+n:end,3)];
    [r,p]=corr(tmp1(:),tmp2(:));
    xlabel('dim1'); ylabel('dim2'); zlabel('dim3');
    title(['peakv (cluster' num2str(i) ')'])

    fr4pc= [squeeze(sFR4PCA_gain(this_ix,4,:)) squeeze(sFR4PCA_gain(this_ix,1,:))];
    [coeff,score,latent] = pca(fr4pc');
    n=length(tix);
    subplot(3,4,(i-1)*4+4)
    plot3(score(1:n,1),score(1:n,2),score(1:n,3),'m')
    hold on
    plot3(score(1+t_index,1),score(1+t_index,2),score(1+t_index,3),'color',[0.5 0 0.5],'LineWidth',2)
    scatter3(score(n,1),score(n,2),score(n,3),40,'^m')
    plot3(score(1+n:end,1),score(1+n:end,2),score(1+n:end,3),'g')
    plot3(score(1+n+t_index,1),score(1+n+t_index,2),score(1+n+t_index,3),'color',[0 .5 0],'LineWidth',2)
    scatter3(score(end,1),score(end,2),score(end,3),40,'^g')
    hold off
    tmp1=[score(1:n,1),score(1:n,2),score(1:n,3)];
    tmp2=[score(1+n:end,1),score(1+n:end,2),score(1+n:end,3)];
    [r,p]=corr(tmp1(:),tmp2(:));
    xlabel('dim1'); ylabel('dim2'); zlabel('dim3');
    title(['decel error (cluster' num2str(i) ')'])
end
%%
dx= 50/1000;
dy=.5;
MAT1=squeeze(sFR4PCA_gain(ix3,1,:));
MAT2=squeeze(sFR4PCA_gain(ix3,4,:));
dax=nax;
day= floor(min([MAT1(:); MAT2(:)])):dy:ceil(max([MAT1(:); MAT2(:)]));
X1=nan(length(dax)-1,length(day));
X2=X1;
for i=1:length(dax)-1
    for ii=1:length(day)
        ix= find(nax>=dax(i) & nax<dax(i+1));
        tmp=MAT1(:,ix);
        X1(i,ii)= length(find(tmp>=day(ii) & tmp<day(ii)+dy));

        tmp=MAT2(:,ix);
        X2(i,ii)= length(find(tmp>=day(ii) & tmp<day(ii)+dy));
    end

end
subplot(2,3,1)
pcolor(dax(1:end-1),day,X1')
shading flat

subplot(2,3,2)
pcolor(dax(1:end-1),day,X2')
shading flat

subplot(2,3,3)
pcolor(dax(1:end-1),day,X1'-X2')
shading flat

subplot(2,2,3)
plot(dax(1:end-1),nanvar(X1,[],2))
hold on
plot(dax(1:end-1),nanvar(X2,[],2))
hold off
for  i=1:size(MAT1,1)
    r1(i)= corr(MAT1(i,find(nax>=-.5 & nax<0))',nanmean(MAT1(:,find(nax>=-.5 & nax<0)),1)');
    r2(i)= corr(MAT2(i,find(nax>=-.5 & nax<0))',nanmean(MAT2(:,find(nax>=-.5 & nax<0)),1)');
end