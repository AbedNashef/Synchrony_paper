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
prob=1/4;
load ../LASSO/FullData_withAutoClusters.mat;
ax=-tbf*1000:bin2use:taf*1000; ax=ax(1:end-1); ax=ax./1000;
%%
nax= [-tbf:1/1000:taf];nax=nax(1:end-1);
t1= find(nax<-0.25); t2= find(nax>0 & nax<.25);
bhv_ax=[-tbf:1/freq:taf];
iterations=100;

%%
all_fr=[];
all_sync=[];
unit_names=[];
all_isi=[];
all_real_isi=[];

all_cv=[];

all_cv2=[];


all_endpoints=[];
all_dec=[];
all_maxV=[];
all_SE_maxV2Dec=[];
all_outcome=[];
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

        %         endpoint_time=ReachS(i).out(ix,1);
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
        if ~vExclude & ~vStim% & ReachS(i).Outcome==0
            c_trials=c_trials+1;
            for cc=1:length(cellData)
                c_this=c_this+1;
                this_trc=cellData(cc).Bin1;
                index= find(this_trc(:,1)>=endpoint_time-tbf & this_trc(:,1)<=endpoint_time+taf);
                FR(c_this,c_trials,:)= this_trc(index,2);
                nTrc=get_null_fr(this_trc(:,2),1);
                nFR(c_this,c_trials,:)= nTrc(index);

                index= find(this_trc(:,1)>=endpoint_time-0.5 & this_trc(:,1)<=endpoint_time+0.5);
                FR2gain(c_this,c_trials,:)= this_trc(index,2);

                endpoints(c_trials)=ReachS(i).out(end,2);
                [h,t]=findpeaks(xV);
                [~,mix]= max(xV(1:find(bhv_ax==0)));
                maxV(c_trials,1)=xV(mix);
                %                     mix=t(find(bhv_ax(t)>-.10,1,'first'));
                time_max(c_trials,1)=bhv_ax(mix);
                h=find(islocalmin(xV));
                minix=h(find(bhv_ax(h)>time_max(c_trials,1),1,'first'));
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
                outcome(c_trials)=ReachS(i).Outcome;
            end
        end
    end
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
            if overlap==0 & all_gain(i)<=2 & all_gain(ii)<=2%((all_gain(i)==1 & all_gain(ii)==2) | (all_gain(i)==2 & all_gain(ii)==1)) % looking at 1 and 2 synchrony has the best results
                c_all=c_all+1;
                cell1= squeeze(FR(i,:,:));
                cell2= squeeze(FR(ii,:,:));
                ncell1= squeeze(nFR(i,:,:));
                ncell2= squeeze(nFR(ii,:,:));
                tix= find(ax>-.15 & ax<.0); %-100 to 300 can be salvaged.

                %                 sim_spikes=nansum((cell1(:,tix).*cell2(:,tix)),2);
                Pr_S1= nanmean(cell1(:,tix),1);
                Pr_S2= nanmean(cell2(:,tix),1);
                tmp=nansum((cell1(:,tix).*cell2(:,tix))-repmat(Pr_S1.*Pr_S2,size(cell1,1),1),2);
                %

                sim_spikes=nanmean((cell1(:,tix).*cell2(:,tix)),2)-nanmean((ncell1(:,tix).*ncell2(:,tix)),2);%nanmean(nanmean(cell1(:,tix),1).*nanmean(cell2(:,tix),1));
                % sim_spikes=tmp;
                fr=(nanmean(cell1(:,tix),2)+nanmean(cell2(:,tix),2))/2;
                mod1=nanmean(cell1(:,tix),2)-nanmean(nanmean(cell1(:,find(ax<-.5))));
                mod1=mod1;%./max(abs(mod1));
                mod2=nanmean(cell2(:,tix),2)-nanmean(nanmean(cell2(:,find(ax<-.5))));
                mod2=mod2;%./max(abs(mod2));
                fr=nanmean([mod1 mod2],2);
                fr=fr./max(abs(fr));
                cell12= cell1+cell2;
                cell12(find(cell12>1))=1;
                for cc=1:size(cell1,1)
                    ix= find(cell1(cc,:)==1);
                    dix1=diff(ix);
                    this_isi(cc,1)=nanmean(dix1);
                    

                    ix= find(cell2(cc,:)==1);
                    dix2=diff(ix);
                    this_isi(cc,2)=nanmean(dix2);


                    ix= find(cell12(cc,:)==1);
                    dix12=diff(ix);
                    this_isi(cc,3)=nanmean(dix12);

                    this_cv(cc,1)=nanstd(dix1)/nanmean(dix1);
                    this_cv(cc,2)=nanstd(dix2)/nanmean(dix2);
                    this_cv(cc,3)=nanstd(dix12)/nanmean(dix12);

                    isi=dix1(:); isi_tilde=circshift(dix1(:),1,1);
                    m1=2.*(abs(isi-isi_tilde))./(isi+isi_tilde);
                    this_cv2(cc,1)= sum(m1)/length(m1);
                    
                    isi=dix2(:); isi_tilde=circshift(dix2(:),1,1);
                    m2=2.*(abs(isi-isi_tilde))./(isi+isi_tilde);
                    this_cv2(cc,2)= sum(m2)/length(m2);

                    isi=dix12(:); isi_tilde=circshift(dix12(:),1,1);
                    m12=2.*(abs(isi-isi_tilde))./(isi+isi_tilde);
                    this_cv2(cc,3)= sum(m12)/length(m12);
                end
                all_real_isi=[all_real_isi; this_isi];
                this_isi=this_isi./repmat(max(this_isi,[],1),size(this_isi,1),1);
                this_cv=this_cv./repmat(max(this_cv,[],1),size(this_cv,1),1);
                this_cv2=this_cv2./repmat(max(this_cv2,[],1),size(this_cv2,1),1);

                all_isi=[all_isi; this_isi];
                all_cv=[all_cv; this_cv];
                all_cv2=[all_cv2; this_cv2];
                
                unit_names=[unit_names; ones(size(this_isi,1),1).*c_all];
                [r,p]=corr(fr(:),sim_spikes(:));
                all_r_fr2sync(c_all)=r;
                all_fr=[all_fr; fr(:)];
                all_sync=[all_sync; sim_spikes./max(abs(sim_spikes))];
               ixnan=find(~isnan(mean(this_isi,2)));
                for cc=1:size(all_isi,2)
                    [r1,p1]= corr(sim_spikes(ixnan),this_isi(ixnan,cc));
                    [r2,p2]= corr(fr(ixnan),this_isi(ixnan,cc));
                    all_r_isi(c_all,cc,:)=[r1 p1 r2 p2];

                    [r1,p1]= corr(sim_spikes(ixnan),this_cv(ixnan,cc));
                    [r2,p2]= corr(fr(ixnan),this_cv(ixnan,cc));
                    all_r_cv(c_all,cc,:)=[r1 p1 r2 p2];
                    
                    [r1,p1]= corr(sim_spikes(ixnan),this_cv2(ixnan,cc));
                    [r2,p2]= corr(fr(ixnan),this_cv2(ixnan,cc));
                    all_r_cv2(c_all,cc,:)=[r1 p1 r2 p2];
                end

                
            end
        end
    end
    clear FR FR2gain this_S nX nY nZ xVelocity yVelocity Speed xDec yDec endpoints
    clear time_min time_max minV maxV Chs nFR se maxV outcome
    clear this_isi this_cv this_cv2
end
%% binning parameters
bin_range=[.1 .9];
nbin= 6;
FRlim=[quantile(all_fr(:),bin_range(1)) quantile(all_fr(:),bin_range(2))];
dFR=linspace(FRlim(1),FRlim(2),nbin);
Synclim=[quantile(all_sync(:),bin_range(1)) quantile(all_sync(:),bin_range(2))];
dSync=linspace(Synclim(1),Synclim(2),nbin);
%%
isi1_matrix=nan(nbin-1,nbin-1);
cv_1_matrix=nan(nbin-1,nbin-1);
cv2_1_matrix=nan(nbin-1,nbin-1);

isi2_matrix=nan(nbin-1,nbin-1);
cv_2_matrix=nan(nbin-1,nbin-1);
cv2_2_matrix=nan(nbin-1,nbin-1);

real_isi_both_matrix=nan(nbin-1,nbin-1);
isi_both_matrix=nan(nbin-1,nbin-1);
cv_both_matrix=nan(nbin-1,nbin-1);
cv2_both_matrix=nan(nbin-1,nbin-1);

fr_values=[];
sync_values=[];
for i=1:nbin-1
    for ii=1:nbin-1
        ix= find(all_fr>=dFR(i) & all_fr<dFR(i+1) & all_sync>=dSync(ii) & all_sync<dSync(ii+1));% & all_outcome==0);
       
        isi1_matrix(i,ii)=nanmean(all_isi(ix,1));
        isi2_matrix(i,ii)=nanmean(all_isi(ix,2));
        isi_both_matrix(i,ii)=nanmean(all_isi(ix,3));

         real_isi_both_matrix(i,ii)=nanmean(all_real_isi(ix,3));
         real_isi_both_matrix_std(i,ii)=nanstd(all_real_isi(ix,3))./sqrt(length(ix));

        cv_1_matrix(i,ii)=nanmean(all_cv(ix,1));
        cv_2_matrix(i,ii)=nanmean(all_cv(ix,2));
        cv_both_matrix(i,ii)=nanmean(all_cv(ix,3));

        cv2_1_matrix(i,ii)=nanmean(all_cv2(ix,1));
        cv2_2_matrix(i,ii)=nanmean(all_cv2(ix,2));
        cv2_both_matrix(i,ii)=nanmean(all_cv2(ix,3));
        
        fr_values(i,ii)=dFR(i);
        sync_values(i,ii)=dSync(ii);
    end
end

%%
ix1= find(all_fr<=0.5 & all_sync>=.2);% & all_outcome==0);
ix2= find(all_fr<=0.5 & all_sync<=-.2);
%%
indx1=find(all_sync<mean(all_sync)-.5*std(all_sync));
indx2=find(all_sync>mean(all_sync)+.5*std(all_sync));
% endpoint
subplot(3,4,1)
imagesc(dFR(1:end-1),dSync(1:end-1),imgaussfilt(isi1_matrix,1)');
set(gca,'YDir','normal')
xlabel('normalized FR')
ylabel('normalized excess synchrony')
set(gca,'TickDir','out'); box off
title('ISI of 1st unit')

subplot(3,4,2)
imagesc(dFR(1:end-1),dSync(1:end-1),imgaussfilt(isi2_matrix,1)');
set(gca,'YDir','normal')
xlabel('normalized FR')
ylabel('normalized excess synchrony')
set(gca,'TickDir','out'); box off
title('ISI of 2nd unit')

subplot(3,4,3)
imagesc(dFR(1:end-1),dSync(1:end-1),imgaussfilt(isi_both_matrix,1)');
set(gca,'YDir','normal')
xlabel('normalized FR')
ylabel('normalized excess synchrony')
set(gca,'TickDir','out'); box off
title('ISI of both units')

subplot(3,4,4)
for i=1:size(isi_both_matrix,2)
   anti_diag(i,1)=isi1_matrix(i,end-i+1);
   anti_diag(i,2)=isi2_matrix(i,end-i+1);
   anti_diag(i,3)=isi_both_matrix(i,end-i+1);
end
plot(1:5,anti_diag)
ylabel('normalize ISI')
set(gca,'TickDir','out'); box off

subplot(3,4,5)
imagesc(dFR(1:end-1),dSync(1:end-1),imgaussfilt(cv_1_matrix,1)');
set(gca,'YDir','normal')
xlabel('normalized FR')
ylabel('normalized excess synchrony')
set(gca,'TickDir','out'); box off
title('CV of 1st unit')

subplot(3,4,6)
imagesc(dFR(1:end-1),dSync(1:end-1),imgaussfilt(cv_2_matrix,1)');
set(gca,'YDir','normal')
xlabel('normalized FR')
ylabel('normalized excess synchrony')
set(gca,'TickDir','out'); box off
title('CV of 2nd unit')

subplot(3,4,7)
imagesc(dFR(1:end-1),dSync(1:end-1),imgaussfilt(cv_both_matrix,1)');
set(gca,'YDir','normal')
xlabel('normalized FR')
ylabel('normalized excess synchrony')
set(gca,'TickDir','out'); box off
title('CV of both unit')

subplot(3,4,8)
for i=1:size(isi_both_matrix,2)
   anti_diag(i,1)=cv_1_matrix(i,end-i+1);
   anti_diag(i,2)=cv_2_matrix(i,end-i+1);
   anti_diag(i,3)=cv_both_matrix(i,end-i+1);
end
plot(1:5,anti_diag)
ylabel('normalize CV')
set(gca,'TickDir','out'); box off

subplot(3,4,9)
imagesc(dFR(1:end-1),dSync(1:end-1),imgaussfilt(cv2_1_matrix,1)');
set(gca,'YDir','normal')
xlabel('normalized FR')
ylabel('normalized excess synchrony')
set(gca,'TickDir','out'); box off
title('CV2 of 1st unit')

subplot(3,4,10)
imagesc(dFR(1:end-1),dSync(1:end-1),imgaussfilt(cv2_2_matrix,1)');
set(gca,'YDir','normal')
xlabel('normalized FR')
ylabel('normalized excess synchrony')
set(gca,'TickDir','out'); box off
title('CV2 of 2nd unit')

subplot(3,4,11)
imagesc(dFR(1:end-1),dSync(1:end-1),imgaussfilt(cv2_both_matrix,1)');
set(gca,'YDir','normal')
xlabel('normalized FR')
ylabel('normalized excess synchrony')
set(gca,'TickDir','out'); box off
title('CV2 of both unit')

subplot(3,4,12)
for i=1:size(isi_both_matrix,2)
   anti_diag(i,1)=cv2_1_matrix(i,end-i+1);
   anti_diag(i,2)=cv2_2_matrix(i,end-i+1);
   anti_diag(i,3)=cv2_both_matrix(i,end-i+1);
end
plot(1:5,anti_diag)
ylabel('normalize CV2')
set(gca,'TickDir','out'); box off

colormap(flipud(slanCM('RdYlBu',265)))
%%
clear ISI_high
uu=unique(unit_names);
for i=1:length(uu)
    ix2= find(all_fr>=dFR(1) & all_fr<dFR(3) & all_sync>=dSync(1) & all_sync<dSync(3) & unit_names==uu(i));
    ISI_high(i,1)=nanmean(all_real_isi(ix2,3));
        ix1= find(all_fr>=dFR(1) & all_fr<dFR(3) & all_sync>dSync(3) & all_sync<dSync(5) & unit_names==uu(i));
    ISI_high(i,2)=nanmean(all_real_isi(ix1,3));
end