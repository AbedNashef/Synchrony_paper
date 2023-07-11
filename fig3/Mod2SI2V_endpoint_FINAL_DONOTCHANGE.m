% Figure 3. For each trial calculate SI, modulation and behavioral kinematics
clear
udir= '../new_data/';
dates = dir([udir 'D*']);
bin2use=1;
addpath ../helper_functions/
tbf=1; taf=.5;
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
move_thresh=0.3;
prob=1/4;
%%
nax= [-0.5:1/1000:0.5];nax=nax(1:end-1);
t1= find(nax<-0.25); t2= find(nax>0 & nax<.25);
bhv_ax=[-tbf:1/freq:taf];

all_dec=[];
all_SI=[];
all_mxV=[];
all_endpoint=[];
all_cell1=[];
all_cell2=[];
date_id=[];
gain1=[]; gain2=[];
cell1_id=[];
cell2_id=[];
all_mice_id=[];
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
        endpoint_time=ReachS(i).out(ix,1);
        endpoint_time=ReachS(i).out(end,1);
        endpoint_v=[ReachS(i).out(end,6) ReachS(i).out(end,7)];
        if isempty(endpoint_time)
            ix=find(ReachS(i).filt_kin(:,2)>move_thresh,1,'first');
            endpoint_time=ReachS(i).filt_kin(ix,1);
        end
        times= ReachS(i).filt_kin(:,1);
        P1= repmat(ReachS(i).filt_kin(1,2:4),size(ReachS(i).filt_kin,1),1);
        tt=times(find(times>=endpoint_time-tbf & times<=endpoint_time+taf));
        
        ix= find(ReachS(i).filt_kin(:,1)>=endpoint_time-tbf & ReachS(i).filt_kin(:,1)<=endpoint_time+taf);
        dt=median(diff(ReachS(i).filt_kin(ix,1)));
        tq=-tbf:1/freq:taf;
        t=ReachS(i).filt_kin(ix,1)-endpoint_time;
        xV=interp1(t,ReachS(i).filt_kin(ix,6),tq);
        yV=interp1(t,ReachS(i).filt_kin(ix,7),tq);
        V=interp1(t,ReachS(i).filt_kin(ix,5),tq);
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
                
                endpoint_XY(c_trials,:)=[ReachS(i).out(end,2) ReachS(i).out(end,3)];
                endpoint_D(c_trials)= sqrt(ReachS(i).out(end,2)^2+ReachS(i).out(end,3)^2);
                
                T= ReachS(i).out(end,1);
                ix= find(ReachS(i).filt_kin(:,1)==T);
                index= find(ReachS(i).filt_kin(:,1)>=T-.3 & ReachS(i).filt_kin(:,1)<=T+.3);
                maxV(c_trials,:)= [max(ReachS(i).filt_kin(index,5)) max(ReachS(i).filt_kin(index,6)) max(ReachS(i).filt_kin(index,7)) max(ReachS(i).filt_kin(index,8))];
                all_gain(c_this)=cellData(cc).gain;
                
                [h,t]=findpeaks(xV);
                maxV(c_trials,1)=h(find(bhv_ax(t)>0,1,'first'));
                mix=t(find(bhv_ax(t)>0,1,'first'));
                time_max(c_trials,1)=bhv_ax(mix);
                h=find(islocalmin(xV));
                minix=h(find(bhv_ax(h)>time_max(c_trials,1),1,'first'));
                minV(c_trials,1)=xV(minix);%endpoint_v(1);%
                time_min(c_trials,1)=bhv_ax(minix);
                xDec(c_trials)=(maxV(c_trials,1)-minV(c_trials,1))/(time_max(c_trials,1)-time_min(c_trials,1));%-endpoint_time);%
                
                [h,t]=findpeaks(yV);
                maxV(c_trials,2)=h(find(bhv_ax(t)>0,1,'first'));
                mix=t(find(bhv_ax(t)>0,1,'first'));
                time_max(c_trials,2)=bhv_ax(mix);
                h=find(islocalmin(yV));
                minix=h(find(bhv_ax(h)>time_max(c_trials,2),1,'first'));
                minV(c_trials,2)=yV(minix);%endpoint_v(1);%
                time_min(c_trials,2)=bhv_ax(minix);
                yDec(c_trials)=(maxV(c_trials,2)-minV(c_trials,2))/(time_max(c_trials,2)-time_min(c_trials,2));%-endpoint_time);%
                
                Chs{c_this}= cellData(cc).Channels;
            end
        end
    end
    endpoint_XY(:,1)=endpoint_XY(:,1);%./max(endpoint_XY(:,1));
    endpoint_XY(:,2)=endpoint_XY(:,2);%./max(endpoint_XY(:,2));
    endpoint_XY=endpoint_XY-repmat(nanmean(endpoint_XY,1),size(endpoint_XY,1),1); 
    maxV=maxV./repmat(max(maxV,[],1),size(maxV,1),1);
    endpoint_D=endpoint_D-nanmean(endpoint_D);
    xDec=abs(xDec); xDec=xDec./max(xDec);
    yDec=abs(yDec); yDec=yDec./max(yDec);
    for i=1:size(FR,1)-1
        for ii=i+1:size(FR,1)
            ch1=Chs{i};
            ch2=Chs{ii};
            if ~isempty(ch1) & ~isempty(ch2)
                overlap= length(intersect(ch1,ch2))/length(unique([ch1,ch2]));
            else
                overlap=0;
            end
            if overlap==0
                c_all=c_all+1;
                cell1= squeeze(FR(i,:,:));
                cell2= squeeze(FR(ii,:,:));
                
                tix= find(ax>-.3 & ax<.3);
                Pr_S1= nanmean(cell1(:,tix),1);
                Pr_S2= nanmean(cell2(:,tix),1);
                Pr_S1S2= nanmean(cell1(:,tix).*cell2(:,tix),1);
                sumSI= nansum((cell1(:,tix).*cell2(:,tix))./repmat(Pr_S1.*Pr_S2,size(all_gain,1),1),2);
                %             sumSI= max((cell1(:,tix).*cell2(:,tix))./repmat(Pr_S1.*Pr_S2,size(all_gain,1),1),[],2);
                %             tmp=(cell1(:,tix).*cell2(:,tix))./repmat(Pr_S1.*Pr_S2,size(all_gain,1),1);
                %             sumSI= tmp(:,find(ax==-.05));
                all_dec=[all_dec; xDec(:) yDec(:)];
                all_SI=[all_SI; sumSI./max(sumSI)];
                all_mxV=[all_mxV; maxV];
                all_endpoint=[all_endpoint; endpoint_XY endpoint_D(:)];
                tmp=(sum(cell1(:,tix),2)-sum(cell1(:,find(ax<-.5)),2));
%                 tmp=(tmp)./max(tmp);
                all_cell1=[all_cell1; tmp./max(abs(tmp))];
                tmp=(sum(cell2(:,tix),2)-sum(cell2(:,find(ax<-.5)),2));
%                 tmp=(tmp)./max(tmp);
                all_cell2=[all_cell2; tmp./max(abs(tmp))];
                date_id= [date_id; ones(size(sumSI)).*j];
                cell1_id=[cell1_id; ones(size(sumSI)).*i];
                cell2_id=[cell2_id; ones(size(sumSI)).*ii];
                
                gain1=[gain1; ones(length(sumSI),1).*all_gain(i)];
                gain2=[gain2; ones(length(sumSI),1).*all_gain(ii)];
                Dist4cells(c_all)= abs(cellData(i).depth-cellData(ii).depth);
                all_mice_id=[all_mice_id; repmat(mat2cell(mouse_id,1),length(sumSI),1)];
            end
        end
    end
    clear FR FR2gain this_S nX nY nZ endpoint_XY endpoint_D maxV xDec yDec Chs
end
%%
orig_SI= all_SI;
all_mod= nanmean([all_cell1 all_cell2],2);
orig_mod= all_mod;
%%
norm_SI= all_SI;%./max(all_SI);%range(all_SI);
norm_mod = all_mod;%./max(all_mod);%range(all_mod);
all_SI=norm_SI;
all_mod=norm_mod;
ix= find(gain1<-.05 & gain2<-.05);
pausers_SI= all_SI(ix);
pausers_mod= all_mod(ix);
pausers_endpoint= all_endpoint(ix);
pausers_maxV= all_mxV(ix,:);
pausers_dec= all_dec(ix,:);
pausers_days= date_id(ix);
pausers_cell1= cell1_id(ix);
pausers_cell2= cell2_id(ix);
pausers_mouse= all_mice_id(ix);
pausers_SI=pausers_SI;%./max(pausers_SI);
pausers_mod=pausers_mod;%./max(pausers_mod);
pausers_maxV=pausers_maxV;%./max(pausers_maxV);
pausers_dec=pausers_dec;%./max(pausers_dec);

% [b,a]= regress(all_endpoint(ix,1),[nanmean([all_cell1(ix) all_cell2(ix)],2),all_SI(ix),ones(length(ix),1)]);
% % occupational matrix
% nbin=50; lims=[-1.2 1.2]; win= range(lims)/nbin;
% dend= [lims(1):range(lims)/nbin:lims(2)];
% MAT1=[]; MAT2=[];
% for i=1:length(dend)
%     for j=1:length(dend)
%         endix = find(all_endpoint(ix,1)>=dend(i) & all_endpoint(ix,1)<dend(i)+win & all_endpoint(ix,2)>=dend(j) & all_endpoint(ix,2)<dend(j)+win);
%         MAT1(i,j)=nanmedian(norm_SI(endix));%./length(endix);
%         MAT2(i,j)=nanmedian(norm_mod(endix));%./length(endix);
%     end
% end

%%
sm=1;
top_qut=.99; bot_qut=.05;
nbin=10;
SIlim=[quantile(pausers_SI,bot_qut) quantile(pausers_SI,top_qut)];%00];
dSI=[SIlim(1):range(SIlim)/nbin:SIlim(2)];
win_SI=range(SIlim)/nbin;
modlim=[quantile(pausers_mod,bot_qut) quantile(pausers_mod,top_qut)];
dMod=[modlim(1):range(modlim)/nbin:modlim(2)];
win_mod=range(modlim)/nbin;
endMAT=[]; Pmat=[]; stdMAT=[];
for i=1:length(dSI)
    for ii=1:length(dMod)
        ix = find(pausers_SI>= dSI(i) & pausers_SI<dSI(i)+win_SI & pausers_mod>= dMod(ii) & pausers_mod<dMod(ii)+win_mod);
        endMAT(i,ii)= nanmedian(pausers_endpoint(ix,1));%./log(length(ix));
        stdMAT(i,ii)= nanstd(pausers_endpoint(ix,1))./sqrt(length(ix));%./log(length(ix));
        [~,Pmat(i,ii)]= ttest(pausers_endpoint(ix,1),0);
%         subplot(length(dSI),length(dMod),(length(dSI)*(i-1))+ii)
%         histogram(pausers_endpoint(ix,1),'Normalization','probability')
%         hold on
%         a=axis;
%         plot([nanmedian(pausers_endpoint(ix,1)) nanmedian(pausers_endpoint(ix,1))],a(3:4))
%         hold off
    end
end
Pausers_mat=imgaussfilt(endMAT,sm);
Pausers_std_mat=imgaussfilt(stdMAT,sm);
pausers_binP_end= Pmat;

subplot 231
scatter(pausers_mod,pausers_SI,20,pausers_endpoint(:,1),'s','fill')
caxis([-0.2 .2])
xlabel('Modulation'); ylabel('\SigmaSI');
set(gca,'TickDir','out'); box off

subplot 234
pcolor(dMod,dSI,imgaussfilt(endMAT,.5))
shading flat
% caxis([-0.3 .3])
xlabel('Modulation'); ylabel('\SigmaSI');
set(gca,'TickDir','out'); box off

%%
% SIlim=[0 .4]%00];
% dSI=[0:range(SIlim)/nbin:SIlim(2)];
% win_SI=range(SIlim)/nbin;
% modlim=[-.5 .5];
% dMod=[modlim(1):range(modlim)/nbin:modlim(2)];
win_mod=range(modlim)/nbin;
endMAT=[]; Pmat=[]; stdMAT=[];
for i=1:length(dSI)
    for ii=1:length(dMod)
        ix = find(pausers_SI>= dSI(i) & pausers_SI<dSI(i)+win_SI & pausers_mod>= dMod(ii) & pausers_mod<dMod(ii)+win_mod);
        endMAT(i,ii)= nanmedian(pausers_maxV(ix,1));%./log(length(ix));
        [~,Pmat(i,ii)]= ttest(pausers_maxV(ix,1),nanmean(pausers_maxV(:,1)));
    end
end
Pausers_mat_V=imgaussfilt(endMAT,sm);
pausers_binP_maxv= Pmat;

subplot 232
scatter(pausers_mod,pausers_SI,10,pausers_maxV(:,1),'o','fill')
caxis([-0.2 .2])
xlabel('Modulation'); ylabel('\SigmaSI');
set(gca,'TickDir','out'); box off

subplot 235
pcolor(dMod,dSI,imgaussfilt(endMAT,.5))
shading flat
% caxis([-0.3 .3])
xlabel('Modulation'); ylabel('\SigmaSI');
set(gca,'TickDir','out'); box off

%%

% SIlim=[0 .4]%00];
% dSI=[0:range(SIlim)/nbin:SIlim(2)];
% win_SI=range(SIlim)/nbin;
% modlim=[-.5 .5];
% dMod=[modlim(1):range(modlim)/nbin:modlim(2)];
win_mod=range(modlim)/nbin;
endMAT=[]; Pmat=[]; stdMAT=[];
for i=1:length(dSI)
    for ii=1:length(dMod)
        ix = find(pausers_SI>= dSI(i) & pausers_SI<dSI(i)+win_SI & pausers_mod>= dMod(ii) & pausers_mod<dMod(ii)+win_mod);
        endMAT(i,ii)= nanmedian(pausers_dec(ix,1));%./log(length(ix));
        stdMAT(i,ii)= nanstd(pausers_dec(ix,1))./sqrt(length(ix));%./log(length(ix));
        [~,Pmat(i,ii)]= ttest(pausers_dec(ix,1),nanmean(pausers_dec(:,1)));
    end
end
Pausers_mat_dec=imgaussfilt(endMAT,sm);
Pausers_stdmat_dec=imgaussfilt(stdMAT,sm);
pausers_binP_dec= Pmat;

subplot 233
scatter(pausers_mod,pausers_SI,10,pausers_dec(:,1),'o','fill')
caxis([-0.2 .2])
xlabel('Modulation'); ylabel('\SigmaSI');
set(gca,'TickDir','out'); box off

subplot 236
pcolor(dMod,dSI,imgaussfilt(endMAT,.5))
shading flat
% caxis([-0.3 .3])
xlabel('Modulation'); ylabel('\SigmaSI');
set(gca,'TickDir','out'); box off

%%
sm=1;
figure
subplot(2,3,1)
pcolor(dMod,dSI,imgaussfilt(Pausers_mat,sm))
shading flat
% caxis([-0.3 .3])
xlabel('Modulation'); ylabel('\SigmaSI');
set(gca,'TickDir','out'); box off

subplot(2,3,2)
pcolor(dMod,dSI,imgaussfilt(Pausers_mat_V,sm))
shading flat
% caxis([-0.3 .3])
xlabel('Modulation'); ylabel('\SigmaSI');
set(gca,'TickDir','out'); box off

subplot(2,3,3)
main_diag=[diag(Pausers_mat_V) diag(Pausers_mat)];
scatter(main_diag(:,1),main_diag(:,2),20,'ob','fill')
ix= find(~isnan(mean(main_diag,2)));
[r,p]= corr(main_diag(ix,1),main_diag(ix,2));
h=fitlm(main_diag(ix,1),main_diag(ix,2),'poly1');
slope= h.Coefficients.Estimate(2);
Lg{1}= ['Diag; ' num2str(slope)];

for i=1:size(Pausers_mat,1)
    anti_diag(i,:)=[Pausers_mat_V(i,end-i+1) Pausers_mat(i,end-i+1)];
end
hold on
scatter(anti_diag(:,1),anti_diag(:,2),20,'ok','fill')
ix= find(~isnan(mean(anti_diag,2)));
[r,p]= corr(anti_diag(ix,1),anti_diag(ix,2));
h=fitlm(anti_diag(ix,1),anti_diag(ix,2),'poly1');
slope= h.Coefficients.Estimate(2);
Lg{2}= ['Anti Diag; ' num2str(slope)];

ix= round(size(Pausers_mat,1)/2);
clampSI= [Pausers_mat_V(ix,:)' Pausers_mat(ix,:)'];
scatter(clampSI(:,1),clampSI(:,2),20,'om','fill')
ix= find(~isnan(mean(clampSI,2)));
[r,p]= corr(clampSI(ix,1),clampSI(ix,2));
h=fitlm(clampSI(ix,1),clampSI(ix,2),'poly1');
slope= h.Coefficients.Estimate(2);
Lg{3}= ['SI clamp; ' num2str(slope)];

ix= round(size(Pausers_mat,1)/2);
clampMod= [Pausers_mat_V(:,ix) Pausers_mat(:,ix)];
scatter(clampMod(:,1),clampMod(:,2),20,'or','fill')
ix= find(~isnan(mean(clampMod,2)));
[r,p]= corr(clampMod(ix,1),clampMod(ix,2));
h=fitlm(clampMod(ix,1),clampMod(ix,2),'poly1');
slope= h.Coefficients.Estimate(2);
Lg{4}= ['mod clamp; ' num2str(slope)];
l=lsline;
h=get(gca,'Children')
h(1).Color='b'; h(2).Color='k';
h(3).Color='m'; h(4).Color='r';
legend(Lg)
h=get(gca,'Legend');
h.String=h.String(1:4)

set(gca,'TickDir','out'); box off
xlabel('max V');
ylabel('\DeltaEndpoint')
%%

subplot(2,3,4)
pcolor(dMod,dSI,imgaussfilt(Pausers_mat,sm))
shading flat
% caxis([-0.3 .3])
xlabel('Modulation'); ylabel('\SigmaSI');
set(gca,'TickDir','out'); box off

subplot(2,3,5)
pcolor(dMod,dSI,imgaussfilt(Pausers_mat_dec,sm))
shading flat
% caxis([-0.3 .3])
xlabel('Modulation'); ylabel('\SigmaSI');
set(gca,'TickDir','out'); box off

subplot(2,3,6)
main_diag=[diag(Pausers_mat_dec) diag(Pausers_mat)];
scatter(main_diag(:,1),main_diag(:,2),20,'ob','fill')
ix= find(~isnan(mean(main_diag,2)));
[r,p]= corr(main_diag(ix,1),main_diag(ix,2));
h=fitlm(main_diag(ix,1),main_diag(ix,2),'poly1');
slope= h.Coefficients.Estimate(2);
Lg{1}= ['Diag; ' num2str(slope)];

for i=1:size(Pausers_mat,1)
    anti_diag(i,:)=[Pausers_mat_dec(i,end-i+1) Pausers_mat(i,end-i+1)];
end
hold on
scatter(anti_diag(:,1),anti_diag(:,2),20,'ok','fill')
ix= find(~isnan(mean(anti_diag,2)));
[r,p]= corr(anti_diag(ix,1),anti_diag(ix,2));
h=fitlm(anti_diag(ix,1),anti_diag(ix,2),'poly1');
slope= h.Coefficients.Estimate(2);
Lg{2}= ['Anti-Diag; ' num2str(slope)];

ix= round(size(Pausers_mat,1)/2);
clampSI= [Pausers_mat_dec(ix,:)' Pausers_mat(ix,:)'];
scatter(clampSI(:,1),clampSI(:,2),20,'om','fill')
ix= find(~isnan(mean(clampSI,2)));
[r,p]= corr(clampSI(ix,1),clampSI(ix,2));
h=fitlm(clampSI(ix,1),clampSI(ix,2),'poly1');
slope= h.Coefficients.Estimate(2);
Lg{3}= ['SI clamp; ' num2str(slope)];

ix= round(size(Pausers_mat,1)/2);
clampMod= [Pausers_mat_dec(:,ix) Pausers_mat(:,ix)];
scatter(clampMod(:,1),clampMod(:,2),20,'or','fill')
ix= find(~isnan(mean(clampMod,2)));
[r,p]= corr(clampMod(ix,1),clampMod(ix,2));
h=fitlm(clampMod(ix,1),clampMod(ix,2),'poly1');
slope= h.Coefficients.Estimate(2);
Lg{4}= ['Mod clamp; ' num2str(slope)];
l=lsline;
h=get(gca,'Children');
h(1).Color='b'; h(2).Color='k';
h(3).Color='m'; h(4).Color='r';
legend(Lg)
h=get(gca,'Legend');
h.String=h.String(1:4)

set(gca,'TickDir','out'); box off
xlabel('dec');
ylabel('\DeltaEndpoint')

colormap((brewermap([],'PRGn')))
%%
%% test significance of matrices
% 1- randomize the matrices
figure
clear Ppausers
Mat2test=[];
iterations=100;
X= Pausers_mat(:);
for i=1:iterations
    [~,rix] = sort(rand(1,length(X)));
    tmp= X(rix);
    Mat2test(i,:,:)= reshape(tmp,size(Pausers_mat));
end
Ppausers=[];
for i=1:size(Pausers_mat,1)
    for ii=1:size(Pausers_mat,2)
        [~,Ppausers(i,ii)]= ttest(squeeze(Mat2test(:,i,ii)),Pausers_mat(i,ii));
    end
end
subplot(3,4,1)
pcolor(dMod,dSI,imgaussfilt(-log(Ppausers),.5))
shading flat
set(gca,'TickDir','out'); box off
xlabel('Mod'); ylabel('SI');

Mat2test=[];
X= Pausers_mat_V(:);
for i=1:iterations
    [~,rix] = sort(rand(1,length(X)));
    tmp= X(rix);
    Mat2test(i,:,:)= reshape(tmp,size(Pausers_mat_V));
end
PV=[];
for i=1:size(Pausers_mat_V,1)
    for ii=1:size(Pausers_mat_V,2)
        [~,PV(i,ii)]= ttest(squeeze(Mat2test(:,i,ii)),Pausers_mat_V(i,ii));
    end
end
subplot(3,4,5)
pcolor(dMod,dSI,imgaussfilt(-log(PV),.5));
shading flat
set(gca,'TickDir','out'); box off
xlabel('Mod'); ylabel('SI');

X= Pausers_mat_dec(:);
for i=1:iterations
    [~,rix] = sort(rand(1,length(X)));
    tmp= X(rix);
    Mat2test(i,:,:)= reshape(tmp,size(Pausers_mat_dec));
end
Pdec=[];
for i=1:size(Pausers_mat_dec,1)
    for ii=1:size(Pausers_mat_dec,2)
        [~,Pdec(i,ii)]= ttest(squeeze(Mat2test(:,i,ii)),Pausers_mat_dec(i,ii));
    end
end
subplot(3,4,9)
pcolor(dMod,dSI,imgaussfilt(-log(Pdec),.5));
shading flat
set(gca,'TickDir','out'); box off
xlabel('Mod'); ylabel('SI');


C=colormap(hot);
colormap(flipud(C));
% above vs. under diag
above=[]; under=[];
for i=1:size(Pausers_mat,1)
    under= [under Pausers_mat(i,i+1:end)];
end
for i=1:size(Pausers_mat,2)
    above= [above; Pausers_mat(i+1:end,i)];
end
subplot(3,4,2)
histogram(above,20,'FaceColor','r','Normalization','probability')
hold on
histogram(under,20,'FaceColor','b','Normalization','probability')
a=axis;
plot([nanmean(above) nanmean(above)],a(3:4),'r','LineWidth',2)
plot([nanmean(under) nanmean(under)],a(3:4),'b','LineWidth',2)
plot([0 0],a(3:4),'--k','LineWidth',1)
hold off
[~,p]= ttest2(above,under);
title(p);
set(gca,'TickDir','out'); box off
xlabel('\DeltaEndpoint'); ylabel('Fraction');

above=[]; under=[];
for i=1:size(Pausers_mat_V,1)
    under= [under Pausers_mat_V(i,i+1:end)];
end
for i=1:size(Pausers_mat_V,2)
    above= [above; Pausers_mat_V(i+1:end,i)];
end
subplot(3,4,6)
histogram(above,20,'FaceColor','r','Normalization','probability')
hold on
histogram(under,20,'FaceColor','b','Normalization','probability')
a=axis;
plot([nanmean(above) nanmean(above)],a(3:4),'r','LineWidth',2)
plot([nanmean(under) nanmean(under)],a(3:4),'b','LineWidth',2)
plot([0 0],a(3:4),'--k','LineWidth',1)
hold off
[~,p]= ttest2(above,under);
title(p);
set(gca,'TickDir','out'); box off
xlabel('\DeltaEndpoint'); ylabel('Fraction');

above=[]; under=[];
for i=1:size(Pausers_mat_dec,1)
    under= [under Pausers_mat_dec(i,i+1:end)];
end
for i=1:size(Pausers_mat_dec,2)
    above= [above; Pausers_mat_dec(i+1:end,i)];
end
subplot(3,4,10)
histogram(above,20,'FaceColor','r','Normalization','probability')
hold on
histogram(under,20,'FaceColor','b','Normalization','probability')
a=axis;
plot([nanmean(above) nanmean(above)],a(3:4),'r','LineWidth',2)
plot([nanmean(under) nanmean(under)],a(3:4),'b','LineWidth',2)
plot([0 0],a(3:4),'--k','LineWidth',1)
hold off
[~,p]= ttest2(above,under);
title(p);
set(gca,'TickDir','out'); box off
xlabel('\DeltaEndpoint'); ylabel('Fraction');

clear Lg;
subplot(3,4,3)
X=[]; Y=[];
for i=1:size(Pausers_mat,1)
    Y=[Y; Pausers_mat(i,:)];
end
scatter(dMod,nanmedian(Y,1),30,'ok','fill');
[r,p]= corr(dMod(:),nanmean(Y,1)');
Lg{1}=['Mod; r=' num2str(r) '; p=' num2str(p)];

X=[]; Y=[];
for i=1:size(Pausers_mat,2)
    Y=[Y; Pausers_mat(:,i)'];
end
hold on
scatter(dSI,nanmedian(Y,1),30,'or','fill');
[r,p]= corr(dSI(:),nanmean(Y,1)');
hold off
Lg{2}=['SI; r=' num2str(r) '; p=' num2str(p)];
legend(Lg)

subplot(3,4,7)
X=[]; Y=[];
for i=1:size(Pausers_mat_V,1)
    Y=[Y; Pausers_mat_V(i,:)];
end
scatter(dMod,nanmedian(Y,1),30,'ok','fill');
[r,p]= corr(dMod(:),nanmean(Y,1)');
Lg{1}=['Mod; r=' num2str(r) '; p=' num2str(p)];

X=[]; Y=[];
for i=1:size(Pausers_mat_V,2)
    Y=[Y; Pausers_mat_V(:,i)'];
end
hold on
scatter(dSI,nanmedian(Y,1),30,'or','fill');
[r,p]= corr(dSI(:),nanmean(Y,1)');
hold off
Lg{2}=['SI; r=' num2str(r) '; p=' num2str(p)];
legend(Lg)

subplot(3,4,11)
X=[]; Y=[];
for i=1:size(Pausers_mat_dec,1)
    Y=[Y; Pausers_mat_dec(i,:)];
end
scatter(dMod,nanmedian(Y,1),30,'ok','fill');
[r,p]= corr(dMod(:),nanmean(Y,1)');
Lg{1}=['Mod; r=' num2str(r) '; p=' num2str(p)];

X=[]; Y=[];
for i=1:size(Pausers_mat_dec,2)
    Y=[Y; Pausers_mat_dec(:,i)'];
end
hold on
scatter(dSI,nanmedian(Y,1),30,'or','fill');
[r,p]= corr(dSI(:),nanmean(Y,1)');
hold off
Lg{2}=['SI; r=' num2str(r) '; p=' num2str(p)];
legend(Lg)

C=brewermap((nbin-1)*2+1,'RdYlBu');
k=[-(nbin-1):(nbin-1)];
for i=1:length(k)
    subplot(3,12,[10 11])
    tmp= diag(Pausers_mat,k(i));
    plot(1:length(tmp),tmp,'Color',C(i,:),'LineWidth',1);
    hold on
    
    subplot(3,12,12)
    scatter(randn(1)/20,nanmean(tmp),30,'Marker','o','MarkerFaceColor',C(i,:),'MarkerEdgeColor','none')
    hold on
    
    subplot(3,12,[22 23])
    tmp= diag(Pausers_mat_V,k(i));
    plot(1:length(tmp),tmp,'Color',C(i,:),'LineWidth',1);
    hold on
    
    subplot(3,12,24)
    scatter(randn(1)/20,nanmean(tmp),30,'Marker','o','MarkerFaceColor',C(i,:),'MarkerEdgeColor','none')
    hold on
    
    subplot(3,12,[34 35])
    tmp= diag(Pausers_mat_dec,k(i));
    plot(1:length(tmp),tmp,'Color',C(i,:),'LineWidth',1);
    hold on
    
    subplot(3,12,36)
    scatter(randn(1)/20,nanmean(tmp),30,'Marker','o','MarkerFaceColor',C(i,:),'MarkerEdgeColor','none')
    hold on
end

    subplot(3,12,[10 11])
set(gca,'TickDir','out'); box off
xlabel('Bin#'); ylabel('\Delta Endpoint');
a=axis;
plot([a(1) a(2)],[0 0],'--k');
hold off

    subplot(3,12,[22 23])
set(gca,'TickDir','out'); box off
xlabel('Bin#'); ylabel('\Delta Endpoint');
a=axis;
plot([a(1) a(2)],[0 0],'--k');
hold off

    subplot(3,12,[34 35])
set(gca,'TickDir','out'); box off
xlabel('Bin#'); ylabel('\Delta Endpoint');
a=axis;
plot([a(1) a(2)],[0 0],'--k');
hold off

    subplot(3,12,12)
set(gca,'TickDir','out','XTick',[]); box off
xlabel(''); ylabel('\Delta Endpoint');
a=axis;
plot([a(1) a(2)],[0 0],'--k');
hold off
h=get(gca,'Children');
[r,p]= corr([1:length(h)-1]',[h(2:end).YData]');
title(['r=' num2str(r) '; p=' num2str(p)]);

    subplot(3,12,24)
set(gca,'TickDir','out','XTick',[]); box off
xlabel(''); ylabel('\Delta Endpoint');
a=axis;
plot([a(1) a(2)],[0 0],'--k');
hold off
h=get(gca,'Children');
[r,p]= corr([1:length(h)-1]',[h(2:end).YData]');
title(['r=' num2str(r) '; p=' num2str(p)]);

subplot(3,12,36)
set(gca,'TickDir','out','XTick',[]); box off
xlabel(''); ylabel('\Delta Endpoint');
a=axis;
plot([a(1) a(2)],[0 0],'--k');
hold off
h=get(gca,'Children');
[r,p]= corr([1:length(h)-1]',[h(2:end).YData]');
title(['r=' num2str(r) '; p=' num2str(p)]);
%% for final figure

subplot(2,3,1)
pcolor(dMod,dSI,imgaussfilt(Pausers_mat_dec,sm))
shading flat
% caxis([-0.3 .3])
xlabel('Modulation'); ylabel('\SigmaSI');
set(gca,'TickDir','out'); box off
title('Deceleration')

subplot(2,3,2)
pcolor(dMod,dSI,imgaussfilt(Pausers_mat,sm))
shading flat
% caxis([-0.3 .3])
xlabel('Modulation'); ylabel('\SigmaSI');
set(gca,'TickDir','out'); box off
title('\DeltaEndpoint')

subplot(4,3,3)
for i=1:size(Pausers_mat,1)
    anti_diag(i,:)=[Pausers_mat_dec(i,end-i+1) Pausers_mat(i,end-i+1)];
    std_anti_diag(i,:)=[Pausers_stdmat_dec(i,end-i+1) Pausers_std_mat(i,end-i+1)];
end
errorbar(dMod,anti_diag(:,1),std_anti_diag(:,1),'Color','r')
hold on
errorbar(dMod,nanmean(Pausers_mat_dec,1),nanstd(Pausers_mat_dec,[],2),'Color',[1 0.7 0.3])
errorbar(dSI,nanmean(Pausers_mat_dec,2),nanstd(Pausers_mat_dec,[],1),'Color',[1 0.3 0.7])
h=fitlm(dSI,anti_diag(:,1),'poly1');
slope= h.Coefficients.Estimate(2);
r2=h.Rsquared.Ordinary;
Lg{1}= ['Anti-Diag; ' num2str(slope) '(' num2str(r2) ')'];

h=fitlm(dMod,nanmean(Pausers_mat_dec,1),'poly1');
slope= h.Coefficients.Estimate(2);
r2=h.Rsquared.Ordinary;
Lg{2}= ['Mod; ' num2str(slope) '(' num2str(r2) ')'];

h=fitlm(dSI,nanmean(Pausers_mat_dec,2),'poly1');
slope= h.Coefficients.Estimate(2);
r2=h.Rsquared.Ordinary;
Lg{3}= ['SI; ' num2str(slope) '(' num2str(r2) ')'];
legend(Lg);
hold off

subplot(4,3,6)
errorbar(dMod,anti_diag(:,2),std_anti_diag(:,2),'Color','k')
hold on
errorbar(dMod,nanmean(Pausers_mat,1),nanstd(Pausers_mat,[],2),'Color',[.7 0.7 0.7])
errorbar(dSI,nanmean(Pausers_mat,2),nanstd(Pausers_mat,[],1),'Color',[.3 0.3 0.3])
h=fitlm(dSI,anti_diag(:,2),'poly1');
slope= h.Coefficients.Estimate(2);
r2=h.Rsquared.Ordinary;
Lg{1}= ['Anti-Diag; ' num2str(slope) '(' num2str(r2) ')'];

h=fitlm(dMod,nanmean(Pausers_mat,1),'poly1');
slope= h.Coefficients.Estimate(2);
r2=h.Rsquared.Ordinary;
Lg{2}= ['Mod; ' num2str(slope) '(' num2str(r2) ')'];

h=fitlm(dSI,nanmean(Pausers_mat,2),'poly1');
slope= h.Coefficients.Estimate(2);
r2=h.Rsquared.Ordinary;
Lg{3}= ['SI; ' num2str(slope) '(' num2str(r2) ')'];
legend(Lg);
hold off
