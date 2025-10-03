%
clear
udir= '/media/personlab/TARDIS/all_data/';
dates = dir([udir 'D*']);
bin2use=1;
addpath ../../helper_functions/
addpath ../../Rev1_analyses/
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
load ../../LASSO/FullData_withAutoClusters.mat;
%%
iterations=500; %per cluster
bin_size=20;
all_corr_coeff=[];
nax= [-0.5:1/1000:0.5];nax=nax(1:end-1);
t1= find(nax<-0.25); t2= find(nax>0 & nax<.25);
new_ax=ax(find(ax>=-tbf & ax<=taf));

% [ res shift_predict psth_pred surprise_mat std_mat NUM_TRIAL pval_mat] = my_JPSTH(cut1, cut2, BIN_SIZE,gauss_filt_std)
%%
tix= find(ax>=-tbf & ax<=taf);
iter=0;
for c=1:4
    iter=0;
    while iter<=iterations
        D= randi(length(AllData),1,2);
        nTr=min([length(AllData(D(1)).Trial) length(AllData(D(2)).Trial)]);
        C(1)= randi(size(AllData(D(1)).Trial(1).realNeu,1),1);
        C(2)= randi(size(AllData(D(2)).Trial(1).realNeu,1),1);
        if AllData(D(1)).Trial(1).clusters(C(1))==AllData(D(2)).Trial(1).clusters(C(2)) & AllData(D(1)).Trial(1).clusters(C(1))==c
            iter=iter+1;
            for i=1:nTr
                FR1(i,:)= AllData(D(1)).Trial(i).realNeu(C(1),tix);
                FR2(i,:)= AllData(D(2)).Trial(i).realNeu(C(2),tix);
            end
            
            gains(iter,1)=AllData(D(1)).Trial(1).clusters(C(1));
            gains(iter,2)=AllData(D(2)).Trial(1).clusters(C(2));
            if gains(iter,1)<gains(iter,2)
                [raw, shift_predict, pred, surprise, pst_std, ~, ~] = my_JPSTH(FR2', FR1', nbin,1);
            else
                [raw, shift_predict, pred, surprise, pst_std, ~, ~] = my_JPSTH(FR1', FR2', nbin,1);
            end
            rand_jpsth(iter,:,:)= (raw-pred)./nTr;
            all_fr(iter,:)= nanmean([FR1; FR2],1);
            clear  D C FR1 FR2;
        end
    end
    %%
    cix = find(gains(:,1)==c & gains(:,2)==c);
    X=smooth(nanmean(all_fr(cix,:),1),5).*1000; Y=smooth(diag(squeeze(nanmean(rand_jpsth(cix,:,:),1))),5);
    d= linspace(min(X),max(X),bin_size);
    for i=1:length(d)-1
        ix = find(X>=d(i) & X<d(i+1));
        x_lim(i)=mean(X(ix));
        y_lim(i)=mean(Y(ix));
        s_lim(i)=std(Y(ix))./sqrt(length(ix));
    end
    subplot(4,4,(c-1)*4+1)
    plot(ax(tix),squeeze(nanmean(all_fr(cix,:),1)));
    
    subplot(4,4,(c-1)*4+2)
    pcolor(ax(tix),ax(tix),squeeze(nanmean(rand_jpsth(cix,:,:),1))); shading flat
    
    subplot(4,4,(c-1)*4+3)
    scatter(X,Y,'.')
    xlabel('FR'); ylabel('cofiring');
    
    subplot(4,4,(c-1)*4+4)
    errorbar(x_lim,y_lim,s_lim,'LineWidth',2,'CapSize',0)
    xlabel('FR'); ylabel('cofiring');
    
    clear gains all_fr rand_jpsth
end