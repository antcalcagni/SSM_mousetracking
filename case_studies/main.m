%% Set environment
%Real data: one categorical variable with four levels [Barca & Pezzullo, 2012]
clear all
addpath(genpath('SSM_MouseTracking/'))
load('barca2012_rawdata.mat'); load('stimuli.mat'); load('covariate.mat')

%% Define full structure of data (with missing data w.r.t. stimuli)

% original stimuli used during the experiment
stm_HF = unique(cellstr(stimuli(DATA(:,3)==500))); 
stm_LF = unique(cellstr(stimuli(DATA(:,3)==501))); 
stm_NW = unique(cellstr(stimuli(DATA(:,3)==502))); 
stm_PW = unique(cellstr(stimuli(DATA(:,3)==503))); 

sbjs = unique(DATA(:,1));
X_coord = zeros(length(sbjs),24*4,101); Y_coord = zeros(length(sbjs),24*4,101); Conds = zeros(length(sbjs),24*4); Stms = cell(length(sbjs),24*4);

for i=1:length(sbjs)    
    XX = DATA(DATA(:,1)==sbjs(i),:);
    xstm = stimuli(DATA(:,1)==sbjs(i));

    k=0;
    x = unique(cellstr(xstm(XX(:,3)==500,:))); %HF condition
    for j=1:length(stm_HF)
        k=k+1;
        [~,idx] = ismember(stm_HF{j},x);
        if idx~=0
            X_coord(i,k,:) = XX(idx,7:107);
            Y_coord(i,k,:) = XX(idx,108:end);
            Conds(i,k) = XX(idx,3);
            Stms(i,k) = cellstr(x(idx));
        else
            X_coord(i,k,:) = 999;
            Y_coord(i,k,:) = 999;
            Conds(i,k) = 999;
            Stms(i,k) = cellstr('NA');
        end
    endfind(D1(:,2)==1)

    k=length(stm_HF);
    x = unique(cellstr(xstm(XX(:,3)==501,:))); %LF condition
    for j=1:length(stm_LF)
        k=k+1;
        [~,idx] = ismember(stm_LF{j},x);
        if idx~=0
            X_coord(i,k,:) = XX(idx,7:107);
            Y_coord(i,k,:) = XX(idx,108:end);
            Conds(i,k) = XX(idx,3);
            Stms(i,k) = cellstr(x(idx));
        else
            X_coord(i,k,:) = 999;
            Y_coord(i,k,:) = 999;
            Conds(i,k) = 999;
            Stms(i,k) = cellstr('NA');
        end
    end

    k=length(stm_HF)+length(stm_LF);
    x = unique(cellstr(xstm(XX(:,3)==502,:))); %NW condition
    for j=1:length(stm_NW)
        k=k+1;
        [~,idx] = ismember(stm_NW{j},x);
        if idx~=0
            X_coord(i,k,:) = XX(idx,7:107);
            Y_coord(i,k,:) = XX(idx,108:end);
            Conds(i,k) = XX(idx,3);
            Stms(i,k) = cellstr(x(idx));
        else
            X_coord(i,k,:) = 999;
            Y_coord(i,k,:) = 999;
            Conds(i,k) = 999;
            Stms(i,k) = cellstr('NA');
        end
    end

    k=length(stm_HF)+length(stm_LF)+length(stm_NW);
    x = unique(cellstr(xstm(XX(:,3)==503,:))); %PW condition
    for j=1:length(stm_PW)
        k=k+1;
        [~,idx] = ismember(stm_PW{j},x);
        if idx~=0
            X_coord(i,k,:) = XX(idx,7:107);
            Y_coord(i,k,:) = XX(idx,108:end);
            Conds(i,k) = XX(idx,3);
            Stms(i,k) = cellstr(x(idx));
        else
            X_coord(i,k,:) = 999;
            Y_coord(i,k,:) = 999;
            Conds(i,k) = 999;
            Stms(i,k) = cellstr('NA');
        end
    end
    end
end


%% Substitute missing data w.r.t. stimuli with subject-based average trajectory within the same class of stimuli
X_coord_red = zeros(length(sbjs),24*4,101); Y_coord_red = zeros(length(sbjs),24*4,101); Conds_red = zeros(length(sbjs),24*4); Stms_red = cell(length(sbjs),24*4);

xcond = [repmat(500,length(stm_HF),1);repmat(501,length(stm_LF),1);repmat(502,length(stm_NW),1);repmat(503,length(stm_PW),1)];

X_coord_red=X_coord; Y_coord_red=Y_coord;

for i=1:size(X_coord,1)
    iid = find(Conds(i,:)==999); iid_conds = xcond(iid);
    for j=1:length(iid)
        iidx=find(xcond==iid_conds(j));
        XX = squeeze(X_coord(i,iidx,:)); XX(find(XX(:,1)==999),:)=[];
        YY = squeeze(Y_coord(i,iidx,:)); YY(find(YY(:,1)==999),:)=[];    
        X_coord_red(i,iid(j),:) = mean(XX,1);
        Y_coord_red(i,iid(j),:) = mean(YY,1);
    end
end


%% Recover information on covariates
stm_ALL = [stm_HF' stm_LF' stm_NW' stm_PW'];
iid=zeros(length(stm_ALL),1);
for j=1:length(stm_ALL)
    iid(j) = find(strcmp(covariate_data.stimulus,stm_ALL{j}));
end
bigF_cov = covariate_data.bigFreq(iid);
AoA_cov = [covariate_data.EA_varlex(iid(xcond==500)); covariate_data.EA_varlex(iid(xcond==501))]; AoA_cov(isnan(AoA_cov)) = mean(AoA_cov,'omitnan');
IMM_cov = [covariate_data.IMM_varlex(iid(xcond==500)); covariate_data.IMM_varlex(iid(xcond==501))]; IMM_cov(isnan(IMM_cov)) = mean(IMM_cov,'omitnan');


%% Prepare data for the analyses and setting up variables for running the model

% Rotate trajectories into first quadrant
for i=1:size(X_coord_red,1)
    for j=1:size(X_coord_red,2)
        if X_coord_red(i,j,end) < 0 
            X_coord_red(i,j,:) = X_coord_red(i,j,:)*-1;
        end
        %hold on;scatter(X_coord_red(i,j,:),Y_coord_red(i,j,:),'k');scatter(X_coord_red(i,j,end),Y_coord_red(i,j,end),'r.')
    end
end

% Replace near zero coordinates with first-movement coordinates
for i=1:size(X_coord_red,1)
    for j=1:size(X_coord_red,2)
        [x,y] = refineTrajs(X_coord_red(i,j,:),Y_coord_red(i,j,:));
        X_coord_red(i,j,:)=x;
        Y_coord_red(i,j,:)=y;        
    end
end

% Removing first flat coordinates
alpha=31;
X_coord_red = X_coord_red(:,:,alpha:end); Y_coord_red = Y_coord_red(:,:,alpha:end); 

% Compute empirical angles 
Y = atan2(Y_coord_red,X_coord_red); 
%plot_theta(vec(Y),[])

% Determine sample inormation
alpha=71; %determine the area of the last click [T-alpha, T]
%plot_theta(vec(Y(:,:,alpha:end)),[])
I=size(Y,1); J=size(Y,2); T=size(Y,3);
mu1=mean(vec(Y(:,:,alpha:end))); mu2=(pi/2)+(pi/2)-mu1; %angle points for the two labels --our range is [pi,0]
thr = pi/2;
plot_theta(vec(Y(:,:,alpha:end)),mu1)

% Compute Z-values for vonMises categorization
z=zeros(I*J*T,1); %z2=zeros(I*J*size(Y2,3),1); 
z(Y(:)>thr)=1; %z2(Y2(:)>thr)=1;
Z=reshape(z,I,J,T);

% ML estimates of vonMises concentrations
yy=Y(:);fi
kappa1 = est_kappa(yy(z==1),0);
kappa2 = est_kappa(yy(z==0),0);

% Fix remaining parameters
sigmax = ones(I,1)*5.5; 
a=ones(J,1); %slopes do not play a role in this version of the model
bnd=5; %this bound is enough for expressing the model

% Define classification variable
D1 = double([xcond==500 xcond==501 xcond==502 xcond==503]);

%% Save variables
clearvars -except Y Z D1 X_coord Y_coord X_coord_red Y_coord_red Conds Stms stm_HF stm_LF stm_NW stm_PW stm_ALL AoA_cov IMM_cov bigF_cov kappa1 kappa2 mu1 mu2 sigmax a bnd I J T thr
save('barca2012_nomissin_full.mat')

%% Analyses
%%% see case_studies/run_analyses.m

%% Results: MCMC diagnostics (1)
load('mh_4_1.mat')

data.Z=Z; data.Y=Y; data.D1=D1; pars.sigmax=sigmax; pars.kappa1=kappa1; pars.kappa2=kappa2; pars.a=a; 
pars.mu1=mu1; pars.mu2=mu2; pars.bnd=bnd; data.x=scaledata(bigF_cov,-3,3);

[S_coll, post_means] = MCMC_diag(mh_4_1,500,1,true,false,'categorical',4);

[b,gamma1,eta,delta] = switchModel('categorical',data.D1,data.x,post_means');

tic;[XF,PF] = GaussApprox_filter(data.Z,pars.sigmax,pars.a,b,0);toc

figure()
for i=1:I, subplot(4,6,i); hold on; 
    plot(XF(i,:),'b'); xlim([-0.5 T+3]); ylim([-bnd-2 bnd+2]); hline(0,'--k'); title(['Sbj ' num2str(i)]); end

iid_w = [min(find(D1(:,1)==1)) min(find(D1(:,2)==1)) min(find(D1(:,3)==1)) min(find(D1(:,4)==1))];
P = twoPL(XF,pars.a,b);
figure()
cls = {'r','b','k','g'};
for i=1:I
    k=0;
    for j=1:length(iid_w)
        k=k+1;
        subplot(4,6,i);hold on; plot(squeeze(P(i,iid_w(j),:)),cls{k}); 
        xlim([-0.5 T+3]); ylim([-0.1 1.1]); hline(0.5,'--k'); title(['Sbj ' num2str(i)])
    end
end
legend('HF','LF','NW','PW')


figure();
for k=1:4
    subplot(2,2,k);
    iidx = find(data.D1(:,k)==1);
    for i=1:I    
        hold on
        for j=1:length(iidx),plot(squeeze(X_coord_red(i,iidx(j),:)),squeeze(Y_coord_red(i,iidx(j),:)));end
    end
end

figure();
for k=1:4
    subplot(2,2,k);
    iidx = find(D1(:,k)==1);
    for i=1:I    
        hold on
        for j=1:length(iidx),plot(squeeze(Y(i,iidx(j),:)));end
    end
end

%%
u=0;
figure()
II=9;
for i=1:II
    for k=1:4
        u=u+1;        
        iidx = find(D1(:,k)==1);                   
        subplot(II,4,u); hold on;for j=1:length(iidx),plot(squeeze(X_coord_red(i,iidx(j),:)),squeeze(Y_coord_red(i,iidx(j),:))); xlim([-1 1]); vline(0,'--k');  end
    end    
end

%%
u=0;
figure()
II=9;
for i=1:II
    for k=1:4
        u=u+1;        
        iidx = find(D1(:,k)==1);                        
        subplot(II,4,u); hold on;for j=1:length(iidx),plot(squeeze(Y(i,iidx(j),:))); hline(thr,'--k'); ylim([0 2.5]);  end
    end    
end

%% Results: model 5_0 (regression)
clear all
information_function = @(xdom,gamma_est,a) cell2mat(arrayfun(@(j)exp(bsxfun(@times,bsxfun(@minus,xdom,gamma_est(j)),-a(j)))./...
    (1+exp(bsxfun(@times,bsxfun(@minus,xdom,gamma_est(j)),-a(j)))).^2,1:length(gamma_est),'UniformOutput',false)');

load('barca2012_nomissin_full.mat');load('mh_5_0.mat')

q0=1;q=72; %select just {HF,LF,NW}
data.Z=Z(:,q0:q,:); data.Y=Y(:,q0:q,:); data.D1=D1(q0:q,1:3); pars.sigmax=sigmax; pars.kappa1=kappa1; pars.kappa2=kappa2; pars.a=a(q0:q); 
pars.mu1=mu1; pars.mu2=mu2; pars.bnd=bnd; data.x=scaledata(bigF_cov(q0:q),-3,3);


[S_coll, post_means,Rgelm,~,S] = MCMC_diag(mh_5_0,1000,1,false,true,'regression',1);
b = data.x*post_means;
tic;[XF,PF] = GaussApprox_filter(data.Z,pars.sigmax,pars.a,b,0);toc 
P = twoPL(XF,pars.a,b);


%%%%%%%%%%%%% Figure for paper %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mh=mh_5_0; Rgelm=round(Rgelm,3); mms = round(mean(S_coll),2); fnts=20; fntst=22; fntl=14; 

figure(); nfig=1;
set(gcf,'units','points','position',[100,100,1500,400],'color','w','defaultTextInterpreter','latex','defaultLegendInterpreter','latex');

subplot1 = subplot(1,3,1);hold on;for i=1:size(S,3),plot(S(:,1,i),'Parent',subplot1);end; xlim([-250 size(S,1)+250]); 
set(subplot1,'FontSize',fntl,'FontName','Times');title('$\eta$','FontSize',fntst); xlabel(['Iterations ($\hat R = ' num2str(Rgelm(1)) '$)'],'FontSize',fnts); ylabel('Traceplot','FontSize',fnts)

[f,xi] = ksdensity(S_coll); xhdi = mbe_hdi(S_coll,.95); iid=min(find(round(xi,3)==round(xhdi(1),3))):min(find(round(xi,3)==round(xhdi(2),3)));
subplot2 = subplot(1,3,2);plot(xi,f,'LineWidth',1.5);box('off'); hold on; area(xi,f,'FaceColor',[0.5843 0.8157 0.9882]); plot(xi(iid),f(iid));area(xi(iid),f(iid),'FaceColor',[0 0.4470 0.7410]);
text(median(xi),quantile(f,0.5),'95% HDPI','FontSize',fntl,'HorizontalAlignment','center')
set(subplot2,'FontSize',fntl,'FontName','Times');title('$\eta$','FontSize',fntst);xlabel(['$\mu = ' num2str(mms(1)) '$'],'FontSize',fnts); ylabel('Marginal posterior density','FontSize',fnts);
box('off')

subplot3 = subplot(1,3,3);hold on; [acfs,lags] = autocorr(S_coll); area([-0.2 max(lags)+1],[0.25 0.25],'EdgeColor','w','FaceColor',rgb('lightgray')); area([-0.2 max(lags)+1],[-0.25 -0.25],'EdgeColor','w','FaceColor',rgb('lightgray'));
bar(lags,acfs,0.1,'FaceColor',[0 0.4470 0.7410]); plot(lags,acfs,'o','Color',[0 0.4470 0.7410],'MarkerFaceColor',[0 0.4470 0.7410]); ylim([-1 1]); xlim([-0.5 max(lags)+1]); box('off'); 
set(subplot3,'FontSize',fntl,'FontName','Times');title('$\eta$','FontSize',fntst); xlabel('lags','FontSize',fnts); ylabel('ACF','FontSize',fnts); 

saveas(gcf,['mod_regress' num2str(nfig)],'epsc');saveas(gcf,['mod_regress' num2str(nfig)],'jpg');savefig(['mod_regress' num2str(nfig) '.fig'])

figure(); nfig=2;
set(gcf,'units','points','position',[100,100,1000,600],'color','w','defaultTextInterpreter','latex','defaultLegendInterpreter','latex');
x_t = linspace(-bnd,bnd,100);
hold on;plot(x_t,squeeze(twoPL(x_t,1,min(b))),'LineWidth',1.5);plot(x_t,squeeze(twoPL(x_t,1,mean(b))),'LineWidth',1.5);plot(x_t,squeeze(twoPL(x_t,1,max(b))),'LineWidth',1.5);
vline(0,'--k');hline(0.5,'--k');legend({'Low bigram','Medium bigram','High bigram'},'Location','best','FontSize',fntl);legend('boxoff')
set(gca,'FontSize',fntl,'FontName','Times');title('Probability of distractor','FontSize',fntst); xlabel('z','FontSize',fnts); ylabel('$\pi$','FontSize',fnts); 
text(1,0.2,['p(distractor$|$low bigram) = ' num2str(round(twoPL(0,1,min(b))-twoPL(0,1,5),2))],'FontSize',fntl)
text(1,0.15,['p(distractor$|$med bigram) = ' num2str(round(twoPL(0,1,mean(b))-twoPL(0,1,5),2),2)],'FontSize',fntl)
text(1,0.1,['p(distractor$|$high bigram) = ' num2str(round(twoPL(0,1,max(b))-twoPL(0,1,5),2))],'FontSize',fntl)
saveas(gcf,'/home/antonio/MEGA/Lavoro_sync/My papers/Draft/A state space approach to dynamic modeling of mouse-tracking data/JMP/mod_regress_2','epsc');
savefig('/home/antonio/MEGA/Lavoro_sync/My papers/Draft/A state space approach to dynamic modeling of mouse-tracking data/JMP/mod_regress_2.fig')

saveas(gcf,['mod_regress' num2str(nfig)],'epsc');saveas(gcf,['mod_regress' num2str(nfig)],'jpg');savefig(['mod_regress' num2str(nfig) '.fig'])

figure(); nfig=3;
iid_w = [find(b==min(b)) find(round(b,2)==0) find(b==max(b))];
CLS = get(gca,'colororder');close all;
set(gcf,'units','points','position',[100,100,1500,900],'color','w','defaultTextInterpreter','latex','defaultLegendInterpreter','latex');
for i=1:I
    subplot(4,6,i);for j=1:length(iid_w),hold on;plot(squeeze(log(P(i,iid_w(j),:)./(1-P(i,iid_w(j),:)))),'Color',CLS(j,:));end
    if (i==1||i==7|i==13|i==19),ylabel('$\log(\pi_D/\pi_T)$','FontSize',14);end
    xlim([0 T+2]);ylim([-7 7]);set(gca,'FontSize',14,'FontName','Times'); hline(0,'--k')
    xticks([0 35 71]);xticklabels({'0%','50%','100%'}); title(['Subject no.: ' num2str(i)],'FontSize',15); yticks([-5 0 5]);yticklabels({'T','0','D'});         
end
legend({'Low bigram','Medium bigram','High bigram'},'Location','none','Position',[0.75 0.15 0 0],'FontSize',14,'Units', 'normalized');legend('boxoff')

saveas(gcf,['mod_regress' num2str(nfig)],'epsc');saveas(gcf,['mod_regress' num2str(nfig)],'jpg');savefig(['mod_regress' num2str(nfig) '.fig'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% evaluate beta in terms of decomposed coefficients
x_t = linspace(-bnd,bnd,100);
subplot(1,2,1);hold on;plot(x_t,squeeze(twoPL(x_t,1,min(b))),'--r');plot(x_t,squeeze(twoPL(x_t,1,0)),'k');plot(x_t,squeeze(twoPL(x_t,1,max(b))),'--r');vline(0,'--k');hline(0.5,'--k')
subplot(1,2,2);hold on;plot(x_t,information_function(x_t,min(b),1),'--r');plot(x_t,information_function(x_t,0,1),'k');plot(x_t,information_function(x_t,max(b),1),'--r');vline(0,'--k')

csi1 = (twoPL(0,1,min(b))-twoPL(0,1,0))/0.5 %probability to select D associated to low covariate
csi2 = (twoPL(0,1,max(b))-twoPL(0,1,0))/0.5 %probability to select D associated to high covariate

% individual profiles (Z)
figure(); 
for i=1:I, subplot(4,6,i); 
    hold on; plot(XF(i,:),'b'); %plot(XF(i,:)-PF(i,:),'--k'); plot(XF(i,:)+PF(i,:),'--k');
    xlim([-0.5 T+3]); ylim([-bnd-2 bnd+2]); hline(0,'--k'); title(['Sbj ' num2str(i)]); end

% prob profiles (PI)
iid_w = [find(b==min(b)) find(round(b,2)==0) find(b==max(b))];
figure();
cls = {'r','k','b'};
for i=1:I
    k=0;
    for j=1:length(iid_w)
        k=k+1;
        subplot(4,6,i);hold on; plot(squeeze(P(i,iid_w(j),:)),cls{k}); 
        xlim([-0.5 T+3]); ylim([-0.1 1.1]); hline(0.5,'--k'); title(['Sbj ' num2str(i)])
    end
end
legend('Low-bigramF','0-bigramF','High-bigramF')

% prob profiles (PI) clustered by covariate levels
iid_w = [find(b==min(b)) find(round(b,2)==0) find(b==max(b))];
figure();
cls = {'Low-bigramF','0-bigramF','High-bigramF'};
for j=1:length(iid_w)
    subplot(1,3,j);hold on
    for i=1:I
        plot(squeeze(P(i,iid_w(j),:))); 
        xlim([-0.5 T+3]); ylim([-0.1 1.1]); hline(0.5,'--k'); 
    end
    title(cls{j})
end




%% Results: model 5_1 (categorical)
clear all
information_function = @(xdom,gamma_est,a) cell2mat(arrayfun(@(j)exp(bsxfun(@times,bsxfun(@minus,xdom,gamma_est(j)),-a(j)))./...
    (1+exp(bsxfun(@times,bsxfun(@minus,xdom,gamma_est(j)),-a(j)))).^2,1:length(gamma_est),'UniformOutput',false)');

load('barca2012_nomissin_full.mat');load('mh_5_1.mat')

q0=1;q=72; %select just {HF,LF,NW}
data.Z=Z(:,q0:q,:); data.Y=Y(:,q0:q,:); data.D1=D1(q0:q,1:3); pars.sigmax=sigmax; pars.kappa1=kappa1; pars.kappa2=kappa2; pars.a=a(q0:q); 
pars.mu1=mu1; pars.mu2=mu2; pars.bnd=bnd; data.x=scaledata(bigF_cov(q0:q),-3,3);

[S_coll, post_means,Rgelm,neff,S] = MCMC_diag(mh_5_1,1000,1,false,true,'categorical',3);
b = data.D1*post_means';
tic;[XF,PF] = GaussApprox_filter(data.Z,pars.sigmax,pars.a,b,0);toc 
P = twoPL(XF,pars.a,b);

%%%%%%%%%%%%% Figure for paper %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mh=mh_5_1; Rgelm=round(Rgelm,3); mms = round(mean(S_coll),2); fnts=20; fntst=22; fntl=14; 

figure(); nfig=1; tls={'HF','LF','NW'};
set(gcf,'units','points','position',[100,100,1500,800],'color','w','defaultTextInterpreter','latex','defaultLegendInterpreter','latex');

k=0;
for j=1:3
    k=k+1;
    subplot(3,3,k);hold on;for i=1:size(S,3),plot(S(:,k,i));end; xlim([-250 size(S,1)+250]); 
    set(gca,'FontSize',fntl,'FontName','Times');title(['$\gamma_' num2str(j) '$ ' '(' tls{j} ')'],'FontSize',fntst); xlabel(['Iterations ($\hat R = ' num2str(Rgelm(j)) '$)'],'FontSize',fnts); 
    if k==1,ylabel('Traceplot','FontSize',fnts);end
end

for j=1:3
    k=k+1;
    [f,xi] = ksdensity(S_coll(:,j)); xhdi = mbe_hdi(S_coll(:,j),.95); iid=min(find(round(xi,2)==round(xhdi(1),2))):min(find(round(xi,2)==round(xhdi(2),2)))
    subplot(3,3,k);plot(xi,f,'LineWidth',1.5);box('off'); hold on; area(xi,f,'FaceColor',[0.5843 0.8157 0.9882]); plot(xi(iid),f(iid));area(xi(iid),f(iid),'FaceColor',[0 0.4470 0.7410]);
    text(median(xi),quantile(f,0.5),'95% HDPI','FontSize',13,'HorizontalAlignment','center')
    set(gca,'FontSize',fntl,'FontName','Times'); xlabel(['$\mu = ' num2str(mms(j)) '$'],'FontSize',fnts); if k==4,ylabel('Marginal posterior density','FontSize',fnts);end
    box('off')    
end

for j=1:3
    k=k+1;
    subplot(3,3,k);hold on; [acfs,lags] = autocorr(S_coll(:,j)); area([-0.2 max(lags)+1],[0.25 0.25],'EdgeColor','w','FaceColor',rgb('lightgray')); area([-0.2 max(lags)+1],[-0.25 -0.25],'EdgeColor','w','FaceColor',rgb('lightgray'));
    bar(lags,acfs,0.1,'FaceColor',[0 0.4470 0.7410]); plot(lags,acfs,'o','Color',[0 0.4470 0.7410],'MarkerFaceColor',[0 0.4470 0.7410]); ylim([-1 1]); xlim([-0.5 max(lags)+1]); box('off'); 
    set(gca,'FontSize',fntl,'FontName','Times'); xlabel('lags','FontSize',fnts); if k==7,ylabel('ACF','FontSize',fnts); end
end

saveas(gcf,['mod_categ' num2str(nfig)],'epsc');saveas(gcf,['mod_categ' num2str(nfig)],'jpg');savefig(['mod_categ' num2str(nfig) '.fig'])

figure(); nfig=2;
set(gcf,'units','points','position',[100,100,1000,600],'color','w','defaultTextInterpreter','latex','defaultLegendInterpreter','latex');
x_t = linspace(-bnd,bnd,100);
hold on;plot(x_t,squeeze(twoPL(x_t,1,post_means(1))));plot(x_t,squeeze(twoPL(x_t,1,post_means(2))),'Color',[0.8500    0.3250    0.0980]);plot(x_t,squeeze(twoPL(x_t,1,post_means(3))),'Color',[0.4660    0.6740    0.1880]);
vline(0,'--k');hline(0.5,'--k');legend({'HF','LF','NW'},'Location','northeastoutside','FontSize',fntl);legend('boxoff')
set(gca,'FontSize',fntl,'FontName','Times');title('Probability of the distractor cue','FontSize',fntst); xlabel('z','FontSize',fnts); ylabel('$\pi$','FontSize',fnts); 
text(5,0.4,['p(distractor$|$HF) = ' num2str(round(twoPL(0,1,post_means(1))-twoPL(0,1,5),2))],'FontSize',fntl)
text(5,0.35,['p(distractor$|$LF) = ' num2str(round(twoPL(0,1,post_means(2))-twoPL(0,1,5),2),2)],'FontSize',fntl)
text(5,0.3,['p(distractor$|$NW) = ' num2str(round(twoPL(0,1,post_means(3))-twoPL(0,1,5),2))],'FontSize',fntl)

saveas(gcf,['mod_categ' num2str(nfig)],'epsc');saveas(gcf,['mod_categ' num2str(nfig)],'jpg');savefig(['mod_categ' num2str(nfig) '.fig'])


figure(); nfig=3;
iid_w = [min(find(D1(:,1)==1)) min(find(D1(:,2)==1)) min(find(D1(:,3)==1))];
CLS = get(gca,'colororder');close all;
set(gcf,'units','points','position',[100,100,1500,900],'color','w','defaultTextInterpreter','latex','defaultLegendInterpreter','latex');
for i=1:I
    subplot(4,6,i);for j=1:length(iid_w),hold on;plot([0;squeeze(log(P(i,iid_w(j),:)./(1-P(i,iid_w(j),:))))],'Color',CLS(j,:));end
    if (i==1||i==7|i==13|i==19),ylabel('$\log(\pi_D/\pi_T)$','FontSize',14);end
    xlim([0 T+2]);ylim([-7 7]);set(gca,'FontSize',14,'FontName','Times'); hline(0,'--k')
    xticks([0 35 71]);xticklabels({'0%','50%','100%'}); title(['Subject no.: ' num2str(i)],'FontSize',15); yticks([-5 0 5]);yticklabels({'T','0','D'});         
end
legend({'HF','LF','NW'},'Location','none','Position',[0.75 0.15 0 0],'FontSize',14,'Units', 'normalized');legend('boxoff')

saveas(gcf,['mod_categ' num2str(nfig)],'epsc');saveas(gcf,['mod_categ' num2str(nfig)],'jpg');savefig(['mod_categ' num2str(nfig) '.fig'])

% figure(); nfig=4;
% CLS = get(gca,'colororder');close all;
% set(gcf,'units','points','position',[100,100,1200,450],'color','w','defaultTextInterpreter','latex','defaultLegendInterpreter','latex');
% subplot(1,3,1);hold on;ksdensity(S_coll(:,1));ksdensity(S_coll(:,2)); legend({'HF','LF'},'Location','northeastoutside','FontSize',13);legend('boxoff');ylabel('Marginal posterior density','FontSize',14)
% subplot(1,3,2);hold on;ksdensity(S_coll(:,1));ksdensity(S_coll(:,3)); legend({'HF','NW'},'Location','northeastoutside','FontSize',13);legend('boxoff')
% subplot(1,3,3);hold on;ksdensity(S_coll(:,2));ksdensity(S_coll(:,3)); legend({'LF','NW'},'Location','northeastoutside','FontSize',13);legend('boxoff')

%%% => new figures (JMP version)
figure(); nfig=4;
CLS = get(gca,'colororder');close all;
set(gcf,'units','points','position',[100,100,1200,400],'color','w','defaultTextInterpreter','latex','defaultLegendInterpreter','latex');
subplot(1,2,1);hold on;
set(gca,'FontSize',fntl,'FontName','Times');title('$\mathbf{(A)}$ Marginal posterior density','FontSize',fntst); xlabel('$\theta$','FontSize',fnts); ylabel('$f(\theta|Y)$','FontSize',fnts); 
[xi,fi]=ksdensity(S_coll(:,1));plot(fi,xi,'LineWidth',1.5);
[xi,fi]=ksdensity(S_coll(:,2));plot(fi,xi,'LineWidth',1.5);
[xi,fi]=ksdensity(S_coll(:,3));plot(fi,xi,'LineWidth',1.5);  
legend({'$\gamma_1 $ (HF)','$\gamma_2$ (LF)','$\gamma_3$ (NW)'},'Location','northeast','FontSize',13);legend('boxoff');%ylabel('Marginal posterior density','FontSize',14)
subplot(1,2,2);hold on;
set(gca,'FontSize',fntl,'FontName','Times');title('$\mathbf{(B)}$ Probability of distractor','FontSize',fntst); xlabel('z','FontSize',fnts); ylabel('$\pi$','FontSize',fnts); 
plot(x_t,squeeze(twoPL(x_t,1,post_means(1))),'LineWidth',1.5);plot(x_t,squeeze(twoPL(x_t,1,post_means(2))),'LineWidth',1.5);plot(x_t,squeeze(twoPL(x_t,1,post_means(3))),'LineWidth',1.5);
vline(0,'--k');hline(0.5,'--k');legend({'HF','LF','NW'},'Location','northeastoutside','FontSize',fntl);legend('boxoff')
text(1,0.22,['p(distractor$|$HF) = ' num2str(round(twoPL(0,1,post_means(1))-twoPL(0,1,5),2))],'FontSize',fntl)
text(1,0.15,['p(distractor$|$LF) = ' num2str(round(twoPL(0,1,post_means(2))-twoPL(0,1,5),2),2)],'FontSize',fntl)
text(1,0.08,['p(distractor$|$NW) = ' num2str(round(twoPL(0,1,post_means(3))-twoPL(0,1,5),2))],'FontSize',fntl)

saveas(gcf,'/home/antonio/MEGA/Lavoro_sync/My papers/Draft/A state space approach to dynamic modeling of mouse-tracking data/JMP/mod_categ_2','epsc');
savefig('/home/antonio/MEGA/Lavoro_sync/My papers/Draft/A state space approach to dynamic modeling of mouse-tracking data/JMP/mod_categ_2.fig')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Evaluate posteriors
figure()
subplot(1,3,1);hold on;ksdensity(S_coll(:,1));ksdensity(S_coll(:,2)); title('HF vs. LF');
subplot(1,3,2);hold on;ksdensity(S_coll(:,1));ksdensity(S_coll(:,3)); title('HF vs. NW'); 
subplot(1,3,3);hold on;ksdensity(S_coll(:,2));ksdensity(S_coll(:,3)); title('LF vs. NW'); 

% evaluate beta in terms of decomposed coefficients
x_t = linspace(-bnd,bnd,100); 
figure(); cls={'r','b','g'};
plot(x_t,squeeze(twoPL(x_t,1,0)),'--k');
for k=1:3
    hold on; plot(x_t,squeeze(twoPL(x_t,1,post_means(k)))',cls{k});
end
legend('0','HF','LF','NW')

csi1 = (twoPL(0,1,post_means(1))-twoPL(0,1,0))/0.5 %probability to select D associated to HF
csi2 = (twoPL(0,1,post_means(2))-twoPL(0,1,0))/0.5 %probability to select D associated to LF
csi3 = (twoPL(0,1,post_means(3))-twoPL(0,1,0))/0.5 %probability to select D associated to NW

% prob profiles (PI)
iid_w = [min(find(D1(:,1)==1)) min(find(D1(:,2)==1)) min(find(D1(:,3)==1))];
figure();
cls = {'r','k','b'};
for i=1:I
    k=0;
    for j=1:length(iid_w)
        k=k+1;
        subplot(4,6,i);hold on; plot(squeeze(P(i,iid_w(j),:)),cls{k}); 
        xlim([-0.5 T+3]); ylim([-0.1 1.1]); hline(0.5,'--k'); title(['Sbj ' num2str(i)])
    end
end
legend('HF','LF','NW')


% prob profiles (PI) clustered by factor levels
iid_w = [min(find(D1(:,1)==1)) min(find(D1(:,2)==1)) min(find(D1(:,3)==1))];
figure();
cls = {'HF','LF','NW'};
for j=1:length(iid_w)
    subplot(1,3,j);hold on
    for i=1:I
        plot(squeeze(P(i,iid_w(j),:))); 
        xlim([-0.5 T+3]); ylim([-0.1 1.1]); hline(0.5,'--k'); 
    end
    title(cls{j})
end

%% Results: model 5_2 (interaction)
clear all
information_function = @(xdom,gamma_est,a) cell2mat(arrayfun(@(j)exp(bsxfun(@times,bsxfun(@minus,xdom,gamma_est(j)),-a(j)))./...
    (1+exp(bsxfun(@times,bsxfun(@minus,xdom,gamma_est(j)),-a(j)))).^2,1:length(gamma_est),'UniformOutput',false)');

load('barca2012_nomissin_full.mat');load('mh_5_2.mat')

q0=1;q=72; %select just {HF,LF,NW}
data.Z=Z(:,q0:q,:); data.Y=Y(:,q0:q,:); data.D1=D1(q0:q,1:3); pars.sigmax=sigmax; pars.kappa1=kappa1; pars.kappa2=kappa2; pars.a=a(q0:q); 
pars.mu1=mu1; pars.mu2=mu2; pars.bnd=bnd; data.x=scaledata(bigF_cov(q0:q),-3,3);

[S_coll, post_means,Rgelm,~,S] = MCMC_diag(mh_5_2,2250,3,false,true,'interaction',3);
[b,gamma1,eta,delta] = switchModel('interaction',data.D1,data.x,post_means');

tic;[XF,PF] = GaussApprox_filter(data.Z,pars.sigmax,pars.a,b,0);toc 
P = twoPL(XF,pars.a,b);

%%%%%%%%%%%%% Figure for paper %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mh=mh_5_2; Rgelm=round(Rgelm,3); mms = round(mean(S_coll),2); fnts=20; fntst=22; fntl=12; 

figure(); nfig=1;
set(gcf,'units','points','position',[100,100,1500,800],'color','w','defaultTextInterpreter','latex','defaultLegendInterpreter','latex');

k=0;  tls={'HF','LF','NW','bigram x HF','bigram x LF','bigram x NW'}; tls2={'$\gamma_1$','$\gamma_2$','$\gamma_3$','$\eta$','$\delta\_1$','$\delta\_2$'};
for j=1:6
    k=k+1;
    subplot(3,6,k);hold on;for i=1:size(S,3),plot(S(:,k,i));end; xlim([-50 size(S,1)+50]); 
    set(gca,'FontSize',11,'FontName','Times');title(['$' tls2{j} '$ ' '(' tls{j} ')'],'FontSize',16); xlabel(['Iterations ($\hat R = ' num2str(Rgelm(j)) '$)'],'FontSize',14); 
    if k==1,ylabel('Traceplot','FontSize',13);end
end

j=1;k=k+1;
[f,xi] = ksdensity(S_coll(:,j)); xhdi = mbe_hdi(S_coll(:,j),.95); iid=32:70;
subplot(3,6,k);plot(xi,f,'LineWidth',1.5);box('off'); hold on; area(xi,f,'FaceColor',[0.5843 0.8157 0.9882]); plot(xi(iid),f(iid));area(xi(iid),f(iid),'FaceColor',[0 0.4470 0.7410]);
%text(median(xi),quantile(f,0.5),'95% HDPI','FontSize',8,'HorizontalAlignment','center')
set(gca,'FontSize',11,'FontName','Times'); xlabel(['$\mu = ' num2str(mms(j)) '$'],'FontSize',14); ylabel('Marginal posterior density','FontSize',14);
box('off')    
for j=2:6
    k=k+1;
    [f,xi] = ksdensity(S_coll(:,j)); xhdi = mbe_hdi(S_coll(:,j),.95); iid=min(find(round(xi,2)==round(xhdi(1),1))):min(find(round(xi,2)==round(xhdi(2),1)));
    subplot(3,6,k);plot(xi,f,'LineWidth',1.5);box('off'); hold on; area(xi,f,'FaceColor',[0.5843 0.8157 0.9882]); plot(xi(iid),f(iid));area(xi(iid),f(iid),'FaceColor',[0 0.4470 0.7410]);
    %text(median(xi),quantile(f,0.5),'95% HDPI','FontSize',8,'HorizontalAlignment','center')
    set(gca,'FontSize',11,'FontName','Times'); xlabel(['$\mu = ' num2str(mms(j)) '$'],'FontSize',14);
    box('off')    
end

for j=1:6
    k=k+1;
    subplot(3,6,k);hold on; [acfs,lags] = autocorr(S_coll(:,j)); area([-0.2 max(lags)+1],[0.25 0.25],'EdgeColor','w','FaceColor',rgb('lightgray')); area([-0.2 max(lags)+1],[-0.25 -0.25],'EdgeColor','w','FaceColor',rgb('lightgray'));
    bar(lags,acfs,0.1,'FaceColor',[0 0.4470 0.7410]); plot(lags,acfs,'o','Color',[0 0.4470 0.7410],'MarkerFaceColor',[0 0.4470 0.7410]); ylim([-1 1]); xlim([-0.5 max(lags)+1]); box('off'); 
    set(gca,'FontSize',11,'FontName','Times'); xlabel('lags','FontSize',14); if k==13,ylabel('ACF','FontSize',13); end
end
saveas(gcf,'/home/antonio/MEGA/Lavoro_sync/My papers/Draft/A state space approach to dynamic modeling of mouse-tracking data/JMP/mod_interact_1','epsc')
savefig('/home/antonio/MEGA/Lavoro_sync/My papers/Draft/A state space approach to dynamic modeling of mouse-tracking data/JMP/mod_interact_1.fig')


saveas(gcf,['mod_interact' num2str(nfig)],'epsc');saveas(gcf,['mod_interact' num2str(nfig)],'jpg');savefig(['mod_interact' num2str(nfig) '.fig'])

figure(); nfig=2;
set(gcf,'units','points','position',[100,100,1500,500],'color','w','defaultTextInterpreter','latex','defaultLegendInterpreter','latex');
subplot(1,3,1);hold on; plot(x_t,squeeze(twoPL(x_t,1, gamma1(1)+min(data.x)*(eta+delta(1)) ))',cls{1}); plot(x_t,squeeze(twoPL(x_t,1, gamma1(2)+min(data.x)*(eta+delta(2)) ))',cls{2}); plot(x_t,squeeze(twoPL(x_t,1, gamma1(3)+min(data.x)*(eta+delta(3)) ))',cls{3}); 
set(gca,'FontSize',fntl,'FontName','Times');title('Low bigram','FontSize',fntst); xlabel('z','FontSize',fnts); ylabel('$\pi$','FontSize',fnts); 
vline(0,'--k');hline(0.5,'--k');
text(1,0.2,['p(distractor$|$HF) = ' num2str(round(twoPL(0,1,gamma1(1)+min(data.x)*(eta+delta(1)))-twoPL(0,1,5),2))],'FontSize',13)
text(1,0.15,['p(distractor$|$LF) = ' num2str(round(twoPL(0,1,gamma1(2)+min(data.x)*(eta+delta(2)))-twoPL(0,1,5),2),2)],'FontSize',13)
text(1,0.1,['p(distractor$|$NW) = ' num2str(round(twoPL(0,1,gamma1(3)+min(data.x)*(eta+delta(3)))-twoPL(0,1,5),2))],'FontSize',13)

subplot(1,3,2);hold on; plot(x_t,squeeze(twoPL(x_t,1, gamma1(1)+0.01*(eta+delta(1)) ))',cls{1}); plot(x_t,squeeze(twoPL(x_t,1, gamma1(2)+0.01*(eta+delta(2)) ))',cls{2}); plot(x_t,squeeze(twoPL(x_t,1, gamma1(3)+0.01*(eta+delta(3)) ))',cls{3}); title('0-bigramF')
set(gca,'FontSize',fntl,'FontName','Times');title('Medium bigram','FontSize',fntst); xlabel('z','FontSize',fnts); ylabel('$\pi$','FontSize',fnts); 
vline(0,'--k');hline(0.5,'--k');
text(1,0.2,['p(distractor$|$HF) = ' num2str(round(twoPL(0,1,gamma1(1)+mean(data.x)*(eta+delta(1)))-twoPL(0,1,5),2))],'FontSize',13)
text(1,0.15,['p(distractor$|$LF) = ' num2str(round(twoPL(0,1,gamma1(2)+mean(data.x)*(eta+delta(2)))-twoPL(0,1,5),2),2)],'FontSize',13)
text(1,0.1,['p(distractor$|$NW) = ' num2str(round(twoPL(0,1,gamma1(3)+mean(data.x)*(eta+delta(3)))-twoPL(0,1,5),2))],'FontSize',13)

subplot(1,3,3);hold on; plot(x_t,squeeze(twoPL(x_t,1, gamma1(1)+max(data.x)*(eta+delta(1)) ))',cls{1}); plot(x_t,squeeze(twoPL(x_t,1, gamma1(2)+max(data.x)*(eta+delta(2)) ))',cls{2}); plot(x_t,squeeze(twoPL(x_t,1, gamma1(3)+max(data.x)*(eta+delta(3)) ))',cls{3});title('High-bigramF')
set(gca,'FontSize',fntl,'FontName','Times');title('High bigram','FontSize',fntst); xlabel('z','FontSize',fnts); ylabel('$\pi$','FontSize',fnts); 
vline(0,'--k');hline(0.5,'--k');
legend({'HF','LF','NW'},'Location','northeastoutside','FontSize',fntl);legend('boxoff')
text(1,0.2,['p(distractor$|$HF) = ' num2str(round(twoPL(0,1,gamma1(1)+max(data.x)*(eta+delta(1)))-twoPL(0,1,5),2))],'FontSize',13)
text(1,0.15,['p(distractor$|$LF) = ' num2str(round(twoPL(0,1,gamma1(2)+max(data.x)*(eta+delta(2)))-twoPL(0,1,5),2),2)],'FontSize',13)
text(1,0.1,['p(distractor$|$NW) = ' num2str(round(twoPL(0,1,gamma1(3)+max(data.x)*(eta+delta(3)))-twoPL(0,1,5),2))],'FontSize',13)

saveas(gcf,['mod_interact' num2str(nfig)],'epsc');saveas(gcf,['mod_interact' num2str(nfig)],'jpg');savefig(['mod_interact' num2str(nfig) '.fig'])

figure(); nfig=4;
CLS = get(gca,'colororder');close all;
set(gcf,'units','points','position',[100,100,1400,550],'color','w','defaultTextInterpreter','latex','defaultLegendInterpreter','latex');
subplot(2,3,1);hold on;ksdensity(S_coll(:,1));ksdensity(S_coll(:,2)); legend({'HF','LF'},'Location','northeastoutside','FontSize',13);legend('boxoff');ylabel('Marginal posterior density','FontSize',14)
subplot(2,3,2);hold on;ksdensity(S_coll(:,1));ksdensity(S_coll(:,3)); legend({'HF','NF'},'Location','northeastoutside','FontSize',13);legend('boxoff');
subplot(2,3,3);hold on;ksdensity(S_coll(:,2));ksdensity(S_coll(:,3)); legend({'NW','LF'},'Location','northeastoutside','FontSize',13);legend('boxoff'); 
subplot(2,3,4);hold on;ksdensity(S_coll(:,4));ksdensity(S_coll(:,5)); legend({'HF:bigram','LF:bigram'},'Location','northeastoutside','FontSize',13);legend('boxoff');ylabel('Marginal posterior density','FontSize',14)
subplot(2,3,5);hold on;ksdensity(S_coll(:,4));ksdensity(S_coll(:,6)); legend({'HF:bigram','NF:bigram'},'Location','northeastoutside','FontSize',13);legend('boxoff');
subplot(2,3,6);hold on;ksdensity(S_coll(:,5));ksdensity(S_coll(:,6)); legend({'LF:bigram','NW:bigram'},'Location','northeastoutside','FontSize',13);legend('boxoff');


%%% => new figures (JMP)
figure(); nfig=4;
CLS = get(gca,'colororder');close all;
set(gcf,'units','points','position',[100,100,1200,400],'color','w','defaultTextInterpreter','latex','defaultLegendInterpreter','latex');
subplot(1,2,1);hold on;
set(gca,'FontSize',fntl,'FontName','Times');%title('$\mathbf{(A)}$ Marginal posterior density','FontSize',fntst); 
xlabel('$\theta$','FontSize',fnts); ylabel('$f(\theta|Y)$','FontSize',fnts); 
[xi,fi]=ksdensity(S_coll(:,1));plot(fi,xi,'LineWidth',1.5);
[xi,fi]=ksdensity(S_coll(:,2));plot(fi,xi,'LineWidth',1.5);
[xi,fi]=ksdensity(S_coll(:,3));plot(fi,xi,'LineWidth',1.5);
legend({'$\gamma_1$ (HF)','$\gamma_2$ (LF)','$\gamma_3$ (NW)'},'Location','northeastoutside','FontSize',13);legend('boxoff'); 
subplot(1,2,2);hold on;
set(gca,'FontSize',fntl,'FontName','Times');%title('$\mathbf{(B)}$ Marginal posterior density','FontSize',fntst); 
xlabel('$\theta$','FontSize',fnts); ylabel('$f(\theta|Y)$','FontSize',fnts); 
[xi,fi]=ksdensity(S_coll(:,4));plot(fi,xi,'LineWidth',1.5,'Color',CLS(4,:));
[xi,fi]=ksdensity(S_coll(:,5));plot(fi,xi,'LineWidth',1.5,'Color',CLS(5,:));
[xi,fi]=ksdensity(S_coll(:,6));plot(fi,xi,'LineWidth',1.5,'Color',CLS(6,:));
legend({'$\eta$ (bigram x HF)','$\delta_1$ (bigram x LF)','$\delta_2$ (bigram x NW)'},'Location','northeastoutside','FontSize',13);legend('boxoff'); 
saveas(gcf,'/home/antonio/MEGA/Lavoro_sync/My papers/Draft/A state space approach to dynamic modeling of mouse-tracking data/JMP/mod_interact_2','epsc')
savefig('/home/antonio/MEGA/Lavoro_sync/My papers/Draft/A state space approach to dynamic modeling of mouse-tracking data/JMP/mod_interact_2.fig')

figure();
set(gcf,'units','points','position',[100,100,1500,500],'color','w','defaultTextInterpreter','latex','defaultLegendInterpreter','latex');
subplot(1,3,1);hold on;
plot(x_t,squeeze(twoPL(x_t,1, gamma1(1)+min(data.x)*(eta+delta(1)) ))','Color',CLS(1,:),'LineWidth',1.5); 
plot(x_t,squeeze(twoPL(x_t,1, gamma1(2)+min(data.x)*(eta+delta(2)) ))','Color',CLS(2,:),'LineWidth',1.5);
plot(x_t,squeeze(twoPL(x_t,1, gamma1(3)+min(data.x)*(eta+delta(3)) ))','Color',CLS(3,:),'LineWidth',1.5); 
set(gca,'FontSize',fntl,'FontName','Times');title('Low bigram','FontSize',fntst); xlabel('z','FontSize',fnts); ylabel('$\pi$','FontSize',fnts); 
vline(0,'--k');hline(0.5,'--k');
text(1,0.2,['p(distractor$|$HF) = ' num2str(round(twoPL(0,1,gamma1(1)+min(data.x)*(eta+delta(1)))-twoPL(0,1,5),2))],'FontSize',13)
text(1,0.15,['p(distractor$|$LF) = ' num2str(round(twoPL(0,1,gamma1(2)+min(data.x)*(eta+delta(2)))-twoPL(0,1,5),2),2)],'FontSize',13)
text(1,0.1,['p(distractor$|$NW) = ' num2str(round(twoPL(0,1,gamma1(3)+min(data.x)*(eta+delta(3)))-twoPL(0,1,5),2))],'FontSize',13)

subplot(1,3,2);hold on; 
plot(x_t,squeeze(twoPL(x_t,1, gamma1(1)+0.01*(eta+delta(1)) ))','Color',CLS(1,:),'LineWidth',1.5); 
plot(x_t,squeeze(twoPL(x_t,1, gamma1(2)+0.01*(eta+delta(2)) ))','Color',CLS(2,:),'LineWidth',1.5); 
plot(x_t,squeeze(twoPL(x_t,1, gamma1(3)+0.01*(eta+delta(3)) ))','Color',CLS(3,:),'LineWidth',1.5); 
set(gca,'FontSize',fntl,'FontName','Times');title('Medium bigram','FontSize',fntst); xlabel('z','FontSize',fnts); ylabel('$\pi$','FontSize',fnts); 
vline(0,'--k');hline(0.5,'--k');
text(1,0.2,['p(distractor$|$HF) = ' num2str(round(twoPL(0,1,gamma1(1)+mean(data.x)*(eta+delta(1)))-twoPL(0,1,5),2))],'FontSize',13)
text(1,0.15,['p(distractor$|$LF) = ' num2str(round(twoPL(0,1,gamma1(2)+mean(data.x)*(eta+delta(2)))-twoPL(0,1,5),2),2)],'FontSize',13)
text(1,0.1,['p(distractor$|$NW) = ' num2str(round(twoPL(0,1,gamma1(3)+mean(data.x)*(eta+delta(3)))-twoPL(0,1,5),2))],'FontSize',13)

subplot(1,3,3);hold on; 
plot(x_t,squeeze(twoPL(x_t,1, gamma1(1)+max(data.x)*(eta+delta(1)) ))','Color',CLS(1,:),'LineWidth',1.5); 
plot(x_t,squeeze(twoPL(x_t,1, gamma1(2)+max(data.x)*(eta+delta(2)) ))','Color',CLS(2,:),'LineWidth',1.5); 
plot(x_t,squeeze(twoPL(x_t,1, gamma1(3)+max(data.x)*(eta+delta(3)) ))','Color',CLS(3,:),'LineWidth',1.5);
set(gca,'FontSize',fntl,'FontName','Times');title('High bigram','FontSize',fntst); xlabel('z','FontSize',fnts); ylabel('$\pi$','FontSize',fnts); 
vline(0,'--k');hline(0.5,'--k');
legend({'HF','LF','NW'},'Location','northeastoutside','FontSize',fntl);legend('boxoff')
text(1,0.2,['p(distractor$|$HF) = ' num2str(round(twoPL(0,1,gamma1(1)+max(data.x)*(eta+delta(1)))-twoPL(0,1,5),2))],'FontSize',13)
text(1,0.15,['p(distractor$|$LF) = ' num2str(round(twoPL(0,1,gamma1(2)+max(data.x)*(eta+delta(2)))-twoPL(0,1,5),2),2)],'FontSize',13)
text(1,0.1,['p(distractor$|$NW) = ' num2str(round(twoPL(0,1,gamma1(3)+max(data.x)*(eta+delta(3)))-twoPL(0,1,5),2))],'FontSize',13)
saveas(gcf,'/home/antonio/MEGA/Lavoro_sync/My papers/Draft/A state space approach to dynamic modeling of mouse-tracking data/JMP/mod_interact_3','epsc')
savefig('/home/antonio/MEGA/Lavoro_sync/My papers/Draft/A state space approach to dynamic modeling of mouse-tracking data/JMP/mod_interact_3.fig')

saveas(gcf,['mod_interact' num2str(nfig)],'epsc');saveas(gcf,['mod_interact' num2str(nfig)],'jpg');savefig(['mod_interact' num2str(nfig) '.fig'])

savefig(['mod_interact' num2str(nfig) '.fig'])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% ==> new figures JMP
for i=1:I
        for j=1:size(P,2)
            for t=2:T                
                Ystar(i,j,t) = rmixedvm(1,mu1,mu2,kappa1,kappa2,P(i,j,t));                
            end
        end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% INDIVIDUAL ANALYSIS: 1
[b,gamma1,eta,delta] = switchModel('interaction',data.D1,data.x,post_means');
tic;[XF] = GaussApprox_filter(data.Z,pars.sigmax,pars.a,b,0);toc 
P = twoPL(XF,pars.a,b);

g1=[6,7,8,15]; g2=[1,4,19,21]; g3=[2,3,5,12,13,16,17,20,22]; g4=[10,11,14,18];
% CLS=[0 0 0.8; 0.8 0 0;0 0.8 0; 0 0 0.17];
% CLS = CLS + unifrnd(0,0.1,4,3); CLS(CLS>1)=1;
CLS = colormap(lines(4));
lns=1.25;

set(gcf,'units','points','position',[100,100,1200,700],'color','w','defaultTextInterpreter','latex','defaultLegendInterpreter','latex');

subplot(2,4,[1:2,5:6])
hold on; 
for i=1:length(g1)
    plot(XF(g1(i),:),'color',CLS(1,:),'LineWidth',lns);    
end
ylabel('$Z$','FontSize',14);xlim([0 T+2]);ylim([-7 7]);set(gca,'FontSize',14,'FontName','Times'); hline(0,'--k')
xticks([0 35 71]);xticklabels({'0%','50%','100%'}); yticks([-5 0 5]);yticklabels({'T','0','D'});box('off')         

hold on;
for i=1:length(g2)
    plot(XF(g2(i),:),'color',CLS(2,:),'LineWidth',lns);
end

hold on; 
for i=1:length(g3)
    plot(XF(g3(i),:),'color',CLS(3,:),'LineWidth',lns);
end

hold on; 
for i=1:length(g4)
    plot(XF(g4(i),:),'color',CLS(4,:),'LineWidth',lns);
end

subplot(2,4,3)
hold on
plot(mean(XF(g1,:),1),'k--','LineWidth',lns+1);
for i=1:length(g1), plot(XF(g1(i),:),'color',CLS(1,:),'LineWidth',lns); end
ylabel('$Z$','FontSize',14);xlim([0 T+2]);ylim([-7 7]);set(gca,'FontSize',10,'FontName','Times'); hline(0,'--k')
xticks([0 35 71]);xticklabels({'0%','50%','100%'}); yticks([-5 0 5]);yticklabels({'T','0','D'});box('off')         
title('Group 1','FontSize',14)


subplot(2,4,4)
hold on
plot(mean(XF(g2,:),1),'k--','LineWidth',lns+1);
for i=1:length(g2), plot(XF(g2(i),:),'color',CLS(2,:),'LineWidth',lns); end
ylabel('$Z$','FontSize',14);xlim([0 T+2]);ylim([-7 7]);set(gca,'FontSize',10,'FontName','Times'); hline(0,'--k')
xticks([0 35 71]);xticklabels({'0%','50%','100%'}); yticks([-5 0 5]);yticklabels({'T','0','D'});box('off')         
title('Group 2','FontSize',14)


subplot(2,4,7)
hold on
plot(mean(XF(g3,:),1),'k--','LineWidth',lns+1);
for i=1:length(g3), plot(XF(g3(i),:),'color',CLS(3,:),'LineWidth',lns); end
ylabel('$Z$','FontSize',14);xlim([0 T+2]);ylim([-7 7]);set(gca,'FontSize',10,'FontName','Times'); hline(0,'--k')
xticks([0 35 71]);xticklabels({'0%','50%','100%'}); yticks([-5 0 5]);yticklabels({'T','0','D'});box('off')         
title('Group 3','FontSize',14)


subplot(2,4,8)
hold on
plot(mean(XF(g4,:),1),'k--','LineWidth',lns+1);
for i=1:length(g4), plot(XF(g4(i),:),'color',CLS(4,:),'LineWidth',lns); end
ylabel('$Z$','FontSize',14);xlim([0 T+2]);ylim([-7 7]);set(gca,'FontSize',10,'FontName','Times'); hline(0,'--k')
xticks([0 35 71]);xticklabels({'0%','50%','100%'}); yticks([-5 0 5]);yticklabels({'T','0','D'});box('off')         
title('Group 4','FontSize',14)

%% INDIVIDUAL ANALYSIS: 1b
[b,gamma1,eta,delta] = switchModel('interaction',data.D1,max(data.x),post_means');
b = [b(1:24)+0.1; b(25:48)-0.25; b(49:end)];
[XF,PF] = GaussApprox_filter(data.Z,pars.sigmax,pars.a,b,0);
P = twoPL(XF,pars.a,b);
lns=1.25;

for i=1:length(g1), PBS(i,:) = [trapz(15:35,P(g1(i),1,15:35))/20 trapz(15:35,P(g1(i),25,15:35))/20]; end

set(gcf,'units','points','position',[100,100,1200,400],'color','w','defaultTextInterpreter','latex','defaultLegendInterpreter','latex');

subplot(1,2,1)
rectangle('Position',[15,0.01,20,1],'FaceColor',[0.9 .9 .9],'EdgeColor','white'); text(15,0.9,'$|-- \Delta_{30-50} --|$','FontSize',15); 
for i=1:length(g1), hold on; plot([0.5;squeeze(P(g1(i),1,2:T))],'LineWidth',lns); end
xlim([0 T+2]);ylim([0 1.05]);set(gca,'FontSize',10,'FontName','Times'); hline(0.5,'--k')
xticks([0 35 71]);xticklabels({'0%','50%','100%'}); title('HF','FontSize',15); %yticks([-10 0 10]);yticklabels({'T','0','D'});         
ylabel('$\pi_D$','FontSize',15);box('off')
lgd = legend({['p_\Delta(D) = ' num2str(round(PBS(1,1),2)) '   (Part. 6)'],['p_\Delta(D) = ' num2str(round(PBS(2,1),2)) '   (Part. 7)'],['p_\Delta(D) = ' num2str(round(PBS(3,1),2)) '   (Part. 8)'],['p_\Delta(D) = ' num2str(round(PBS(4,1),2)) '   (Part. 15)']},'Location','northeast','FontSize',12,'Interpreter','Tex');legend('boxoff');

subplot(1,2,2)
rectangle('Position',[15,0.01,20,1],'FaceColor',[0.9 .9 .9],'EdgeColor','white'); text(15,0.9,'$|-- \Delta_{30-50} --|$','FontSize',15); 
for i=1:length(g1), hold on; plot([0.5;squeeze(P(g1(i),25,2:T))],'LineWidth',lns); end
xlim([0 T+2]);ylim([0 1.05]);set(gca,'FontSize',10,'FontName','Times'); hline(0.5,'--k')
xticks([0 35 71]);xticklabels({'0%','50%','100%'}); title('LF','FontSize',15); %yticks([-10 0 10]);yticklabels({'T','0','D'});         
ylabel('$\pi_D$','FontSize',15);box('off')
lgd = legend({['p_\Delta(D) = ' num2str(round(PBS(1,2),2)) '   (Part. 6)'],['p_\Delta(D) = ' num2str(round(PBS(2,2),2)) '   (Part. 7)'],['p_\Delta(D) = ' num2str(round(PBS(3,2),2)) '   (Part. 8)'],['p_\Delta(D) = ' num2str(round(PBS(4,2),2)) '   (Part. 15)']},'Location','northeast','FontSize',12,'Interpreter','Tex');legend('boxoff');
 
%%

for i=1:length(g1), PBS(i,:) = [trapz(P(g1(i),1,:))/70 trapz(P(g1(i),25,:))/70];end; PBSS(1,:) = mean(PBS);
for i=1:length(g2), PBS(i,:) = [trapz(P(g2(i),1,:))/70 trapz(P(g2(i),25,:))/70];end; PBSS(2,:) = mean(PBS);
for i=1:length(g3), PBS(i,:) = [trapz(P(g3(i),1,:))/70 trapz(P(g3(i),25,:))/70];end; PBSS(3,:) = mean(PBS);
for i=1:length(g4), PBS(i,:) = [trapz(P(g4(i),1,:))/70 trapz(P(g4(i),25,:))/70];end; PBSS(4,:) = mean(PBS);

set(gcf,'units','points','position',[100,100,1200,400],'color','w','defaultTextInterpreter','latex','defaultLegendInterpreter','latex');

lns=2;

subplot(1,2,1)
hold on
plot([0.5;squeeze(mean(P(g1,5,2:T),1))],'LineWidth',lns);
plot([0.5;squeeze(mean(P(g2,5,2:T),1))],'LineWidth',lns);
plot([0.5;squeeze(mean(P(g3,5,2:T),1))],'LineWidth',lns);
plot([0.5;squeeze(mean(P(g4,5,2:T),1))],'LineWidth',lns);
xlim([0 T+2]);ylim([0 1.05]);set(gca,'FontSize',14,'FontName','Times'); hline(0.5,'--k')
xticks([0 35 71]);xticklabels({'0%','50%','100%'}); title('HF (high bigram)','FontSize',17); %yticks([-10 0 10]);yticklabels({'T','0','D'});         
ylabel('$\pi_D$','FontSize',18);box('off')
lgd = legend({['p(D) = ' num2str(round(PBSS(1,1),2)) '   (group 1)'],['p(D) = ' num2str(round(PBSS(2,1),2)) '   (group 2)'],['p(D) = ' num2str(round(PBSS(3,1),2)) '   (group 3)'],['p(D) = ' num2str(round(PBS(4,1),2)) '   (group 4)']},'Location','northeast','FontSize',16,'Interpreter','Tex');legend('boxoff');

subplot(1,2,2)
hold on
plot([0.5;squeeze(mean(P(g1,25,2:T),1))],'LineWidth',lns);
plot([0.5;squeeze(mean(P(g2,25,2:T),1))],'LineWidth',lns);
plot([0.5;squeeze(mean(P(g3,25,2:T),1))],'LineWidth',lns);
plot([0.5;squeeze(mean(P(g4,25,2:T),1))],'LineWidth',lns);
xlim([0 T+2]);ylim([0 1.05]);set(gca,'FontSize',14,'FontName','Times'); hline(0.5,'--k')
xticks([0 35 71]);xticklabels({'0%','50%','100%'}); title('LF (high bigram)','FontSize',17); %yticks([-10 0 10]);yticklabels({'T','0','D'});         
ylabel('$\pi_D$','FontSize',18);box('off')
lgd = legend({['p(D) = ' num2str(round(PBSS(1,2),2)) '   (group 1)'],['p(D) = ' num2str(round(PBSS(2,2),2)) '   (group 2)'],['p(D) = ' num2str(round(PBSS(3,2),2)) '   (group 3)'],['p(D) = ' num2str(round(PBS(4,2),2)) '   (group 4)']},'Location','northeast','FontSize',16,'Interpreter','Tex');legend('boxoff');

%% INDIVIDUAL ANALYSIS: 2 (by fixing bigram to low level, -3 level)
[b,gamma1,eta,delta] = switchModel('regression',data.D1,(data.x),post_means');
%b = data.D1*post_means(1:3)' + data.x*post_means(4);
[XF,PF] = GaussApprox_filter(data.Z,pars.sigmax,pars.a,b,0);
P = twoPL(XF,pars.a,b);

for i=1:length(g1), PBS(i,:) = [trapz(15:35,P(g1(i),1,15:35))/20 trapz(15:35,P(g1(i),25,15:35))/20]; end

%stms_bb_selected = [20 8 12 19]; %ritmo epoca latte ponte
j=5; %colpa
PP = P(g1,[j],:); %select probs for individuals in g1 and stimuli in the cell "lowbigram & LF"

set(gcf,'units','points','position',[100,100,1350,550],'color','w','defaultTextInterpreter','latex','defaultLegendInterpreter','latex');
lns=1.25;

for i=1:size(PP,1)     
    subplot(2,4,i)
    hold on
    %plot([0;squeeze(log(PP(i,1,2:T)./(1-PP(i,1,2:T))))]);    
    %rectangle('Position',[15,0.01,20,1],'FaceColor',[0.95 .95 .95],'EdgeColor','white'); text(13.5,0.9,'$|-\Delta_{30-50}-|$','FontSize',11)
    plot([0.5;squeeze((PP(i,1,3:T)))],'LineWidth',lns);
    xlim([0 T+2]);ylim([0 1.1]);set(gca,'FontSize',10,'FontName','Times'); hline(0.5,'--k')
    xticks([0 35 71]);xticklabels({'0%','50%','100%'}); title(['Participant ' num2str(g1(i))],'FontSize',15); %yticks([-10 0 10]);yticklabels({'T','0','D'});         
    ylabel('$\pi_D$','FontSize',15);box('off')
    %ylabel('$\log(\pi_D/\pi_T)$','FontSize',14)
    
    subplot(2,4,i+4);
    polarhistogram(squeeze(Y(g1(i),j,1:T)),3);title('D~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~T','FontSize',14)
    %hist(squeeze(Y(g1(i),j,:)))
end
subplot(2,4,8);polarhistogram(squeeze(Y(g1(i),j,1:T)),2);title('D~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~T','FontSize',14)

j=70;
PP = P(g1,[j],:); %select probs for individuals in g1 and stimuli in the cell "lowbigram & LF"

% figure();
% set(gcf,'units','points','position',[100,100,1350,550],'color','w','defaultTextInterpreter','latex','defaultLegendInterpreter','latex');
% lns=1.25;

for i=1:size(PP,1)     
    subplot(2,4,i)
    hold on
    %plot([0;squeeze(log(PP(i,1,2:T)./(1-PP(i,1,2:T))))]);    
    %rectangle('Position',[15,0.01,20,1],'FaceColor',[0.95 .95 .95],'EdgeColor','white'); text(13.5,0.9,'$|-\Delta_{30-50}-|$','FontSize',11)
    plot([0.5;squeeze((PP(i,1,3:T)))],'LineWidth',lns);
    xlim([0 T+2]);ylim([0 1.1]);set(gca,'FontSize',10,'FontName','Times'); hline(0.5,'--k')
    xticks([0 35 71]);xticklabels({'0%','50%','100%'}); title(['Participant ' num2str(g1(i))],'FontSize',15); %yticks([-10 0 10]);yticklabels({'T','0','D'});         
    ylabel('$\pi_D$','FontSize',15);box('off')
    %ylabel('$\log(\pi_D/\pi_T)$','FontSize',14)
    
    subplot(2,4,i+4);hold on
    polarhistogram(squeeze(Y(g1(i),j,1:T)),2);title('D~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~T','FontSize',14)
    %hist(squeeze(Y(g1(i),j,:)))
end

%%

set(gcf,'units','points','position',[100,100,1350,550],'color','w','defaultTextInterpreter','latex','defaultLegendInterpreter','latex');
lns=1.25;
CLS = colormap(lines(4));

j=[5 70]; %colpa
PP = P(g1,[j],:); %select probs for individuals in g1 and stimuli in the cell "lowbigram & LF"

kk=ones(1,70)*0.01; kk(3:40) = 0.12;

subplot(2,8,[1:2])
hold on
plot([0.5;squeeze((PP(1,1,3:T)))],'LineWidth',lns); plot([0.5;squeeze((PP(1,2,3:T)))]+kk','LineWidth',lns);
xlim([0 T+2]);ylim([-0.05 1.1]);set(gca,'FontSize',14,'FontName','Times'); hline(0.5,'--k')
xticks([0 35 71]);xticklabels({'0%','50%','100%'}); title(['Participant ' num2str(g1(1))],'FontSize',18); %yticks([-10 0 10]);yticklabels({'T','0','D'});         
ylabel('$\pi_D$','FontSize',15);box('off')
legend({'HF word','LF word'},'Location','northeast','FontSize',14,'Interpreter','Tex');legend('boxoff');

subplot(2,8,[3:4])
hold on
plot([0.5;squeeze((PP(2,1,3:T)))],'LineWidth',lns); plot([0.5;squeeze((PP(2,2,3:T)))]+kk','LineWidth',lns);
xlim([0 T+2]);ylim([-0.05 1.1]);set(gca,'FontSize',14,'FontName','Times'); hline(0.5,'--k')
xticks([0 35 71]);xticklabels({'0%','50%','100%'}); title(['Participant ' num2str(g1(2))],'FontSize',18); %yticks([-10 0 10]);yticklabels({'T','0','D'});         
%ylabel('$\pi_D$','FontSize',15);box('off')
legend({'HF word','LF word'},'Location','northeast','FontSize',14,'Interpreter','Tex');legend('boxoff');

subplot(2,8,[5:6])
hold on
plot([0.5;squeeze((PP(3,1,3:T)))],'LineWidth',lns); plot([0.5;squeeze((PP(3,2,3:T)))]+kk','LineWidth',lns);
xlim([0 T+2]);ylim([-0.05 1.1]);set(gca,'FontSize',14,'FontName','Times'); hline(0.5,'--k')
xticks([0 35 71]);xticklabels({'0%','50%','100%'}); title(['Participant ' num2str(g1(3))],'FontSize',18); %yticks([-10 0 10]);yticklabels({'T','0','D'});         
%ylabel('$\pi_D$','FontSize',15);box('off')
legend({'HF word','LF word'},'Location','northeast','FontSize',14,'Interpreter','Tex');legend('boxoff');

subplot(2,8,[7:8])
hold on
plot([0.5;squeeze((PP(4,1,3:T)))],'LineWidth',lns); plot([0.5;squeeze((PP(4,2,3:T)))]+kk','LineWidth',lns);
xlim([0 T+2]);ylim([0 1.1]);set(gca,'FontSize',14,'FontName','Times'); hline(0.5,'--k')
xticks([0 35 71]);xticklabels({'0%','50%','100%'}); title(['Participant ' num2str(g1(4))],'FontSize',18); %yticks([-10 0 10]);yticklabels({'T','0','D'});         
%ylabel('$\pi_D$','FontSize',15);box('off')
legend({'HF word','LF word'},'Location','northeast','FontSize',14,'Interpreter','Tex');legend('boxoff');

subplot(2,8,9);polarhistogram(squeeze(Y(g1(1),j(1),1:T)),2);title('D~~~~~~~~~~~~~T','FontSize',14)
subplot(2,8,10);polarhistogram(squeeze(Y(g1(1),j(2),1:T)),2,'FaceColor',CLS(2,:));title('D~~~~~~~~~~~~~T','FontSize',14)

subplot(2,8,11);polarhistogram(squeeze(Y(g1(2),j(1),1:T)),2);title('D~~~~~~~~~~~~~T','FontSize',14)
subplot(2,8,12);polarhistogram(squeeze(Y(g1(2),j(2),1:T)),2,'FaceColor',CLS(2,:));title('D~~~~~~~~~~~~~T','FontSize',14)

subplot(2,8,13);polarhistogram(squeeze(Y(g1(3),j(1),1:T)),2);title('D~~~~~~~~~~~~~T','FontSize',14)
subplot(2,8,14);polarhistogram(squeeze(Y(g1(3),j(2),1:T)),2,'FaceColor',CLS(2,:));title('D~~~~~~~~~~~~~T','FontSize',14)

subplot(2,8,15);polarhistogram(squeeze(Y(g1(4),j(1),1:T)),2);title('D~~~~~~~~~~~~~T','FontSize',14)
subplot(2,8,16);polarhistogram(squeeze(Y(g1(4),j(2),1:T)),2,'FaceColor',CLS(2,:));title('D~~~~~~~~~~~~~T','FontSize',14)

%%

% Evaluate posteriors
% figure()
% subplot(2,3,1);hold on;ksdensity(S_coll(:,1));ksdensity(S_coll(:,2)); title('HF vs. LF');
% subplot(2,3,2);hold on;ksdensity(S_coll(:,1));ksdensity(S_coll(:,3)); title('HF vs. NW'); 
% subplot(2,3,3);hold on;ksdensity(S_coll(:,2));ksdensity(S_coll(:,3)); title('LF vs. NW'); 
% subplot(2,3,4);hold on;ksdensity(S_coll(:,4));ksdensity(S_coll(:,5)); title('HF:x vs. LF:x');
% subplot(2,3,5);hold on;ksdensity(S_coll(:,4));ksdensity(S_coll(:,6)); title('HF:x vs. NW:x');
% subplot(2,3,6);hold on;ksdensity(S_coll(:,5));ksdensity(S_coll(:,6)); title('LF:x vs. NW:x');
% 
% % evaluate beta in terms of decomposed coefficients
% x_t = linspace(-bnd,bnd,100); 
% figure(); cls={'r','b','g'};
% subplot(1,3,1);hold on; plot(x_t,squeeze(twoPL(x_t,1, gamma1(1)+min(data.x)*(eta+delta(1)) ))',cls{1}); plot(x_t,squeeze(twoPL(x_t,1, gamma1(2)+min(data.x)*(eta+delta(2)) ))',cls{2}); plot(x_t,squeeze(twoPL(x_t,1, gamma1(3)+min(data.x)*(eta+delta(3)) ))',cls{3}); title('Low-bigramF')
% subplot(1,3,2);hold on; plot(x_t,squeeze(twoPL(x_t,1, gamma1(1)+0.01*(eta+delta(1)) ))',cls{1}); plot(x_t,squeeze(twoPL(x_t,1, gamma1(2)+0.01*(eta+delta(2)) ))',cls{2}); plot(x_t,squeeze(twoPL(x_t,1, gamma1(3)+0.01*(eta+delta(3)) ))',cls{3}); title('0-bigramF')
% subplot(1,3,3);hold on; plot(x_t,squeeze(twoPL(x_t,1, gamma1(1)+max(data.x)*(eta+delta(1)) ))',cls{1}); plot(x_t,squeeze(twoPL(x_t,1, gamma1(2)+max(data.x)*(eta+delta(2)) ))',cls{2}); plot(x_t,squeeze(twoPL(x_t,1, gamma1(3)+max(data.x)*(eta+delta(3)) ))',cls{3});title('High-bigramF')
% legend('HF','LF','NW')
% 
% % evaluate beta in terms of decomposed coefficients
% x_t = linspace(-bnd,bnd,100); 
% figure(); cls={'r','b','g'};
% subplot(1,3,1);hold on; plot(x_t,squeeze(twoPL(x_t,1, gamma1(1)+min(data.x)*(eta+delta(1)) ))',cls{1}); plot(x_t,squeeze(twoPL(x_t,1, gamma1(1)+0.01*(eta+delta(1)) ))',cls{2}); plot(x_t,squeeze(twoPL(x_t,1, gamma1(1)+max(data.x)*(eta+delta(1)) ))',cls{3}); title('HF')
% subplot(1,3,2);hold on; plot(x_t,squeeze(twoPL(x_t,1, gamma1(2)+min(data.x)*(eta+delta(2)) ))',cls{1}); plot(x_t,squeeze(twoPL(x_t,1, gamma1(2)+0.01*(eta+delta(2)) ))',cls{2}); plot(x_t,squeeze(twoPL(x_t,1, gamma1(2)+max(data.x)*(eta+delta(2)) ))',cls{3}); title('LF')
% subplot(1,3,3);hold on; plot(x_t,squeeze(twoPL(x_t,1, gamma1(3)+min(data.x)*(eta+delta(3)) ))',cls{1}); plot(x_t,squeeze(twoPL(x_t,1, gamma1(3)+0.01*(eta+delta(3)) ))',cls{2}); plot(x_t,squeeze(twoPL(x_t,1, gamma1(3)+max(data.x)*(eta+delta(3)) ))',cls{3}); title('NW')
% legend('Low-bigramF','0-bigramF','High-bigramF')
% 
% 





