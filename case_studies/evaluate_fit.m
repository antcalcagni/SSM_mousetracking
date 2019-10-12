%% Set environment
addpath(genpath('SSM_MouseTracking/'))
clear all; close all

load('barca2012_nomissin_full.mat');
ncores=2; M=4; 

parpool;

%% Model 5_0 (categorical)
tic
disp(' ');disp('@@ Evaluate model: 5_0')

load('mh_5_0.mat')

q0=1;q=72; %select just {HF,LF,NW}
data.Z=Z(:,q0:q,:); data.Y=Y(:,q0:q,:); data.D1=D1(q0:q,1:3); pars.sigmax=sigmax; pars.kappa1=kappa1; pars.kappa2=kappa2; pars.a=a(q0:q); 
pars.mu1=mu1; pars.mu2=mu2; pars.bnd=bnd; data.x=scaledata(bigF_cov(q0:q),-3,3);

[S_coll, post_means,Rgelm,~,S] = MCMC_diag(mh_5_0,1000,1,false,true,'regression',1);
b = data.x*post_means;
[XF,PF] = GaussApprox_filter(data.Z,pars.sigmax,pars.a,b,0);
P = twoPL(XF,pars.a,b);

% Generated data
if mod(M/ncores,1)==0, iid = repmat(M/ncores,1,ncores);else iid = repmat(floor(M/ncores),1,ncores-1); iid = [iid M-sum(iid)]; end
disp(['@@ Datasets by core: ',num2str(iid)])

parfor n=1:ncores
    datagen{n} = arrayfun(@(i)generateData(T,I,size(data.Y,2),mu1,mu2,kappa1,kappa2,b,ones(size(data.Y,2),1),thr,sigmax,bnd),1:iid(n),'UniformOutput',false);
end
datagen_agg=cell(0); for n=1:ncores, datagen_agg = [datagen_agg datagen{n}]; end

% Evaluate data
[PA_ov,AOR_time,AOR_sbj,PA_sbj] = compute_fit(M,datagen_agg,data);

% Save results
save('mh_5_0_fit.mat','PA_ov','AOR_time','PA_sbj','datagen_agg','-v7.3')
disp(['@@ Evaluate model: done. Elapsed time: ' num2str(toc)])

%% Model 5_1 (regression)
tic
disp(' ');disp('@@ Evaluate model: 5_1')

load('mh_5_1.mat')

q0=1;q=72; %select just {HF,LF,NW}
data.Z=Z(:,q0:q,:); data.Y=Y(:,q0:q,:); data.D1=D1(q0:q,1:3); pars.sigmax=sigmax; pars.kappa1=kappa1; pars.kappa2=kappa2; pars.a=a(q0:q); 
pars.mu1=mu1; pars.mu2=mu2; pars.bnd=bnd; data.x=scaledata(bigF_cov(q0:q),-3,3);

[S_coll, post_means,Rgelm,neff,S] = MCMC_diag(mh_5_1,1000,1,false,true,'categorical',3);
b = data.D1*post_means';
[XF,PF] = GaussApprox_filter(data.Z,pars.sigmax,pars.a,b,0);
P = twoPL(XF,pars.a,b);

% Generated data
if mod(M/ncores,1)==0, iid = repmat(M/ncores,1,ncores);else iid = repmat(floor(M/ncores),1,ncores-1); iid = [iid M-sum(iid)]; end
disp(['@@ Datasets by core: ',num2str(iid)])

parfor n=1:ncores
    datagen{n} = arrayfun(@(i)generateData(T,I,size(data.Y,2),mu1,mu2,kappa1,kappa2,b,ones(size(data.Y,2),1),thr,sigmax,bnd),1:iid(n),'UniformOutput',false);
end
datagen_agg=cell(0); for n=1:ncores, datagen_agg = [datagen_agg datagen{n}]; end

% Evaluate data
[PA_ov,AOR_time,PA_sbj] = compute_fit(M,datagen_agg,data);

% Save results
save('mh_5_1_fit.mat','PA_ov','AOR_time','PA_sbj','datagen_agg','-v7.3')
disp(['@@ Evaluate model: done. Elapsed time: ' num2str(toc)])

%% Model 5_2 (interaction)
tic
disp(' ');disp('@@ Evaluate model: 5_2')

load('mh_5_2.mat')

q0=1;q=72; %select just {HF,LF,NW}
data.Z=Z(:,q0:q,:); data.Y=Y(:,q0:q,:); data.D1=D1(q0:q,1:3); pars.sigmax=sigmax; pars.kappa1=kappa1; pars.kappa2=kappa2; pars.a=a(q0:q); 
pars.mu1=mu1; pars.mu2=mu2; pars.bnd=bnd; data.x=scaledata(bigF_cov(q0:q),-3,3);

[S_coll, post_means,Rgelm,~,S] = MCMC_diag(mh_5_2,2250,3,false,true,'interaction',3);
[b,gamma1,eta,delta] = switchModel('interaction',data.D1,data.x,post_means');

[XF,PF] = GaussApprox_filter(data.Z,pars.sigmax,pars.a,b,0);
P = twoPL(XF,pars.a,b);

% Generated data
if mod(M/ncores,1)==0, iid = repmat(M/ncores,1,ncores);else iid = repmat(floor(M/ncores),1,ncores-1); iid = [iid M-sum(iid)]; end
disp(['@@ Datasets by core: ',num2str(iid)])

parfor n=1:ncores
    datagen{n} = arrayfun(@(i)generateData(T,I,size(data.Y,2),mu1,mu2,kappa1,kappa2,b,ones(size(data.Y,2),1),thr,sigmax,bnd),1:iid(n),'UniformOutput',false);
end
datagen_agg=cell(0); for n=1:ncores, datagen_agg = [datagen_agg datagen{n}]; end

% Evaluate data
[PA_ov,AOR_time,PA_sbj] = compute_fit(M,datagen_agg,data);

% Save results
save('mh_5_2_fit.mat','PA_ov','AOR_time','PA_sbj','datagen_agg','-v7.3')
disp(['@@ Evaluate model: done. Elapsed time: ' num2str(toc)])

%%
delete(gcp('nocreate')); %for matlab2017


%% Evaluate fit
clear all
load('barca2012_nomissin_full.mat');

load('mh_5_0.mat')
load('/home/antonio/Scrivania/barca2012_fit/mh_5_0_fit.mat')

q0=1;q=72; %select just {HF,LF,NW}
data.Z=Z(:,q0:q,:); data.Y=Y(:,q0:q,:); data.D1=D1(q0:q,1:3); pars.sigmax=sigmax; pars.kappa1=kappa1; pars.kappa2=kappa2; pars.a=a(q0:q); 
pars.mu1=mu1; pars.mu2=mu2; pars.bnd=bnd; data.x=scaledata(bigF_cov(q0:q),-3,3);
[S_coll, post_means,Rgelm,~,S] = MCMC_diag(mh_5_0,1000,1,false,true,'regression',1);
b = data.x*post_means;
[XF,PF] = GaussApprox_filter(data.Z,pars.sigmax,pars.a,b,0);
P = twoPL(XF,pars.a,b);
[PA_ov1,AOR_time1,AOR_sbj1,PA_sbj1] = compute_fit(350,datagen_agg,data);

load('mh_5_1.mat')
load('/home/antonio/Scrivania/barca2012_fit/mh_5_1_fit.mat')

q0=1;q=72; %select just {HF,LF,NW}
data.Z=Z(:,q0:q,:); data.Y=Y(:,q0:q,:); data.D1=D1(q0:q,1:3); pars.sigmax=sigmax; pars.kappa1=kappa1; pars.kappa2=kappa2; pars.a=a(q0:q); 
pars.mu1=mu1; pars.mu2=mu2; pars.bnd=bnd; data.x=scaledata(bigF_cov(q0:q),-3,3);
[S_coll, post_means,Rgelm,neff,S] = MCMC_diag(mh_5_1,1000,1,false,true,'categorical',3);
b = data.D1*post_means';
[XF,PF] = GaussApprox_filter(data.Z,pars.sigmax,pars.a,b,0);
P = twoPL(XF,pars.a,b);
[PA_ov2,AOR_time2,AOR_sbj2,PA_sbj2] = compute_fit(350,datagen_agg,data);

load('mh_5_2.mat')
load('/home/antonio/Scrivania/barca2012_fit/mh_5_2_fit.mat')

q0=1;q=72; %select just {HF,LF,NW}
data.Z=Z(:,q0:q,:); data.Y=Y(:,q0:q,:); data.D1=D1(q0:q,1:3); pars.sigmax=sigmax; pars.kappa1=kappa1; pars.kappa2=kappa2; pars.a=a(q0:q); 
pars.mu1=mu1; pars.mu2=mu2; pars.bnd=bnd; data.x=scaledata(bigF_cov(q0:q),-3,3);
[S_coll, post_means,Rgelm,~,S] = MCMC_diag(mh_5_2,2250,3,false,true,'interaction',3);
[b,gamma1,eta,delta] = switchModel('interaction',data.D1,data.x,post_means');
[XF,PF] = GaussApprox_filter(data.Z,pars.sigmax,pars.a,b,0);
P = twoPL(XF,pars.a,b);
[PA_ov3,AOR_time3,AOR_sbj3,PA_sbj3] = compute_fit(350,datagen_agg,data);

%%
CLS = get(gca,'colororder');close all
fnts=20; fntst=22; fntl=18; 
figure(); nfig=1;
set(gcf,'units','points','position',[100,100,1200,325],'color','w','defaultTextInterpreter','latex','defaultLegendInterpreter','latex');

subplot1 = subplot(1,3,[1 1]); histogram(PA_ov1); vline(mean(PA_ov1),'r--'); box('off'); title('(A) overall AoR','FontSize',18)
%set(subplot1,'FontSize',fntl,'FontName','Times'); xlabel('overall AoR','FontSize',fnts,'FaceColor');

subplot2 = subplot(1,3,[2 3]); boxplot(AOR_sbj1,'PlotStyle','compact','Color',CLS(1,:)); hline(mean(mean(PA_ov1)),'k--'); box('off'); title('(B) by-subject AoR','FontSize',18)
%set(subplot2,'FontSize',fntl,'FontName','Times'); xlabel('by-subject AoR','FontSize',fnts,'FaceColor');

saveas(gcf,['fit_model' num2str(nfig)],'epsc');

figure(); nfig=2;
set(gcf,'units','points','position',[100,100,1200,325],'color','w','defaultTextInterpreter','latex','defaultLegendInterpreter','latex');

subplot1 = subplot(1,3,[1 1]); histogram(PA_ov2); vline(mean(PA_ov2),'r--'); box('off'); title('(A) overall AoR','FontSize',18)
%set(subplot1,'FontSize',fntl,'FontName','Times'); xlabel('overall AoR','FontSize',fnts,'FaceColor');

subplot2 = subplot(1,3,[2 3]); boxplot(AOR_sbj2,'PlotStyle','compact','Color',CLS(1,:)); hline(mean(mean(PA_ov2)),'k--'); box('off'); title('(B) by-subject AoR','FontSize',18)
%set(subplot2,'FontSize',fntl,'FontName','Times'); xlabel('by-subject AoR','FontSize',fnts,'FaceColor');


saveas(gcf,['fit_model' num2str(nfig)],'epsc');

figure(); nfig=3;
set(gcf,'units','points','position',[100,100,1200,325],'color','w','defaultTextInterpreter','latex','defaultLegendInterpreter','latex');

subplot1 = subplot(1,3,[1 1]); histogram(PA_ov3); vline(mean(PA_ov3),'r--'); box('off'); title('(A) overall AoR','FontSize',18)
%set(subplot1,'FontSize',fntl,'FontName','Times'); xlabel('overall AoR','FontSize',fnts,'FaceColor');

subplot2 = subplot(1,3,[2 3]); boxplot(AOR_sbj3,'PlotStyle','compact','Color',CLS(1,:)); hline(mean(mean(PA_ov3)),'k--'); box('off'); title('(B) by-subject AoR','FontSize',18)
%set(subplot2,'FontSize',fntl,'FontName','Times'); xlabel('by-subject AoR','FontSize',fnts,'FaceColor');

saveas(gcf,['fit_model' num2str(nfig)],'epsc');













