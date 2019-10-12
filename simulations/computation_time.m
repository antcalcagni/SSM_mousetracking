%%
clear all
addpath(genpath('SSM_MouseTracking/'))
ncores_tot=2;  %total number of logical cores
ncores_ch=1;
parpool;

%% Define simulation scenario
S = allcomb({'categorical','interaction'},{12 25 50},{61},{12 28},{2 4}); %typeModel x I x T x J x K

%% Set simulation parameters
nsamples=1; nchains=1000;
mu1=2.75;mu2=0.75;kappa1=200;kappa2=200; bnd=5;

prior_mu_b=cell(0); prior_var_b=cell(0);
for s=1:size(S,1)
    K=cell2mat(S(s,5));
    if cell2mat(S(s,1))=="categorical"
        prior_mu_b{s} = {unifrnd(-2,2,K,1)};
        prior_var_b{s} = {ones(1,K,1)*0.5};
    else
        prior_mu_b{s} = {[unifrnd(-2,2,K,1);unifrnd(-1,1,1,1);unifrnd(-2,2,K-1,1)]};
        prior_var_b{s} = {[ones(1,K,1)*0.5 0.25 ones(1,K-1,1)*0.45]'};
    end        
end

%% Run Simulations
for s=1:size(S,1)
    tic
    disp(' ');disp(['@ Simulation: ' num2str(s) '/' num2str(size(S,1))])
    
    I=cell2mat(S(s,2)); T=cell2mat(S(s,3)); J=cell2mat(S(s,4)); K=cell2mat(S(s,5));     
    D1 = kron(eye(K),ones(J/K,1)); %generate a well-distribuited classification matrix
    x = sort(unifrnd(-3,3,J,1));
    sigmax = repmat(5.5,[1 I]);

    typeModel=cell2mat(S(s,1));
        
    disp('@@ Generate data... ')
    [datagen,B,eta,gamma,delta] = generateModels(ncores_tot,nsamples,T,I,J,K,mu1,mu2,kappa1,kappa2,typeModel,prior_mu_b{s}{1},prior_var_b{s}{1},x,D1,sigmax,bnd);    
    disp('@@ Generate data... done.')
    
    disp(['@@ Estimate the model..'])   
    for n=1:nsamples        
        datacurr = datagen{n}; data.Z=datacurr.Z; data.Y=datacurr.Y; data.D1=D1; data.x=x; pars.sigmax=sigmax'; pars.kappa1=kappa1; pars.kappa2=kappa2; pars.a=ones(J,1); pars.mu1=mu1; pars.mu2=mu2; pars.bnd=bnd; 
        W = length(prior_mu_b{s}{1}); %number of parameters to sample        
        mh{n} = MH_adapt_multicore(W,data,pars,nchains,1000,25,[],false,'type1',ncores_ch,typeModel,true);
    end
    elapst(s)=toc;
    disp(['@ Elapsed time: ' num2str(elapst(s))])    
end
save('comput_time.mat','mh','s','-v7.3')    

%%
delete(gcp('nocreate')); %for matlab2017

%% Analysis of results
load('comptime.mat')

xlbls = cell(0);
for i=1:size(S,1), xlbls(i) = {['I=',num2str(cell2mat(S(i,2))), ',J=',num2str(cell2mat(S(i,4))), ',K=',num2str(cell2mat(S(i,5)))]}; end

fnts=20; fntst=22; fntl=15; 
figure(); nfig=1;
set(gcf,'units','points','position',[100,100,1500,425],'color','w','defaultTextInterpreter','latex','defaultLegendInterpreter','latex');

subplot1 = subplot(1,2,1); plot(CTIME(1:12,1),'o--','LineWidth',2,'Parent',subplot1); hold on; 
plot(CTIME(1:12,2),'o--','LineWidth',2,'Parent',subplot1); box('off'); ylim([0 430]); xlim([0 13]);
hline([60 120 180 240 300 360 420],'k--'); legend({'Warming-up','Sampling'},'Location','NorthWest','FontSize',fntl);legend('boxoff'); 
set(subplot1,'FontSize',fntl,'FontName','Times','XTick',1:12); ylabel('Elapsed Time (sec)','FontSize',fnts); xlabel('scenario','FontSize',fnts);
text(0,450,'(A) Categorical model','FontSize',22,'FontWeight','bold'); 

subplot2 = subplot(1,2,2); plot(CTIME(13:end,1),'o--','LineWidth',2,'Parent',subplot2); hold on; 
plot(CTIME(1:12,2),'o--','LineWidth',2,'Parent',subplot2); box('off'); ylim([0 430]); xlim([0 13]);
hline([60 120 180 240 300 360 420],'k--'); legend({'Warming-up','Sampling'},'Location','NorthWest','FontSize',fntl);legend('boxoff'); 
set(subplot2,'FontSize',fntl,'FontName','Times','XTick',1:12); xlabel('scenario','FontSize',fnts); %ylabel('Elapsed Time (sec)','FontSize',fnts); 
text(0,450,'(B) Interaction model','FontSize',22,'FontWeight','bold'); 

saveas(gcf,['ctime.eps' num2str(nfig)],'epsc');















