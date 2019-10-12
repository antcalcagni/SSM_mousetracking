%% Set environment
addpath(genpath('SSM_MouseTracking/'))
clear all
close all

%% Generate datasets
nsamples=1000; ncores=2;

mu1=2.75;mu2=0.75;kappa1=200;kappa2=200; 
I=12; T=61; J=10; K=2; 
sigmax = repmat(5.5,[1 I]); bnd=5;
D1 = kron(eye(K),ones(J/K,1)); %generate a well-distribuited classification matrix
x = sort(unifrnd(-3,3,10,1));

typeModel="interaction";
prior_mu_b = [1.2;-0.9;1.5;-0.25]; prior_var_b = [0.5;0.5;0.25;0.5];

tic
[datagen,B,eta,gamma,delta] = generateModels(ncores,nsamples,T,I,J,K,mu1,mu2,kappa1,kappa2,typeModel,prior_mu_b,prior_var_b,x,D1,sigmax,bnd);
toc


%% Estimate the model
datacurr = datagen{1};
data.Z=datacurr.Z; data.Y=datacurr.Y; data.D1=D1; data.x=x;
pars.sigmax=sigmax'; pars.kappa1=kappa1; pars.kappa2=kappa2; pars.a=ones(J,1); 
pars.mu1=mu1; pars.mu2=mu2; pars.bnd=bnd; 

numcores=2;
[mh] = MH_adapt_multicore(length(prior_mu_b),data,pars,3000,200,25,[],true,'type1',numcores,typeModel,false,prior_mu_b+normrnd(0,0.1,length(prior_mu_b),1))

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Analyse results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
clear all; %parpool;

PA_tot = zeros(16,1);
Lambda = zeros(16,8);
Mu_post = zeros(16,8);
Var_post = zeros(16,8);
Hdpi_post = zeros(16,8,2);
ncores=2;M=250;

%s=1; %set scenario
for s=10:16
    disp(['@ Scenario: ' num2str(s)])

    %%% Load data
    load('/home/antonio/Scrivania/mhsamples/design.mat')
    design=S;

    load(['/home/antonio/Scrivania/mhsamples/mhsamples_' num2str(s) '.mat'])
    load(['/home/antonio/Scrivania/mhsamples/datagen_' num2str(s) '.mat'])

    %%% Re-arrange posteriors
    mhh=mh;
    clear S_coll_means S_coll_tot datagen_mod datagen_single Y1_vec Y0_vec
    for p=1:250
        mh=mhh{p};
        try 
            nsample=size(mh.samples{1}.samples,2);
            JJ=size(mh.samples{1}.samples,1);
            Q=length(mh.samples);

            S = zeros(nsample-1+1,JJ,Q);
            for q=1:Q, S(:,:,q) = mh.samples{q}.samples(:,1:end)'; end %arrange samples into a [Nsample x J x Q] matrix (for 'mcmcdiag' package)

            burn=1500; thin=1; %burning and thinning
            iid_thin=burn:thin:size(S,1);
            S = S(iid_thin,:,:);
            S_coll = []; for q=1:Q, S_coll=[S_coll;S(:,:,q)]; end %mixing-up MCMC chains

            S_coll_tot{p} = S_coll;
        catch
            S_coll_tot{p} = S_coll_tot{p-1};
        end        
    end
    for p=1:250, S_coll_means(:,p)=mean(S_coll_tot{p});end %compute means
    Mu_post(s,1:size(S_coll_means,1)) = mean(S_coll_means,2);
    Var_post(s,1:size(S_coll_means,1)) = var(S_coll_means,[],2);    
    Hdpi_post(s,1:size(S_coll_means,1),:) = hpdi(S_coll_means',95)';

    %%% Compute outcome measures: lambda
    xsup=linspace(-10,10,1000);
    for k=1:cell2mat(design(s,5))    
        [fy1] =ksdensity(S_coll_means(k,:),xsup);
        [fy0] = ksdensity(normrnd(prior_mu_b{s}{1}(k),prior_var_b{s}{1}(k),10000,1),xsup);
        %lambda(k) = trapz(min(fy1,fy0))/trapz(max(fy1,fy0));
        lambda(k) = trapz(min(fy1,fy0))/max(trapz(fy1),trapz(fy0));
    end
    Lambda(s,1:cell2mat(design(s,5))) = lambda(k); 

    %%% Compute outcome measures: AoR
    D1 = kron(eye(cell2mat(design(s,5))),ones(cell2mat(design(s,4))/cell2mat(design(s,5)),1)); %generate a well-distribuited classification matrix
    x = sort(unifrnd(-3,3,cell2mat(design(s,4)),1));
    disp(['@@ datagen..']);

    if mod(M/ncores,1)==0, iid = repmat(M/ncores,1,ncores);else iid = repmat(floor(M/ncores),1,ncores-1); iid = [iid M-sum(iid)]; end
    disp(['@@ Datasets by core: ',num2str(iid)])

    if strcmp(cell2mat(design(s,1)),'interaction')==1
        parfor w=1:ncores
            datagen_single{w} = arrayfun(@(p)generateData(61,cell2mat(design(s,2)),cell2mat(design(s,4)),2.75,0.75,200,200,D1*S_coll_means(1:cell2mat(design(s,5)),p) + x.*(S_coll_means(cell2mat(design(s,5))+1,p) + D1*[0;S_coll_means((2+cell2mat(design(s,5))):end,p)]),ones(cell2mat(design(s,4)),1),(2.75+0.75)/2,repmat(5.5,[1 cell2mat(design(s,2))]),5),1:250,'UniformOutput',false);
        end
    else
        parfor w=1:ncores
            datagen_single{w} = arrayfun(@(p)generateData(61,cell2mat(design(s,2)),cell2mat(design(s,4)),2.75,0.75,200,200,D1*S_coll_means(:,p),ones(cell2mat(design(s,4)),1),(2.75+0.75)/2,repmat(5.5,[1 cell2mat(design(s,2))]),5),1:250,'UniformOutput',false);
        end
    end
    datagen_mod=cell(0); for w=1:ncores, datagen_mod = [datagen_mod datagen_single{w}]; end

    for p=1:250    
        for t=1:cell2mat(design(s,2)), Y1vec(t,:) = vec(datagen_mod{p}.Y(:,:,t)')'; end
        for t=1:cell2mat(design(s,2)), Y0vec(t,:) = vec(datagen{p}.Y(:,:,t)')'; end

        PA_ov(p) = 1 - (norm(Y0vec - Y1vec)^2 / (norm(Y0vec)^2));
    end

    PA_tot(s) = mean(PA_ov);
    disp('@@')
    save(['sim_results_' num2str(s) '.mat'],'datagen_mod','PA_ov','lambda','S_coll_tot','-v7.3')
end
%save('sim_results.mat','Lambda','PA_tot','Mu_post','Var_post','-v7.3')


%% Analyse results (using partial results, for each s=1..16)
clear all
addpath(genpath('SSM_MouseTracking/'))
addpath('/home/antonio/MEGA/Lavoro_sync/MATLAB/mcmcdiag')

PA_tot = zeros(16,1); Lambda = zeros(16,8); Mu_post = zeros(16,8);
Var_post = zeros(16,8); Hdpi_post = zeros(16,8,2); Mu_true = zeros(16,8);

for s=1:16
    clear S_coll_means
    disp(['Design s: ' num2str(s)])
    load(['/home/antonio/Scrivania/mhsamples/datagen_' num2str(s) '.mat'])
    load(['/home/antonio/Scrivania/sim_results/sim_results_' num2str(s) '.mat'])
    load('/home/antonio/Scrivania/mhsamples/design.mat'); design=S;
    
    for p=1:250, S_coll_means(:,p)=mean(S_coll_tot{p});end 
    Mu_post(s,1:size(S_coll_means,1)) = mean(S_coll_means,2); 
    bx = [mean(gamma,2);mean(eta);mean(delta,2)];
    Mu_true(s,1:length(bx)) = bx;
    Var_post(s,1:size(S_coll_means,1)) = var(S_coll_means,[],2);    
    Hdpi_post(s,1:size(S_coll_means,1),:) = hpdi(S_coll_means',95)';
    
    xsup=linspace(-10,10,1000); lambda = zeros(1,cell2mat(design(s,5)));
    for k=1:size(S_coll_means,1)
        [fy1] = ksdensity(S_coll_means(k,:),xsup);
        [fy0] = ksdensity(normrnd(prior_mu_b{s}{1}(k),prior_var_b{s}{1}(k),10000,1),xsup);
        %lambda(k) = trapz(min(fy1,fy0))/trapz(max(fy1,fy0));
        lambda(k) = trapz(min(fy1,fy0))/max(trapz(fy1),trapz(fy0));
    end
    
    Lambda(s,1:length(lambda)) = lambda; 
    PA_tot(s) = mean(PA_ov);
    
    %input('stop')
end


%% Plots
y=Mu_true(:,1);
x=Mu_post(:,1);
vx=Var_post(:,1)+100;
scatter(x,y,vx+1000,'filled')






