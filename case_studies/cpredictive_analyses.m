
%% Set environment and load data
addpath(genpath('SSM_MouseTracking'))
clear all
load('barca2012_nomissin_full.mat')
ncores=feature('numcores');
nsamples=1000;


%% Models

disp(' ');
disp('@@@ RUN MODEL 1_0 @@@')

q=length([stm_HF' stm_LF']);
data.Z=Z(:,1:q,:); data.Y=Y(:,1:q,:); data.D1=D1(1:q,1:2); pars.sigmax=sigmax; pars.kappa1=kappa1; pars.kappa2=kappa2; pars.a=a(1:q); 
pars.mu1=mu1; pars.mu2=mu2; pars.bnd=bnd; data.x=scaledata(AoA_cov,-3,3);

load('mh_1_0.mat')
S_coll = MCMC_diag(mh_1_0,burnin,1,false,false,'regression',0);
iid = round(unifrnd(1,size(S_coll,1),nsamples,1));
dic_1_0 = compute_dic(data,pars,S_coll(iid,:),ncores,'regression');
save('dic_1_0.mat','dic_1_0')

disp(' ');
disp('@@@ RUN MODEL 1_1 @@@')
load('mh_1_1.mat')
S_coll = MCMC_diag(mh_1_1,burnin,1,false,false,'categorical',2);
iid = round(unifrnd(1,size(S_coll,1),nsamples,1));
dic_1_1 = compute_dic(data,pars,S_coll(iid,:),ncores,'categorical');
save('dic_1_1.mat','dic_1_1')

disp(' ');
disp('@@@ RUN MODEL 1_2 @@@')
load('mh_1_2.mat')
S_coll = MCMC_diag(mh_1_2,burnin,1,false,false,'interaction',4);
iid = round(unifrnd(1,size(S_coll,1),nsamples,1));
dic_1_2 = compute_dic(data,pars,S_coll(iid,:),ncores,'interaction');
save('dic_1_2.mat','dic_1_2')


disp(' ');
disp('@@@ RUN MODEL 2_0 @@@')

data.x=scaledata(IMM_cov,-3,3);

load('mh_2_0.mat')
S_coll = MCMC_diag(mh_2_0,burnin,1,false,false,'regression',0);
iid = round(unifrnd(1,size(S_coll,1),nsamples,1));
dic_2_0 = compute_dic(data,pars,S_coll(iid,:),ncores,'regression');
save('dic_2_0.mat','dic_2_0')

disp(' ');
disp('@@@ RUN MODEL 2_2 @@@')
load('mh_2_2.mat')
S_coll = MCMC_diag(mh_2_2,burnin,1,false,false,'interaction',4);
iid = round(unifrnd(1,size(S_coll,1),nsamples,1));
dic_2_2 = compute_dic(data,pars,S_coll(iid,:),ncores,'interaction');
save('dic_2_2.mat','dic_2_2')


disp(' ');
disp('@@@ RUN MODEL 3_0 @@@')

q0=length([stm_HF' stm_LF'])+1; q=q0+length([stm_NW' stm_PW'])-1;
data.Z=Z(:,q0:q,:); data.Y=Y(:,q0:q,:); data.D1=D1(q0:q,3:4); pars.sigmax=sigmax; pars.kappa1=kappa1; pars.kappa2=kappa2; pars.a=a(q0:q); 
pars.mu1=mu1; pars.mu2=mu2; pars.bnd=bnd; data.x=scaledata(bigF_cov(q0:q),-3,3);

load('mh_3_0.mat')
S_coll = MCMC_diag(mh_3_0,burnin,1,false,false,'regression',0);
iid = round(unifrnd(1,size(S_coll,1),nsamples,1));
dic_3_0 = compute_dic(data,pars,S_coll(iid,:),ncores,'regression');
save('dic_3_0.mat','dic_3_0')

disp(' ');
disp('@@@ RUN MODEL 3_1 @@@')
load('mh_3_1.mat')
S_coll = MCMC_diag(mh_3_1,burnin,1,false,false,'categorical',2);
iid = round(unifrnd(1,size(S_coll,1),nsamples,1));
dic_3_1 = compute_dic(data,pars,S_coll(iid,:),ncores,'categorical');
save('dic_3_1.mat','dic_3_1')

disp(' ');
disp('@@@ RUN MODEL 3_2 @@@')
load('mh_3_2.mat')
S_coll = MCMC_diag(mh_3_2,burnin,1,false,false,'interaction',4);
iid = round(unifrnd(1,size(S_coll,1),nsamples,1));
dic_3_2 = compute_dic(data,pars,S_coll(iid,:),ncores,'interaction');
save('dic_3_2.mat','dic_3_2')

disp(' ');
disp('@@@ RUN MODEL 4_0 @@@')

data.Z=Z; data.Y=Y; data.D1=D1; pars.sigmax=sigmax; pars.kappa1=kappa1; pars.kappa2=kappa2; pars.a=a; 
pars.mu1=mu1; pars.mu2=mu2; pars.bnd=bnd; data.x=scaledata(bigF_cov,-3,3);

load('mh_4_0.mat')
S_coll = MCMC_diag(mh_4_0,burnin,1,false,false,'regression',0);
iid = round(unifrnd(1,size(S_coll,1),nsamples,1));
dic_4_0 = compute_dic(data,pars,S_coll(iid,:),ncores,'regression');
save('dic_4_0.mat','dic_4_0')

disp(' ');
disp('@@@ RUN MODEL 4_1 @@@')
load('mh_4_1.mat')
S_coll = MCMC_diag(mh_4_1,burnin,1,false,false,'categorical',4);
iid = round(unifrnd(1,size(S_coll,1),nsamples,1));
dic_4_1 = compute_dic(data,pars,S_coll(iid,:),ncores,'categorical');
save('dic_4_1.mat','dic_4_1')

disp(' ');
disp('@@@ RUN MODEL 4_2 @@@')
load('mh_4_2.mat')
S_coll = MCMC_diag(mh_4_2,burnin,1,false,false,'interaction',8);
iid = round(unifrnd(1,size(S_coll,1),nsamples,1));
dic_4_2 = compute_dic(data,pars,S_coll(iid,:),ncores,'categorical');
save('dic_4_2.mat','dic_4_2')

