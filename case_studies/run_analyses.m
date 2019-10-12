
%% Set environment and load data
addpath(genpath('SSM_MouseTracking'))
clear all
load('barca2012_nomissin_full.mat')

%% Run models

maxIterPre=50; maxIterMH=200;
ncores=feature('numcores');

disp(' ');
disp('@@@ RUN MODEL 1_0 @@@')

q=length([stm_HF' stm_LF']);
data.Z=Z(:,1:q,:); data.Y=Y(:,1:q,:); data.D1=D1(1:q,1:2); pars.sigmax=sigmax; pars.kappa1=kappa1; pars.kappa2=kappa2; pars.a=a(1:q); 
pars.mu1=mu1; pars.mu2=mu2; pars.bnd=bnd; data.x=scaledata(AoA_cov,-3,3);

[mh_1_0] = MH_adapt_multicore(1,data,pars,maxIterMH,maxIterPre,25,[],true,'type1',ncores,'regression',true,[]);
save('mh_1_0.mat','mh_1_0')

disp(' ');
disp('@@@ RUN MODEL 1_1 @@@')

[mh_1_1] = MH_adapt_multicore(2,data,pars,maxIterMH,maxIterPre,25,[],true,'type1',ncores,'categorical',true,[]);
save('mh_1_1.mat','mh_1_1')

disp(' ');
disp('@@@ RUN MODEL 1_2 @@@')

[mh_1_2] = MH_adapt_multicore(4,data,pars,maxIterMH,maxIterPre,25,[],true,'type1',ncores,'interaction',true,[]);
save('mh_1_2.mat','mh_1_2')

disp(' ');
disp('@@@ RUN MODEL 2_0 @@@')

data.x=scaledata(IMM_cov,-3,3);

[mh_2_0] = MH_adapt_multicore(1,data,pars,maxIterMH,maxIterPre,25,[],true,'type1',ncores,'regression',true,[]);
save('mh_2_0.mat','mh_2_0')

disp(' ');
disp('@@@ RUN MODEL 2_2 @@@')

[mh_2_2] = MH_adapt_multicore(4,data,pars,maxIterMH,maxIterPre,25,[],true,'type1',ncores,'interaction',true,[]);
save('mh_2_2.mat','mh_2_2')

disp(' ');
disp('@@@ RUN MODEL 3_0 @@@')

q0=length([stm_HF' stm_LF'])+1; q=q0+length([stm_NW' stm_PW'])-1;
data.Z=Z(:,q0:q,:); data.Y=Y(:,q0:q,:); data.D1=D1(q0:q,3:4); pars.sigmax=sigmax; pars.kappa1=kappa1; pars.kappa2=kappa2; pars.a=a(q0:q); 
pars.mu1=mu1; pars.mu2=mu2; pars.bnd=bnd; data.x=scaledata(bigF_cov(q0:q),-3,3);

[mh_3_0] = MH_adapt_multicore(1,data,pars,maxIterMH,maxIterPre,25,[],true,'type1',ncores,'regression',true,[]);
save('mh_3_0.mat','mh_3_0')

disp(' ');
disp('@@@ RUN MODEL 3_1 @@@')

[mh_3_1] = MH_adapt_multicore(2,data,pars,maxIterMH,maxIterPre,25,[],true,'type1',ncores,'categorical',true,[]);
save('mh_3_1.mat','mh_3_1')

disp(' ');
disp('@@@ RUN MODEL 3_2 @@@')

[mh_3_2] = MH_adapt_multicore(4,data,pars,maxIterMH,maxIterPre,25,[],true,'type1',ncores,'interaction',true,[]);
save('mh_3_2.mat','mh_3_2')


disp(' ');
disp('@@@ RUN MODEL 4_0 @@@')

data.Z=Z; data.Y=Y; data.D1=D1; pars.sigmax=sigmax; pars.kappa1=kappa1; pars.kappa2=kappa2; pars.a=a; 
pars.mu1=mu1; pars.mu2=mu2; pars.bnd=bnd; data.x=scaledata(bigF_cov,-3,3);

[mh_4_0] = MH_adapt_multicore(1,data,pars,maxIterMH,maxIterPre,25,[],true,'type1',ncores,'regression',true,[]);
save('mh_4_0.mat','mh_4_0')

disp(' ');
disp('@@@ RUN MODEL 4_1 @@@')

[mh_4_1] = MH_adapt_multicore(4,data,pars,maxIterMH,maxIterPre,25,[],true,'type1',ncores,'categorical',true,[]);
save('mh_4_1.mat','mh_4_1')

disp(' ');
disp('@@@ RUN MODEL 4_2 @@@')

[mh_4_2] = MH_adapt_multicore(8,data,pars,maxIterMH,maxIterPre,25,[],true,'type1',ncores,'interaction',true,[]);
save('mh_4_2.mat','mh_4_2')

%% New models

q0=1;q=72; %select just {HF,LF,NW}
data.Z=Z(:,q0:q,:); data.Y=Y(:,q0:q,:); data.D1=D1(q0:q,1:3); pars.sigmax=sigmax; pars.kappa1=kappa1; pars.kappa2=kappa2; pars.a=a(q0:q); 
pars.mu1=mu1; pars.mu2=mu2; pars.bnd=bnd; data.x=scaledata(bigF_cov(q0:q),-3,3);

maxIterPre=100; maxIterMH=500; ncores=feature('numcores');

disp(' ');
disp('@@@ RUN MODEL 5_0 @@@')

q=length([stm_HF' stm_LF']);
data.Z=Z(:,1:q,:); data.Y=Y(:,1:q,:); data.D1=D1(1:q,1:2); pars.sigmax=sigmax; pars.kappa1=kappa1; pars.kappa2=kappa2; pars.a=a(1:q); 
pars.mu1=mu1; pars.mu2=mu2; pars.bnd=bnd; data.x=scaledata(AoA_cov,-3,3);

[mh_5_0] = MH_adapt_multicore(1,data,pars,maxIterMH,maxIterPre,25,[],true,'type1',ncores,'regression',true,[]);
save('mh_5_0.mat','mh_5_0')

[mh_5_1] = MH_adapt_multicore(3,data,pars,maxIterMH,maxIterPre,25,[],true,'type1',ncores,'categorical',true,[]);
save('mh_5_1.mat','mh_5_1')

[mh_5_2] = MH_adapt_multicore(6,data,pars,maxIterMH,maxIterPre,25,[],true,'type1',ncores,'interaction',true,[]);
save('mh_5_2.mat','mh_5_2')







