%% Set environment
addpath(genpath('SSM_MouseTracking/'))
clear all
close all

%% Generate true model
T=61; I=2; J=18; rescaled=1; bnd=5;
X=zeros(I,T); P=zeros(I,J,T); Y=zeros(I,J,T); Z=zeros(I,J,T);

mu1=2.755; mu2=0.755; kappa1=200; kappa2=200; 
thr=(2.755+0.755)/2;
sigmax = repmat(5.5,[1 I]);
a=unifrnd(1,1,J,1);

K=3; %number of levels of D1
D1 = kron(eye(K),ones(J/K,1)); %generate a well-distribuited classification matrix

gamma1=unifrnd(-2,2,K,1);
delta=unifrnd(1,3,K-1,1)*0;

x = sort(unifrnd(-2,3,J,1));
eta = 1.4*0;

b = stimuli_eq(D1,x,gamma1,eta,[0;delta*0]);

X(:,1) = 0;     
P(:,:,1) = twoPL(X(:,1),a,b);
for i=1:I
    for j=1:J
        Y(i,j,1) = rmixedvm(1,mu1,mu2,kappa1,kappa2,P(i,j,1));
        if Y(i,j,1)>=thr, Z(i,j,1)=1; end        
    end
end

for i=1:I
    for t=2:T  
        X(i,t) = X(i,t-1) + gauss_rnd(0,sigmax(i));
    end
end

if rescaled==1 
    for i=1:I
        X(i,2:T)=scaledata(X(i,2:T),-bnd,bnd);
    end
end

for i=1:I
    for j=1:J
        for t=2:T
            P(i,j,t) = twoPL(X(i,t),a(j),b(j));
            Y(i,j,t) = rmixedvm(1,mu1,mu2,kappa1,kappa2,P(i,j,t));
            if Y(i,j,t)>=thr, Z(i,j,t)=1; end        
        end
    end
end

data.Z=Z; data.Y=Y; data.D1=D1; data.x=x;
pars.sigmax=sigmax'; pars.kappa1=kappa1; pars.kappa2=kappa2; pars.a=a; 
pars.mu1=mu1; pars.mu2=mu2; pars.bnd=bnd; 
pars.gamma1=gamma1; pars.delta=delta; pars.eta=eta;


%% Filtering
tic
bnds = 0;
[XF,PF] = GaussApprox_filter(Z,sigmax,a,b,bnds);
toc

%% Rough approximation of X
X_rou = rough_approx_X(Z,bnd);

for i=1:I, subplot(5,6,i); hold on;     
    plot(X(i,:)); 
    plot(X_rou(i,:),'r-'); 
    plot(XF(i,:),'g-'); 
end


%% MH adapt [new algorithm]
w=[gamma1];
[mh] = MH_adapt_multicore(length(w),data,pars,3000,200,25,[],true,'type1',1,'categorical',false,w+2)


%% Preliminary plot for MH runnig
for j=1:length(w),subplot(2,3,j);hold on;for q=1:2,plot(mh.samples_pre(j,:,q));end;hline(w(j),'k');end


%% Diagnostics MCMC
addpath(genpath('/home/antonio/MEGA/Lavoro_sync/MATLAB/mcmcdiag/'))
close all
clear S S_coll

nsample=size(mh.samples{1}.samples,2);
J=size(mh.samples{1}.samples,1);
Q=length(mh.samples);

S = zeros(nsample-1+1,J,Q);
for q=1:Q, S(:,:,q) = mh.samples{q}.samples(:,1:end)'; end %arrange samples into a [Nsample x J x Q] matrix (for 'mcmcdiag' package)

burn=1000; thin=1; %burning and thinning
iid_thin=burn:thin:size(S,1);
S = S(iid_thin,:,:);

S_coll = []; for q=1:Q, S_coll=[S_coll;S(:,:,q)]; end %mixing-up MCMC chains

[R,neff,Vh,W,B,tau,thin] = psrf(S); %Gelman-Rubin Potential Scale Reduction Factor

disp(' ');
disp(['R index: ' num2str(R)])
disp(['Neff: ' num2str(neff)])
disp(['post mean: ' num2str(mean(S_coll,1))])

figure();for j=1:J, subplot(1,J,j); for q=1:Q, hold on;plot(S(:,j,q),'Color',[rand rand rand]); title('\gamma'); end; end    
figure();for j=1:J, subplot(1,J,j); hist(S_coll(:,j)); title('\gamma'); end

figure();u=0;for q=1:Q, for j=1:J, u=u+1; subplot(Q,J,u); autocorr(S(:,j,q)); title('\gamma'); end; end







