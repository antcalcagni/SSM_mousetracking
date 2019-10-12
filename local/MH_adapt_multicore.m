function mh = MH_adapt_multicore(M,data,pars,maxIter,maxIter_pre,N,C_0,verbose,type_adapt,cores,typeModel,presampling,x0)

%if isempty(gcp('nocreate'))==0,delete(gcp('nocreate')); end %check for running instances of parpool


%% Presampling phase
if presampling
    disp(' ')
    Q=2; %number of parallel walkers
    C = 0.5*eye(M);
    if isempty(maxIter_pre), maxIter_pre=200; end
    expdelta=0.010; %exp decay factor: expdelta => 0 (low influence), expdelta => 1 (high influence)
    delta = [ones(1,15) exp( -linspace(15,maxIter_pre,maxIter_pre-15)/(1/expdelta) )]; %exp decay factor for vanishing cov adaptation
    delta(maxIter_pre+1)=delta(end);

    % Determine starting point by rough-X approximation
    X_rou = rough_approx_X(data.Z,pars.bnd);
    options = optimset('Display','off','algorithm','active-set');
    KK = size(data.D1,2);
    switch typeModel
        case "regression"
            x0 = fminunc(@(parsx)-observation_eq(data.Z,X_rou,data.Y,pars.a,0,0,data.x,parsx,0,pars.kappa1,pars.kappa2,pars.mu1,pars.mu2),ones(M,1),options);
            x0 = repmat(x0,1,Q);
        case "categorical"
            x0 = fminunc(@(parsx)-observation_eq(data.Z,X_rou,data.Y,pars.a,parsx,data.D1,data.x,0,zeros(KK,1),pars.kappa1,pars.kappa2,pars.mu1,pars.mu2),ones(M,1),options);
            x0 = repmat(x0,1,Q);
        case "additive"
            x0 = fminunc(@(parsx)-observation_eq(data.Z,X_rou,data.Y,pars.a,parsx(1:KK),data.D1,data.x,parsx(end),zeros(KK,1),pars.kappa1,pars.kappa2,pars.mu1,pars.mu2),ones(M,1),options);
            x0 = repmat(x0,1,Q);
        case "interaction"
            x0 = fminunc(@(parsx)-observation_eq(data.Z,X_rou,data.Y,pars.a,parsx(1:KK),data.D1,data.x,parsx(KK+1),[0;parsx(KK+2:end)],pars.kappa1,pars.kappa2,pars.mu1,pars.mu2),ones(M,1),options);
            x0 = repmat(x0,1,Q);
    end

    % Pre-sampling routine
    samples_pre = zeros(M,maxIter_pre,Q); lpdf = zeros(maxIter_pre,Q);
    samples_pre(:,1,1:Q) = x0; %initialize walkers

    for q=1:Q
        [~,lpdfv]=evaluate_sample_pre(data,pars,samples_pre(:,1,q),0,typeModel); 
        lpdf(1,q) = lpdfv;
    end

    upd=textprogressbar(maxIter_pre,' Pre-sampling phase ');
    m=1;
    while m <= maxIter_pre
        m=m+1;

        for q=1:Q        
            samples_pre(:,m,q) = gauss_rnd(samples_pre(:,m-1,q),C*delta(m));
            [res,lpdfv]=evaluate_sample_pre(data,pars,samples_pre(:,m,q),lpdf(m-1,q),typeModel); 
            lpdf(m,q)=lpdfv;
            if res==1
                samples_pre(:,m,q)=samples_pre(:,m-1,q);
                lpdf(m,q)=lpdf(m-1,q);
            end

        end

        if verbose,upd(m);end    
    end    
    x0 = mean(samples_pre(:,:),2); C_0 = topdm(cov(mean(samples_pre,3)')); %starting values and initial covariance for MH running
end

%% MH multicores
disp(' ');disp('@ MHadapt multicores started.');disp(' ');
if verbose
    disp('@@ Starting values for MH:');disp(' ')
    disp(x0');disp(' ')
    disp('@@ Initial proposal cov for MH:');disp(' ')
    disp(tril(C_0));disp(' ')
end

%parpool('local',cores);
parfor s=1:cores    
    mhs{s} = MH_adapt_singlecore(x0,C_0,maxIter,N,data,pars,true,type_adapt,typeModel);
end
%delete(gcp('nocreate')); %for matlab2017
disp(' ');disp('@ Done.')

if presampling, mh.samples_pre=samples_pre; end
mh.samples = mhs;

end


%% MH algorithm
function mh = MH_adapt_singlecore(x0,C_0,maxIter,N,data,pars,verbose,type_adapt,typeModel)

M=length(x0);

if isempty(C_0),C=eye(M)*0.20;else C=C_0;end
rej=0;

delta = [linspace(1e-4,0.2,10) ones(1,maxIter-10) 1];
%delta = [ones(1,10) ones(1,maxIter-10) 1];

samples=zeros(M,maxIter); lpdf = zeros(maxIter);

meann=zeros(M,1);
tempSample=zeros(M,N);
tempSample(:,1) = x0;
k=0;n=0;

samples(:,1) = x0;

upd=textprogressbar(maxIter,' Sampling phase ');
m=1;
while m <= maxIter
    m=m+1;
    
    samples(:,m) = gauss_rnd(samples(:,m-1),delta(m)*C);
    
    [res,lpdfv]=evaluate_sample(data,pars,samples(:,m),samples(:,m-1),lpdf(m-1),typeModel);
    lpdf(m)=lpdfv;
    if res==1
        samples(:,m)=samples(:,m-1);
        lpdf(m)=lpdf(m-1);
        rej=rej+1;
    end
    
    n=n+1;
    tempSample(:,n) = samples(:,m);
    if mod(m,N)==0
        [meann,C,k] = adaptCov(tempSample,meann,C,k,type_adapt);        
        tempSample = zeros(M,N);
        n=0;
    end
    
   if verbose,upd(m);end         
end

mh.samples=samples;
mh.x0=x0;
mh.rej=rej;
mh.acc=100*(maxIter-rej)/maxIter;


end