function [S_coll,post_means,R,neff,S] = MCMC_diag(mh,burn,thin,graph,mcmc_index,typeModel,numK)

if isempty(burn),burn=1000; end
if isempty(burn),thin=1; end
if isempty(graph),graph=1; end
if isempty(mcmc_index),mcmc_index=1; end


%% Prepare data and mixing up chains
nsample=size(mh.samples{1}.samples,2);
J=size(mh.samples{1}.samples,1);
Q=length(mh.samples); %no. of chains

S = zeros(nsample-1+1,J,Q);
for q=1:Q, S(:,:,q) = mh.samples{q}.samples(:,1:end)'; end %arrange samples into a [Nsample x J x Q] matrix (for 'mcmcdiag' package)

iid_thin=burn:thin:size(S,1); %burn-in and thinning
S = S(iid_thin,:,:);

S_coll = []; for q=1:Q, S_coll=[S_coll;S(:,:,q)]; end %mixing-up MCMC chains
post_means = mean(S_coll,1);

R=ones(1,J)*NaN;
if mcmc_index
    [R,neff,Vh,W,B,tau,thin] = psrf(S); %Gelman-Rubin Potential Scale Reduction Factor
end

%% Type of model
switch typeModel
    case "regression"
        lbls = {'\eta'};
    case "categorical"
        for j=1:numK,lbls{j} = ['\gamma' num2str(j)];end
    case "additive"
        for j=1:numK,lbls{j} = ['\gamma' num2str(j)];end
        lbls{numK+1} = '\eta';
    case "interaction"
        for j=1:numK,lbls{j} = ['\gamma' num2str(j)];end
        lbls{numK+1} = '\eta';
        u=0;for j=(numK+2):(numK*2),u=u+1;lbls{j} = ['\delta' num2str(u)];end
end


%% MCMC diagnostic plots
if graph
    pre=1;
    if length(fieldnames(mh))==1, pre=0; end    
    figure(); set(gcf,'units','points','position',[100,100,2000,700])
    u=0;
    if pre==1
        QQ=size(mh.samples_pre,3);        
        % pre-sampling phase
        for j=1:J,u=u+1;subplot(4,J,u);hold on;for q=1:QQ, if j==1, ylabel('Pre-sampling phase');end,plot(mh.samples_pre(j,:,q));end; title(lbls{j}); end
    end
    % chains
    for j=1:J,u=u+1;subplot(4,J,u); if j==1, ylabel('Chains');end, for q=1:Q, hold on;plot(S(:,j,q),'Color',[rand rand rand]); end; 
        xlabel(['R index: ' num2str(R(j))]); end    
    % marginal posteriors    
    for j=1:J,u=u+1; subplot(4,J,u); xhdi = mbe_hdi(S_coll(:,j),.95); ksdensity(S_coll(:,j)); if j==1, ylabel('Marginal Posteriors');end, vline(post_means(j),'k--'); vline(xhdi,'r-.'); 
        xlabel(['\mu: ' num2str(round(post_means(j),3)) ' \sigma^2: ' num2str(round(var(S_coll(:,j)),5))]); end
    % autocorrelation
    for j=1:J, u=u+1; subplot(4,J,u); [acfs,lags] = autocorr(S_coll(:,j)); plot(lags,acfs,'ob-'); if j==1, ylabel('ACFs');end; hline([0 0.25 0.5],'k--'); ylim([-0.1 1.1]); xlabel('lag');end
    
end



end