function [datagen_agg,B,eta,gamma,delta] = generateModels(ncores,M,T,I,J,K,mu1,mu2,kappa1,kappa2,typeModel,prior_mu_b,prior_var_b,x,D1,sigmax,bnd)

thr=(mu1+mu2)/2; a=unifrnd(1,1,J,1);

disp(['@@ Generate data for the model: ',typeModel])
a=unifrnd(1,1,J,1);

switch typeModel
    case "regression"
        gamma=0;delta=0;
        eta = gauss_rnd(prior_mu_b,diag(prior_var_b),M);
        B = x.*eta;
    case "categorical"
        eta=0;delta=0;
        gamma = gauss_rnd(prior_mu_b,diag(prior_var_b),M);
        B = D1*gamma;
    case "interaction"
        gamma = gauss_rnd(prior_mu_b(1:K),diag(prior_var_b(1:K)),M);
        eta = gauss_rnd(prior_mu_b(K+1),diag(prior_var_b(K+1)),M);
        delta = gauss_rnd(prior_mu_b(K+2:end),diag(prior_var_b(K+2:end)),M);
        B = D1*gamma + x.*(eta + D1*[zeros(1,M);delta]);
end

% split data for easy parallelization
if mod(M/ncores,1)==0, iid = repmat(M/ncores,1,ncores);else iid = repmat(floor(M/ncores),1,ncores-1); iid = [iid M-sum(iid)]; end
disp(['@@ Datasets by core: ',num2str(iid)])
B_split = mat2cell(B,J,iid);

parfor q=1:ncores, datagen{q} = arrayfun(@(m)generateData(T,I,J,mu1,mu2,kappa1,kappa2,B_split{q}(:,m),a,thr,sigmax,bnd),1:iid(q),'UniformOutput',false); end        

datagen_agg=cell(0); for q=1:ncores, datagen_agg = [datagen_agg datagen{q}]; end

      
end
