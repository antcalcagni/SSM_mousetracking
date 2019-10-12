function res = compute_dic(data,pars,thetapost,ncores,typeModel)

disp('');disp('@ DIC routine started')

K=size(data.D1,2);
[Q,M]=size(thetapost);

if mod(Q/ncores,1)==0, iid = repmat(Q/ncores,1,ncores);else iid = repmat(floor(Q/ncores),1,ncores-1); iid = [iid Q-sum(iid)]; end
thetapost_split = mat2cell(thetapost,iid,M);

disp('@@ Computing correction term p..')
llSum = zeros(length(iid),1);

switch typeModel
    case 'regression'
        parfor k=1:ncores,llSum(k) = sum(arrayfun(@(q)compute_marginal_likelihood(data,pars,zeros(K,1),zeros(K,1),thetapost_split{k}(q,1)',ncores),1:iid(k)));end
        
    case 'categorical'
        parfor k=1:ncores,llSum(k) = sum(arrayfun(@(q)compute_marginal_likelihood(data,pars,thetapost_split{k}(q,1:K)',zeros(K,1),0,ncores),1:iid(k)));end
        
    case 'additive'        
        parfor k=1:ncores,llSum(k) = sum(arrayfun(@(q)compute_marginal_likelihood(data,pars,thetapost_split{k}(q,1:K)',zeros(K,1),thetapost_split{k}(q,end)',ncores),1:iid(k)));end
       
    case 'interaction'
        parfor k=1:ncores,llSum(k) = sum(arrayfun(@(q)compute_marginal_likelihood(data,pars,thetapost_split{k}(q,1:K)',[0;thetapost_split{k}(q,K+2:end)'],thetapost_split{k}(q,K+1)',ncores),1:iid(k)));end
end


disp('@@ Computing final DIC ..')
thetapost_mean = mean(thetapost);

switch typeModel
    case 'regression'
        l0 = compute_marginal_likelihood(data,pars,zeros(K,1),zeros(K,1),thetapost_mean(1)',ncores);
    
    case 'categorical'
        l0 = compute_marginal_likelihood(data,pars,thetapost_mean(1:K)',zeros(K,1),0,ncores);
        
    case 'additive'        
        l0 = compute_marginal_likelihood(data,pars,thetapost_mean(1:K)',zeros(K,1),thetapost_mean(end)',ncores);
        
    case 'interaction'
        l0 = compute_marginal_likelihood(data,pars,thetapost_mean(1:K)',[0;thetapost_mean(K+2:end)'],thetapost_mean(K+1)',ncores);
end

p0 = 2*(l0-(1/Q*sum(llSum)));
dic = -2*(l0-p0);
disp(['@@ DIC: ' num2str(dic)])

res.llSum=llSum;
res.l0=l0;
res.p0=p0;
res.dic=dic;

end