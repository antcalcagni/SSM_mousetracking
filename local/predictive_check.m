function [Y_star,pred_stats] = predictive_check(S_coll,data,pars,ncores,typeModel)

[M,K] = size(S_coll);

if mod(M/ncores,1)==0, iid = repmat(M/ncores,1,ncores);else iid = repmat(floor(M/ncores),1,ncores-1); iid = [iid M-sum(iid)]; end
thetapost_split = mat2cell(S_coll,iid,K);

parfor q=1:ncores
    [Y_star{q,1},pred_stats{q}] = generate_data(data,pars,thetapost_split{q},typeModel);    
end

Y_star = cat(1,Y_star{:});
pred_stats = cat(1,pred_stats{:});
end


function [Y,stats] = generate_data(data,pars,thetapost,typeModel)
Y=cell(size(thetapost,1),1); stats=zeros(size(thetapost,1),6);

for i=1:size(thetapost,1)
    b = switchModel(typeModel,data.D1,data.x,thetapost(i,:)');
    XF = GaussApprox_filter(data.Z,pars.sigmax,pars.a,b,pars.bnd);
    P = twoPL(XF,pars.a,b);p=P(:);
    y = arrayfun(@(j)rmixedvm(1,pars.mu1,pars.mu2,pars.kappa1,pars.kappa2,p(j)),1:length(p));    
    
    Y{i}=reshape(y,size(data.Z,1),size(data.Z,2),size(data.Z,3));
    stats(i,:) = compute_stats(Y{i},pars);
end

end

function stats = compute_stats(Y,pars)

y=Y(:); z=ones(length(y),1);
z(y<(pars.mu1+pars.mu2)*0.5)=0;

stats(1)=circ_mean(y(z==1));  stats(2)=circ_mean(y(z==0));
stats(3)=sum(z)/length(z); stats(4)=1-(sum(z)/length(z));
stats(5)=est_kappa(y(z==1),1); stats(6)=est_kappa(y(z==0),1);


end
