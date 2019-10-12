function [rej,lpdf] = evaluate_sample_pre(data,pars,xcur,lpdf_pre,typeModel)

rej=0;

[b,gamma1,eta,delta] = switchModel(typeModel,data.D1,data.x,xcur);

X = GaussApprox_filter(data.Z,pars.sigmax,pars.a,b,0);
lpdf = observation_eq(data.Z,X,data.Y,pars.a,gamma1,data.D1,data.x,eta,delta,pars.kappa1,pars.kappa2,pars.mu1,pars.mu2);
%lpdf = compute_marginal_likelihood(data,pars,gamma1,delta,eta,2);

rho = lpdf - lpdf_pre; 
    
if rho < log(rand) %reject
    rej=1;
end


end