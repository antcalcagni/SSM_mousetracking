function xval = compute_marginal_likelihood(data,pars,gamma,delta,eta,ncores)

[I,J,T] = size(data.Z);

% define densities
g = @(x,a,b) 1./(1+exp(-a.*(x-b))); fvm = @(y,k,mu) exp(cos(y-mu).^k);
fmarg = @(y,x,z,b,k1,k2,mu1,mu2,xf,pf) (((g(x,1,b).*fvm(y,k1,mu1)).^z) .* (((1-g(x,1,b)).*fvm(y,k2,mu2)).^(1-z))).*(exp(-0.5*(((x-xf).^2)./pf.^2)));

% split data for easy parallelization
if mod(I/ncores,1)==0, iid = repmat(I/ncores,1,ncores);else iid = [repmat(floor(I/ncores),1,ncores-1) floor(I/ncores)+1]; end
Y_split = mat2cell(data.Y,iid,J,T); Z_split = mat2cell(data.Z,iid,J,T); 
sigmax_split = mat2cell(pars.sigmax,iid,1);

% compute marginal likelihood
parfor k=1:size(Z_split,1)    
    mhs{k} = computeLikeMarg_singlecore(Z_split{k},Y_split{k},data,pars.a,gamma,delta,eta,pars.kappa1,pars.kappa2,pars.mu1,pars.mu2,sigmax_split{k},pars.bnd,fmarg)
end
xval = sum(cell2mat(mhs));

end

function xval = computeLikeMarg_singlecore(Z,Y ,data,a,gamma,delta,eta,k1,k2,mu1,mu2,sigmax,bnd,fmarg)
[I,J,T] = size(Z);

b = stimuli_eq(data.D1,data.x,gamma,eta,delta);

[XF,PF] = GaussApprox_filter(Z,sigmax,a,b,0);
XF=[zeros(I,1) XF(:,1:(T-1))]; PF=[PF(:,1) PF(:,1:T-1)];

xval = sum(arrayfun(@(i)...
    -sum(log(integral(@(x)...
    fmarg(cw_reshape(squeeze(Y(i,:,:))),x,cw_reshape(squeeze(Z(i,:,:))),cw_reshape(kron(ones(T,1),b')),k1,k2,mu1,mu2,repmat(XF(i,:)',J,1),repmat(PF(i,:)',J,1)),...
    -bnd,bnd,'ArrayValued',true))),1:I));

end


function Y = cw_reshape(X)
Y=[]; for j=1:size(X,2),Y=[Y;squeeze(X(:,j,:))];end
end