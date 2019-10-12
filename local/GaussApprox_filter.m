function [XF,PF,XF0,PF0] = GaussApprox_filter(Z,sigmax,a,b,bnds)

[XF,PF,XF0,PF0] = arrayfun(@(i)filter(squeeze(Z(i,:,:)),sigmax(i),a,b,bnds),1:size(Z,1),'UniformOutput', false);

XF = reshape(cell2mat(XF'),size(Z,1),size(Z,3));
PF = reshape(cell2mat(PF'),size(Z,1),size(Z,3));

end


function [xf,pf,xf0,pf0] = filter(Z,sigmax,a,b,bnds)

xf0 = zeros(1,size(Z,2)); pf0 = zeros(1,size(Z,2));
xf = zeros(1,size(Z,2)); pf = zeros(1,size(Z,2));
xf(1) = 0; pf(1) = 0.5;

for t=2:size(Z,2)        
    xf0(t) = xf(t-1);
    pf0(t) = pf(t-1)+sigmax;
    
    xf(t) = broyden(0.1,bnds,a,b,xf0(t),pf0(t),sigmax,Z(:,t));    
    %xf(t) = fsolve(@(mu)eq1(a,b,mu,xf0(t),pf0(t),sigmax,Z(:,t)),1);
    
    pf(t) = -inv(eq2(a,b,pf0(t),sigmax,xf(t))); 
    
end

end

% function [xs,ps] = smoother() %to be implemented
% 
% 
% end
