function X_rou = rough_approx_X(Z,bnd)

I=size(Z,1); T=size(Z,3); J=size(Z,2);

W=zeros(I,T); X_rou=zeros(I,T);

for i=1:I, W(i,:) = sum(squeeze(Z(i,:,:)),1); end

% Just pertube the original Z data in order to fit a smoothing model
if sum(sum(sum(W)))<=1
    for k=1:3, Z(round(unifrnd(1,I)),round(unifrnd(1,J)),round(unifrnd(1,T)))=1; end
end

W(:,2:T) = scaledata(W(:,2:T),-bnd,bnd);

xdom=2:1:T;
for i=1:I
    fitobj = fit(xdom',W(i,2:T)','smoothingspline','SmoothingParam',0.55);
    X_rou(i,2:T) = feval(fitobj,xdom);
%     fitobj = splinefit(xdom,W(i,2:T),20);
%     X_rou(i,2:T) = ppval(fitobj,xdom);
end

end