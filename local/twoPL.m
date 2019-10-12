function [P] = twoPL(X,a,b)

for i=1:size(X,1)
   %P(i,:,:) = 1./(1+exp(bsxfun(@minus,-bsxfun(@times,X(i,:),a),b)));
   P(i,:,:) = 1./(1+exp(bsxfun(@times,bsxfun(@minus,X(i,:),b),-a)));
end


end