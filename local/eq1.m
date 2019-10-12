function v = eq1(aj,bj,mu,sigma,sigma_theta,theta,z)



v = sum((mu.*2.0-theta.*2.0)./(sigma.*2.0+sigma_theta.*2.0)+(aj.*z.*exp(aj.*(bj-theta))) ./ (exp(aj.*(bj-theta))+1.0)-(aj.*exp(aj.*(bj-theta)).*1.0./(exp(aj.*(bj-theta))+1.0).^2.*(z-1.0))./(1.0./(exp(aj.*(bj-theta))+1.0)-1.0));

% v = sum((mu.*2.0-theta.*2.0)./(sigma.*2.0+sigma_theta.*2.0)+(aj.*z.*exp(-bj-aj.*theta))./...
% (exp(-bj-aj.*theta)+1.0)-(aj.*exp(-bj-aj.*theta).*1.0./(exp(-bj-aj.*theta)+1.0).^2.*(z-1.0))./(1.0./(exp(-bj-aj.*theta)+1.0)-1.0));




end