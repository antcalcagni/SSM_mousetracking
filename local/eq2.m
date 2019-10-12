function v = eq2(aj,bj,sigma,sigma_theta,theta)

v = sum(aj.^2.*1.0./cosh(aj.*(bj-theta).*(1.0./2.0)).^2.*(-1.0./4.0)-1.0./(sigma+sigma_theta));
%to get expression without cosh, type: rewrite(aj.^2.*1.0./cosh(aj.*(bj-theta).*(1.0./2.0)).^2.*(-1.0./4.0)-1.0./(sigma+sigma_theta),'exp')

%v = sum(aj.^2.*1.0./cosh(bj.*(1.0./2.0)+aj.*theta.*(1.0./2.0)).^2.*(-1.0./4.0)-1.0./(sigma+sigma_theta));




end