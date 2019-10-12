%___________________________________________________________________________________________________________________________________________
% von Mises Density Function
% Returns the von Mises density function evaluated at a particular value.
% 
% INPUT
% theta: vector of angles measured in radians at which to evaluate the density function.
% mu: mean direction in radians of the von Mises distribution.
% kappa: (non-negative) concentration parameter of the von Mises distribution.
%
% OUTPUT
% Returns the von Mises density function evaluated at theta.
%
% REFERENCES
% Jammalamadaka, S. Rao and SenGupta, A. (2001). Topics in Circular Statistics, Section 2.2.4, World Scientific Press, Singapore.
%
% INFO
% From R-package: CircStats version 0.2-4
%___________________________________________________________________________________________________________________________________________

function [dtheta] = dvm(theta,mu,kappa) 
    
C = exp(-kappa)*besseli(0,kappa); %exponential-scaled besselI function
dtheta = (1/(2 *pi*C)) .* (exp(cos(theta - mu) - 1)).^kappa;
        
end