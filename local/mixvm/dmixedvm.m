%___________________________________________________________________________________________________________________________________________
% Mixture of von Mises Density Function
% Returns the mixture of von Mises density function evaluated at a particular value.
% 
% INPUT
% theta: vector of angles measured in radians at which to evaluate the density function.
% mu1: mean direction in radians of the first von Mises distribution.
% mu2: mean direction in radians of the second von Mises distribution.
% kappa1: (non-negative) concentration parameter of the first von Mises distribution.
% kappa2: (non-negative) concentration parameter of the second von Mises distribution.
% p: mixing proportion between 0 and 1.
%
% OUTPUT
% Evaluates the density function of a mixture of two von Mises distributions at a given vector of angles.
%
% INFO
% From R-package: CircStats version 0.2-4
%___________________________________________________________________________________________________________________________________________


function [dtheta] = dmixedvm(theta,mu1,mu2,kappa1,kappa2,p) 

C1 = exp(-kappa1)*besseli(0,kappa1); %exponential-scaled besselI function
C2 = exp(-kappa2)*besseli(0,kappa2); 

dtheta = p/(2*pi*C1) .* (exp(cos(theta - mu1) - 1)).^kappa1 + (1 - p)/(2 *pi*C2) .* (exp(cos(theta - mu2) - 1)).^kappa2;
end