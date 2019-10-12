%___________________________________________________________________________________________________________________________________________
% Random Generation from the mixed von Mises Distribution
% Simulation from the von Mises distribution is done via the algorithm due to Best and Fisher (1979).
% 
% INPUT
% n: number of random variables to generate.
% mu1: mean direction in radians of the first von Mises distribution.
% mu2: mean direction in radians of the second von Mises distribution.
% kappa1: concentration parameter of the first von Mises distribution.
% kappa2: concentration parameter of the first von Mises distribution.
%
% OUTPUT
% Returns a vector of n independent random variables generated from a mixture of two von Mises distributions.
% Values are between 0 and 2 pi.
%
% REFERENCES
% Best, D. and Fisher, N. (1979). Efficient simulation of the von Mises distribution. Applied Statistics, 24, 152-157.
%
% INFO
% From R-package: CircStats version 0.2-4
%___________________________________________________________________________________________________________________________________________


function [theta] = rmixedvm(n,mu1,mu2,kappa1,kappa2,p)

for i=1:n
    test = unifrnd(0,1);
    if test < p, theta(i) = rvm(1,mu1,kappa1);
    else theta(i) = rvm(1,mu2,kappa2);
    end
end