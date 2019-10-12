%___________________________________________________________________________________________________________________________________________
% von Mises Maximum Likelihood Estimates
% Computes the maximum likelihood estimates for the parameters of a von
% Mises distribution.
% 
% INPUT
% theta: vector of angular measurements in radians.
% bias: logical flag: if 1, the estimate for kappa is computed with a bias corrected method.
%
% OUTPUT
% Returns a data frame with two components: mu and kappa, the MLES of the respective parameters.
%
% REFERENCES
% Jammalamadaka, S. Rao and SenGupta, A. (2001). Topics in Circular Statistics, Section 4.2.1, World Scientific Press, Singapore.
%
% INFO
% From R-package: CircStats version 0.2-4
%___________________________________________________________________________________________________________________________________________

function [kappa] = est_kappa(theta,bias) 

if isempty(bias),bias=0; end

    mu = circ_mean(theta);
    kappa = A1inv(mean(cos(theta - mu))); %estimate kappa without bias correction
    
    if bias == 1 %estimate kappa with bias correction
        kappa_ml = kappa;
        n = length(theta);
        if (kappa_ml < 2), kappa = max(kappa_ml - 2 * (n * kappa_ml)^-1, 0); end
        if (kappa_ml >= 2), kappa = ((n - 1)^3 * kappa_ml)/(n^3 + n); end
    end
    
end



