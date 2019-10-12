%___________________________________________________________________________________________________________________________________________
% Random Generation from the von Mises Distribution
% Simulation from the von Mises distribution is done via the algorithm due to Best and Fisher (1979).
% 
% INPUT
% n: number of random variables to generate.
% mean: mean direction in radians of the von Mises distribution.
% k: concentration parameter of the von Mises distribution.
%
% OUTPUT
% Returns a vector of n independent random variables generated from a von Mises distribution. 
% Values are between 0 and 2 pi.
%
% REFERENCES
% Best, D. and Fisher, N. (1979). Efficient simulation of the von Mises distribution. Applied Statistics, 24, 152-157.
%
% INFO
% From R-package: CircStats version 0.2-4
%___________________________________________________________________________________________________________________________________________

function [theta] = rvm(n,mean,k) 

    theta = 1:n;
    a = 1 + (1 + 4 * (k^2))^0.5;
    b = (a - (2 * a)^0.5)/(2 * k);
    r = (1 + b^2)/(2 * b);
    obs = 1;
    
    while obs <= n
        U1 = unifrnd(0, 1);
        z = cos(pi * U1);
        f = (1 + r * z)/(r + z);
        c = k * (r - f);
        U2 = unifrnd(0,1);
        if c * (2 - c) - U2 > 0
            U3 = unifrnd(0,1);
            theta(obs) = sign(U3 - 0.5) * acos(f) + mean;
            theta(obs) = mod(theta(obs),(2 * pi));
            obs = obs + 1;
        elseif log(c/U2) + 1 - c >= 0
                U3 = unifrnd(0, 1);
                theta(obs) = sign(U3 - 0.5) * acos(f) + mean;
                theta(obs) = mod(theta(obs),(2 * pi));
                obs = obs + 1;
        end
    end
end
