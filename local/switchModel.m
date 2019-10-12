function [b,gamma1,eta,delta] = switchModel(typeModel,D1,x,parx)

K=size(D1,2);

switch typeModel
    case "regression"
        gamma1 = zeros(K,1); eta = parx(1); delta = zeros(K,1);
        b = stimuli_eq(D1,x,gamma1,eta,delta);
    case "categorical"
        gamma1 = parx; eta = 0; delta = zeros(K,1);
        b = stimuli_eq(D1,x,gamma1,eta,delta);
    case "additive"
        gamma1 = parx(1:K); eta = parx(end); delta = zeros(K,1);
        b = stimuli_eq(D1,x,gamma1,eta,delta);
    case "interaction"
        gamma1 = parx(1:K); eta = parx(K+1); delta = [0;parx(K+2:end)];
        b = stimuli_eq(D1,x,gamma1,eta,delta);
end



end