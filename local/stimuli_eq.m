function b = stimuli_eq(D,x,gamma,eta,delta)

b = D*gamma + x.*(eta + D*delta);

end