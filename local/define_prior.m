%% Function to be editable for prior specification: 
% Users can edit the body of this function to include their own density priors. 
% @ INPUT: 'x' a vector or scalar. 'K' indices to navigate x (in case 'x' is
% an array of parameters)
% @ OUTPUT: 'dx' (scalar) density evaluated at x
function dx = define_prior(x,K)

dx = sum(log(normpdf(x,0,100)));

end