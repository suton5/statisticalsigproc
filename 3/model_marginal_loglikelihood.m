function model_marg_LL = model_marginal_loglikelihood(N,P,G,error_var,theta_prior,prior_var,y)
%   Compute the marginal likelihood needed for Bayesian model selection
[theta_MAP,phi,big_theta,post_var] = post_param(G,error_var,theta_prior,prior_var,y);
exp_term = transpose(y)*y + error_var*transpose(theta_prior)*inv(prior_var)*theta_prior - transpose(big_theta)*theta_MAP;
den = ((2*pi)^(P/2))*sqrt(det(prior_var))*sqrt(det(phi))*((2*pi*error_var)^((N-P)/2));
model_marg_LL = -log(den) - exp_term/(2*error_var);
end

