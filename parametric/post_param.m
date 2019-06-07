function [theta_MAP,phi,big_theta,post_var] = post_param(G,error_var,...
    theta_prior,prior_var,y)
%POSTERIOR_PARAM
%   Compute the parameters of the posterior
phi = transpose(G)*G + error_var*inv(prior_var);
big_theta = transpose(G)*y + error_var*inv(prior_var)*theta_prior;
theta_MAP = inv(phi)*big_theta;
post_var = error_var*inv(phi);
end