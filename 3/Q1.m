clear all;
rng(1,'twister');

% Select a value for theta and data size
theta=5;
N=5;

% Generate observations
error_var = 1;
G = ones(N,1);
y = G*theta + sqrt(error_var)*randn(N,1);

% ML estimate for theta
theta_ML = inv(transpose(G)*G)*transpose(G)*y

% Prior distribution parameters
theta_prior = 4
prior_var = 0.5;

% Likelihood parameters
likelihood_var = error_var/(transpose(G)*G)

% Posterior parameters
[theta_MAP,phi,big_theta,post_var] = post_param(G,error_var,theta_prior,prior_var,y);

% Posterior distribution
x = [2:0.001:8];
y_post = normpdf(x,theta_MAP,sqrt(post_var));
y_likelihood = normpdf(x,theta_ML,sqrt(likelihood_var));
y_prior = normpdf(x,theta_prior, sqrt(prior_var));
plot(x, y_post) 
hold on
plot(x, y_likelihood)
plot(x, y_prior)
xline(theta_ML,'-.','ML')
xline(theta_MAP,'-.','MAP')
legend('Posterior','Likelihood','Prior')
hold off