clear all;
rng(2,'twister');

% Select a value for theta and data size
theta=5;
N=100;

% Generate observations
error_var = 10;
G = ones(N,1);
y = G*theta + sqrt(error_var)*randn(N,1);
% figure(1)
% plot(y)
% hold on
% yline(theta, 'r', 'theta')
% hold off
% xlabel('n')
% ylabel('y_n')
% title(['Example data with theta=5, errorvar=1, N=', num2str(N)])

% ML estimate for theta
theta_ML = inv(transpose(G)*G)*transpose(G)*y

% Prior distribution parameters
theta_prior = 5.2
prior_var = 0.2;

% Likelihood parameters
likelihood_var = error_var/(transpose(G)*G)

% Posterior parameters
[theta_MAP,phi,big_theta,post_var] = post_param(G,error_var, ...
    theta_prior,prior_var,y);

% Posterior distribution
x = [2:0.001:8];
y_post = normpdf(x,theta_MAP,sqrt(post_var));
y_likelihood = normpdf(x,theta_ML,sqrt(likelihood_var));
y_prior = normpdf(x,theta_prior, sqrt(prior_var));
figure(2)
plot(x, y_post) 
hold on
plot(x, y_likelihood)
plot(x, y_prior)
xline(theta_ML,'-.','ML','LabelHorizontalAlignment','left')
xline(theta_MAP,'-.','MAP','LabelHorizontalAlignment','left')
legend('Posterior','Likelihood','Prior')
hold off
xlabel('theta')
ylabel('Probability Density')
title(sprintf('The parameter distributions with theta=5, for N=%s, errorvar=%s, priormean=%s, priorvar=%s', num2str(N), num2str(error_var), num2str(theta_prior), num2str(prior_var)))