clear all;
rng(1,'twister');

% Select a value for theta and data size
theta=[5;2];
N=50;

% Generate observations
error_var = 1;
g1 = ones(N,1);
g2 = [1:N]';
G = [g1 g2];
y = G*theta + sqrt(error_var)*randn(N,1);

% ML estimate for theta
theta_ML = inv(transpose(G)*G)*transpose(G)*y

% Prior distribution parameters
theta_prior = [4;2];
prior_var = [1 0 ; 0 0.1];

% Likelihood parameters
likelihood_var = error_var*inv((transpose(G)*G));

% Posterior parameters
[theta_MAP,phi,big_theta,post_var] = post_param(G,error_var,theta_prior,prior_var,y);

% Generate plots
plot(y);
text(2,y(end),...
    sprintf('Actual theta1: %4.2f \nActual theta2: %4.2f \nML theta1: %4.2f \nML theta2: %4.2f \nMAP theta1: %4.2f \nMAP theta2: %4.2f',...
    theta(1), theta(2), theta_ML(1), theta_ML(2), theta_MAP(1), theta_MAP(2)));