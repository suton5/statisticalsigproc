clear all;

% Select a value for theta and data size
theta2=5;
theta3=[5;2];
N=50;

% Generate observations
error_var = 1;
e = sqrt(error_var)*randn(N,1);
G2 = ones(N,1);
G3 = [ones(N,1) [1:N]'];
y1 = e;
y2 = G2*theta2 + e;
y3 = G3*theta3 + e;

% ML estimate for theta
theta_ML2 = inv(transpose(G2)*G2)*transpose(G2)*y2;
theta_ML3 = inv(transpose(G3)*G3)*transpose(G3)*y3;

% Prior distribution parameters
theta_prior2 = 0;
theta_prior3 = [0;0];
theta_var = 10^(15);
prior_var2 = theta_var;
prior_var3 = theta_var*eye(2);

% Model marginal log likelihoods (prevent underflow issues)
model_11 = logmultigaussian(N,y1,0,error_var*eye(N));
model_21 = model_marginal_loglikelihood(N,1,G2,error_var,theta_prior2,prior_var2,y1);
model_31 = model_marginal_loglikelihood(N,2,G3,error_var,theta_prior3,prior_var3,y1);

model_12 = logmultigaussian(N,y2,0,error_var*eye(N));
model_22 = model_marginal_loglikelihood(N,1,G2,error_var,theta_prior2,prior_var2,y2);
model_32 = model_marginal_loglikelihood(N,2,G3,error_var,theta_prior3,prior_var3,y2);

model_13 = logmultigaussian(N,y3,0,error_var*eye(N));
model_23 = model_marginal_loglikelihood(N,1,G2,error_var,theta_prior2,prior_var2,y3);
model_33 = model_marginal_loglikelihood(N,2,G3,error_var,theta_prior3,prior_var3,y3);
BF23=model_23-model_33