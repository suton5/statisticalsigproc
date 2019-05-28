clear all

N=100;
P=1;
sigma_n=2;
sigma_theta=1;

signal = [-3 5 -2 4 1 3 5 -1 2 4 6 5 -2 -2 1];

% noise = sigma_n*randn(N,1);
% theta = sigma_theta*randn(1);
% offset = round(90*rand(1));
% y=noise;
% 
% signal_offset=0*noise;
% signal_offset(offset:offset+14)=signal'*theta;
% 
% y=y+signal_offset;
% plot(y)
% hold on
% plot(signal_offset)
% hold off

str=load('hidden_data.mat');
y=str.y;

model_LL_list=[];
for i = 1:90
    G = zeros(N,1);
    G(i:i+14) = signal';
    G = G(1:N,:);
    model_LL = model_marginal_loglikelihood(N,P,G,sigma_n^2,0,sigma_theta^2,y);
    model_LL_list = [model_LL_list; model_LL];
end

[best_value, best_offset] = max(model_LL_list);
probs = exp(model_LL_list);
probs_norm = probs/sum(probs);
figure(1)
plot(probs_norm)

G = zeros(N,1);
G(best_offset:best_offset+14) = signal';

theta_ML = inv(transpose(G)*G)*transpose(G)*y;

error_var = sigma_n^2;

% Prior distribution parameters
theta_prior = 0;
prior_var = sigma_theta^2;

% Likelihood parameters
likelihood_var = error_var/(transpose(G)*G);

% Posterior parameters
[theta_MAP,phi,big_theta,post_var] = post_param(G,error_var,theta_prior,prior_var,y);

% Posterior distribution
x = [-2:0.001:2];
y_post = normpdf(x,theta_MAP,sqrt(post_var));
y_likelihood = normpdf(x,theta_ML,sqrt(likelihood_var));
y_prior = normpdf(x,theta_prior, sqrt(prior_var));
figure(2)
plot(x, y_post) 
hold on
plot(x, y_likelihood)
plot(x, y_prior)
xline(theta_ML,'-.','ML');
xline(theta_MAP,'-.','MAP');
legend('Posterior','Likelihood','Prior')
hold off

% Null hypothesis
null_hypothesis = y_post(2001)