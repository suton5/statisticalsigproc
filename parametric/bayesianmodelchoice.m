clear all;
rng(3,'twister');

% Select a value for theta and data size
theta2=5;
theta3=[5;2];
N=100;

% Generate observations
error_var = 1;
e = sqrt(error_var)*randn(N,1);
G2 = ones(N,1);
G3 = [ones(N,1) [1:N]'];
y1 = e;
y2 = G2*theta2 + e;
y3 = G3*theta3 + e;

% ML estimate for theta
theta_ML2 = inv(transpose(G2)*G2)*transpose(G2)*y3;
theta_ML3 = inv(transpose(G3)*G3)*transpose(G3)*y3;

% Prior distribution parameters
theta_prior2 = 0;
theta_prior3 = [0;0];
theta_var = 10^(0);
prior_var2 = theta_var;
prior_var3 = theta_var*eye(2);

[theta_MAP2,phi2,big_theta2,post_var2] = post_param(G2,error_var,...
    theta_prior2,prior_var2,y3);
[theta_MAP3,phi3,big_theta3,post_var3] = post_param(G3,error_var,...
    theta_prior3,prior_var3,y3);

% plot(y3)
% hold on
% plot(zeros(N,1))
% plot(G2*theta_MAP2)
% plot(G3*theta_MAP3)
% hold off
% xlabel('n')
% ylabel('x_n')
% ylim([-1, 120])
% legend('Real data','Model 1','Model 2','Model 3')
% title(sprintf('Data generated from Model 3, with N=%s and errorvar=%s', num2str(N), num2str(error_var)))


% bar(exp([model_13, model_23, model_33]))
% xlabel('Model')
% ylabel('Model likelihood')
% title(sprintf('Data generated from Model 3, with N=%s and errorvar=%s', num2str(N), num2str(error_var)))

list1=[];
list2=[];
list3=[];

for i=0:15
    theta_var=10^i;
    prior_var2 = theta_var;
    prior_var3 = theta_var*eye(2);

    [theta_MAP2,phi2,big_theta2,post_var2] = post_param(G2,error_var,...
        theta_prior2,prior_var2,y3);
    [theta_MAP3,phi3,big_theta3,post_var3] = post_param(G3,error_var,...
        theta_prior3,prior_var3,y3);
    model_1 = logmultigaussian(N,y3,0,error_var*eye(N));
    model_2 = model_marginal_loglikelihood(N,1,G2,error_var,theta_prior2,prior_var2,y3);
    model_3 = model_marginal_loglikelihood(N,2,G3,error_var,theta_prior3,prior_var3,y3);
    sum=exp(model_1)+exp(model_2)+exp(model_3);
    model_1=exp(model_1)/sum;
    model_2=exp(model_2)/sum;
    model_3=exp(model_3)/sum;
    list1=[list1 model_1];
    list2=[list2 model_2];
    list3=[list3 model_3];
end
figure(1)
plot(0:15, list1)
title(sprintf('Model likelihood for Model 1, data generated from Model 3, with N=%s and errorvar=%s', num2str(N), num2str(error_var)))
xlabel('logpriorvar')
ylabel('normalise likelihood')
figure(2)
plot(0:15, list2)
title(sprintf('Model likelihood for Model 2, data generated from Model 3, with N=%s and errorvar=%s', num2str(N), num2str(error_var)))
xlabel('logpriorvar')
ylabel('normalise likelihood')
figure(3)
plot(0:15, list3)
title(sprintf('Model likelihood for Model 3, data generated from Model 3, with N=%s and errorvar=%s', num2str(N), num2str(error_var)))
xlabel('logpriorvar')
ylabel('normalise likelihood')
    