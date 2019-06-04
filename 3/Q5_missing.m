% % Missing file 1 (AR )
% [x, Fs]=audioread('../audio/missing/armst_37_missing.wav');
% xt=x(:,1)';
% error_var=4.8e-07;

% Missing file 2 (AR )
[x, Fs]=audioread('../audio/missing/grosse_40_percent_missing.wav');
xt=x(:,1)';
error_var=3.12e-06;

N=size(xt,2);

% model_LL_list=[];
% model_MLMSE_list=[];
% model_MAPMSE_list=[];
% for i = 1:50
%     G = fliplr(buffer(xt(1:end-1), N-i, N-1-i, 'nodelay'));
%     
%     y=xt(i+1:end)';
%     
%     % ML estimate for theta
%     theta_ML = inv(transpose(G)*G)*transpose(G)*y;
% 
%     % Prior distribution parameters
%     theta_prior = zeros(i,1);
%     prior_var = eye(i);
% 
%     % Posterior parameters
%     [theta_MAP,phi,big_theta,post_var] = post_param(G,error_var,theta_prior,prior_var,y);
%     model_LL = model_marginal_loglikelihood(N,i,G,error_var,zeros(i,1),eye(i),xt(i+1:end)');
%     model_LL_list = [model_LL_list; model_LL];
%     y_ML = G*theta_ML;
%     y_MAP = G*theta_MAP;
%     model_MLMSE_list = [model_MLMSE_list immse(y, y_ML)];
%     model_MAPMSE_list = [model_MAPMSE_list immse(y, y_MAP)];
% end
% 
% [best_value, best_hyperparam] = max(model_LL_list);
% probs_norm = exp(model_LL_list - logsumexp(model_LL_list));
% figure(1)
% plot(probs_norm)
% figure(2)
% plot(model_MLMSE_list)
% hold on 
% plot(model_MAPMSE_list)
% hold off
% legend('ML MSEs', 'MAP MSEs')

% Chosen AR order
P=50;
G = fliplr(buffer(xt(1:end-1), N-P, N-1-P, 'nodelay'));
y=xt(P+1:end)';
% ML estimate for theta
theta_ML = inv(transpose(G)*G)*transpose(G)*y;
% Prior distribution parameters
theta_prior = zeros(P,1);
prior_var = eye(P);
% Posterior parameters (inspecting values, closer to 0 due to the prior)
[theta_MAP,phi,big_theta,post_var] = post_param(G,error_var,theta_prior,prior_var,y);

x_f = xt;
x_b = xt;
x_predict = xt;

% Number of missing packets
L=25000;
% Starting point
p=41000;

% Forward prediction mode (w/o adding noise)
for packet = p:p+L-1
    x_f(packet) = x_f(packet-1:-1:packet-P)*theta_MAP;
end

% Backward prediction mode (w/o adding noise)
for packet = p+L-1:-1:p
    x_b(packet) = x_b(packet+1:packet+P)*theta_MAP;
end

% Weighted sum
for i = 1:L
    alpha = (L-i)/(L-1);
    x_predict(p-1+i) = alpha*x_f(p-1+i) + (1-alpha)*x_b(p-1+i);
end

% Bayesian interpolation
temp = [-flip(theta_MAP') 1];
A = convmtx(temp, N-P);
A_i = A(:,p:p+L-1);
A_neg_i = [A(:,1:p-1), A(:,p+L:end)];
x_neg_i = [xt(1:p-1), xt(p+L:end)];
x_LS = -inv(A_i'*A_i)*A_i'*A_neg_i*x_neg_i';
x_predict2 = [xt(1:p-1), x_LS', xt(p+L:end)];

plot(xt)
hold on
% plot(x_f)
plot(x_predict2)
% plot(x_predict)
hold off
% legend('lossy','for','back','weighted')