% Full 5th note in piano extract (AR 33)
[x, Fs]=audioread('../audio/piano_clean.wav');
xt=x(11800:12900,1)';
error_var=4.8e-07;

% % Full 6th note in piano extract (AR 42)
% x=audioread('../audio/piano_clean.wav');
% xt=x(14800:end,1)';
% error_var=2.45e-07;

% % The consonant 'J' in the word 'John' (AR 50)
% [x, Fs]=audioread('../audio/f1lcapae.wav');
% xt=x(4500:6001,1)';
% error_var=1.49e-04;

% % The vowel 'I' in the word 'resigned' (AR 50)
% [x, Fs]=audioread('../audio/f1lcapae.wav');
% xt=x(11000:13501,1)';
% error_var=8.21e-05;

% % Full organ extract (AR 49)
% [x, Fs]=audioread('../audio/organ.wav');
% xt=x(:,1)';
% error_var=2.35e-06;

% % Truncated organ extract (AR 49)
% [x, Fs]=audioread('../audio/organ.wav');
% xt=x(18000:23000,1)';
% error_var=2e-06;

% N=size(xt,2);
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

N=size(xt,2);
% Chosen AR order
P=33;
G = fliplr(buffer(xt(1:end-1), N-P, N-1-P, 'nodelay'));
y=xt(P+1:end)';
% ML estimate for theta
theta_ML = inv(transpose(G)*G)*transpose(G)*y;
% Prior distribution parameters
theta_prior = zeros(P,1);
prior_var = eye(P);
% Posterior parameters (inspecting values, closer to 0 due to the prior)
[theta_MAP,phi,big_theta,post_var] = post_param(G,error_var,theta_prior,prior_var,y);

% Delete L packets starting from p
p=300;
L=100;
xt_lossy = xt;
xt_lossy(p:p+L-1)=0;

x_f = xt_lossy;
x_b = xt_lossy;
x_predict1 = xt_lossy;
x_predict2 = xt_lossy;

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
    x_predict1(p-1+i) = alpha*x_f(p-1+i) + (1-alpha)*x_b(p-1+i);
end

% Bayesian interpolation
temp = [-flip(theta_MAP') 1];
A = convmtx(temp, N-P);
A_i = A(:,p:p+L-1);
A_neg_i = [A(:,1:p-1), A(:,p+L:end)];
x_neg_i = [xt_lossy(1:p-1), xt_lossy(p+L:end)];
x_LS = -inv(A_i'*A_i)*A_i'*A_neg_i*x_neg_i';
x_predict2 = [xt_lossy(1:p-1), x_LS', xt_lossy(p+L:end)];

immse(xt, xt_lossy)
immse(xt, x_f)
immse(xt, x_b)
immse(xt, x_predict1)
immse(xt, x_predict2)