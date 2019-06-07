clear all

% % Full 5th note in piano extract (AR 2)
% [x, Fs]=audioread('../audio/piano_clean.wav');
% xt=x(11800:12900,1)';
% error_var=1.6105e-05;

% % The consonant 'J' in the word 'John' (AR 2)
% [x, Fs]=audioread('../audio/f1lcapae.wav');
% xt=x(5380:6001,1)';
% error_var=9.7499e-05;

% % The vowel 'I' in the word 'resigned' (AR 3)
% [x, Fs]=audioread('../audio/f1lcapae.wav');
% xt=x(11230:13501,1)';
% error_var=1.0315e-04;

% % The vowel 'E' (AR 7)
% [x, Fs]=audioread('../audio/alphabet/E.wav');
% xt=x(1000:4300,1)';
% error_var=3.2507e-04;

% % Missing file 1 (AR 16)
% [x, Fs]=audioread('../audio/missing/armst_37_missing.wav');
% xt=x(:,1)';
% error_var=9.4588e-04;

% Missing file 2 (AR 2)
[x, Fs]=audioread('../audio/missing/grosse_40_percent_missing.wav');
xt=x(:,1)';
error_var=3.7607e-06;

N=size(xt,2);


% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

% model_LL_list=[];
% % model_MLMSE_list=[];
% % model_MAPMSE_list=[];
% for i = 1:20
%     G = fliplr(buffer(xt(1:end-1), N-i, N-1-i, 'nodelay'));
%     
%     y=xt(i+1:end)';
%     theta_ML = inv(transpose(G)*G)*transpose(G)*y;
%     % Prior distribution parameters
%     theta_prior = zeros(i,1);
%     prior_var = eye(i);
%     % Posterior parameters
%     [theta_MAP,phi,big_theta,post_var] = post_param(G,error_var,theta_prior,prior_var,y);
%     model_LL = model_marginal_loglikelihood(N,i,G,error_var,zeros(i,1),eye(i),xt(i+1:end)');
%     model_LL_list = [model_LL_list; model_LL];
% end
% 
% [best_value, best_hyperparam] = max(model_LL_list);
% probs_norm = exp(model_LL_list - logsumexp(model_LL_list));
% figure(1)
% plot(probs_norm)

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 


P=16;
% Chosen AR order
G = fliplr(buffer(xt(1:end-1), N-P, N-1-P, 'nodelay'));
y=xt(P+1:end)';
% ML estimate for theta
theta_ML = inv(transpose(G)*G)*transpose(G)*y;
% Prior distribution parameters
theta_prior = zeros(P,1);
prior_var = eye(P);
% Posterior parameters (inspecting values, closer to 0 due to the prior)
[theta_MAP,phi,big_theta,post_var] = post_param(G,error_var,...
    theta_prior,prior_var,y);

y_ML = G*theta_ML;
y_MAP = G*theta_MAP;

% figure(2)
% plot(y)
% hold on
% plot(y_ML)
% plot(y_MAP)
% hold off
% legend('Real','ML','MAP')

% specML = abs(fft(y_ML)).^2;
% specMAP = abs(fft(y_MAP)).^2;
% figure(1)
% plot(specML)
% hold on
% plot(specMAP)
% hold off
% xlabel('Frequency')
% ylabel('|DFT|^2')
% title('Squared magnitude of the frequency spectrum')
% legend('ML','MAP')

theta_ML;
theta_MAP;
var(y-y_MAP)