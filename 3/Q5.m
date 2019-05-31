clear all;
rng(2,'twister');

N=100;
error_var=0.5G;
pole1=0.5+0.1j;
pole2=0.5-0.1j;
a1=pole1+pole2;
a2=-pole1*pole2;
% a1=0.1;
% a2=-0.2;
a3=0.3;

signal_AR=[0.1 0.2];
for i=2:N-1
    next=a1*signal_AR(i) + a2*signal_AR(i-1) + sqrt(error_var)*randn(1,1);
    signal_AR=[signal_AR next];
end

% signal_AR=[0.1 0.2 0.1];
% for i=3:N-1
%     next=a1*signal_AR(i) + a2*signal_AR(i-1) + a3*signal_AR(i-2) + sqrt(error_var)*randn(1,1);
%     signal_AR=[signal_AR next];
% end


% spec = abs(fft(signal_AR));
% figure(2)
% plot(spec)
% xlim([-N/2, 3*N/2])

% % P=2
% y=signal_AR(3:N)';
% G=fliplr(buffer(signal_AR(1:end-1), N-2, N-3, 'nodelay'));
% 
% % ML estimate for theta
% theta_ML = inv(transpose(G)*G)*transpose(G)*y
% 
% % Prior distribution parameters
% 
% % Setting prior as 0s and being very confident about it basically means 
% % that MAP tries to balance between the clearly non-zero data and the zero
% % prior estimation.
% theta_prior = [0;0];
% prior_var = eye(2);
% 
% % Likelihood parameters
% likelihood_var = error_var*inv((transpose(G)*G));
% 
% % Posterior parameters
% [theta_MAP,phi,big_theta,post_var] = post_param(G,error_var,theta_prior,prior_var,y);
% theta_MAP
% 
% y_ML = G*theta_ML;
% y_MAP = G*theta_MAP;
% 
% figure(1)
% plot(y)
% hold on
% plot(y_ML)
% plot(y_MAP)
% hold off
% title(["P=2, (ML MSE, MAP MSE)=" num2str(immse(y, y_ML)) ...
%     num2str(immse(y, y_MAP))])
% legend('Original', 'ML', 'MAP')
% 
% % P=10
% y=signal_AR(11:N)';
% G=fliplr(buffer(signal_AR(1:end-1), N-10, N-11, 'nodelay'));
% 
% % ML estimate for theta
% theta_ML = inv(transpose(G)*G)*transpose(G)*y
% 
% % Prior distribution parameters
% theta_prior = zeros(10,1);
% prior_var = eye(10);
% 
% % Likelihood parameters
% likelihood_var = error_var*inv((transpose(G)*G));
% 
% % Posterior parameters
% [theta_MAP,phi,big_theta,post_var] = post_param(G,error_var,theta_prior,prior_var,y);
% theta_MAP
% 
% y_ML = G*theta_ML;
% y_MAP = G*theta_MAP;
% 
% figure(2)
% plot(y)
% hold on
% plot(y_ML)
% plot(y_MAP)
% hold off
% title(["P=10, (ML MSE, MAP MSE)=" num2str(immse(y, y_ML)) ...
%     num2str(immse(y, y_MAP))])
% legend('Original', 'ML', 'MAP')

% % P=4
% y=signal_AR(5:N)';
% G=fliplr(buffer(signal_AR(1:end-1), 96, 95, 'nodelay'));
% 
% % ML estimate for theta
% theta_ML = inv(transpose(G)*G)*transpose(G)*y
% 
% % Prior distribution parameters
% theta_prior = zeros(4,1);
% prior_var = eye(4);
% 
% % Likelihood parameters
% likelihood_var = error_var*inv((transpose(G)*G));
% 
% % Posterior parameters
% [theta_MAP,phi,big_theta,post_var] = post_param(G,error_var,theta_prior,prior_var,y);
% theta_MAP
% 
% y_real = signal_AR(4:N-1)';
% y_ML = G*theta_ML;
% y_MAP = G*theta_MAP;
% 
% figure(2)
% plot(y_real)
% hold on
% plot(y_ML)
% plot(y_MAP)
% hold off
% title(["P=4, (ML MSE, MAP MSE)=" num2str(immse(y_real, y_ML)) ...
%     num2str(immse(y_real, y_MAP))])
% legend('Original', 'ML', 'MAP')

% % If using proper priors
% model_LL_list=[model_marginal_loglikelihood(N,1,signal_AR(1:99)',error_var,1,1,signal_AR(2:end)')];

model_LL_list=[];
model_MLMSE_list=[];
model_MAPMSE_list=[];
for i = 1:10
    G = fliplr(buffer(signal_AR(1:end-1), N-i, N-1-i, 'nodelay'));
    
    y=signal_AR(i+1:end)';
    
    % ML estimate for theta
    theta_ML = inv(transpose(G)*G)*transpose(G)*y;

    % Prior distribution parameters
    theta_prior = zeros(i,1);
    prior_var = eye(i);

    % Posterior parameters
    [theta_MAP,phi,big_theta,post_var] = post_param(G,error_var,theta_prior,prior_var,y);
    model_LL = model_marginal_loglikelihood(N,i,G,error_var,zeros(i,1),eye(i),signal_AR(i+1:end)');
%     model_LL = model_marginal_loglikelihood(N,i,G,error_var,[1; -0.26; zeros(i-2,1)],eye(i),signal_AR(i+1:end)');
    model_LL_list = [model_LL_list; model_LL];
    y_ML = G*theta_ML;
    y_MAP = G*theta_MAP;
    model_MLMSE_list = [model_MLMSE_list immse(y, y_ML)];
    model_MAPMSE_list = [model_MAPMSE_list immse(y, y_MAP)];
end

[best_value, best_hyperparam] = max(model_LL_list);
probs = exp(model_LL_list);
probs_norm = probs/sum(probs);
figure(1)
plot(1:10,probs_norm)
figure(2)
plot(1:10,model_MLMSE_list)
hold on 
plot(1:10,model_MAPMSE_list)
hold off
legend('ML MSEs', 'MAP MSEs')