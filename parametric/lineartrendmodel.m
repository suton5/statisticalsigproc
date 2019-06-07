% Select a value for theta and data size
theta=[5;2];
N=50;
error_var = 1;
% Prior distribution parameters
theta_prior = [4.5;1.8];
prior_var = 0.1*eye(2);

%Generate observations
g1 = ones(N,1);
g2 = [1:N]';
G = [g1 g2];
y = G*theta + sqrt(error_var)*randn(N,1);

% ML estimate
theta_ML = inv(transpose(G)*G)*transpose(G)*y

% Posterior parameters
[theta_MAP,phi,big_theta,post_var] = post_param(G,error_var,...
    theta_prior,prior_var,y);
theta_MAP

% Generate plots
% plot(y);
% text(2,y(end)-10,...
%     sprintf('Actual theta1: %4.2f \nActual theta2: %4.2f \nPrior theta1: %4.2f \nPrior theta2: %4.2f \nML theta1: %4.2f \nML theta2: %4.2f \nMAP theta1: %4.2f \nMAP theta2: %4.2f',...
%     theta(1), theta(2), theta_prior(1), theta_prior(2), theta_ML(1), theta_ML(2), theta_MAP(1), theta_MAP(2)), 'FontSize', 18);
% xlabel('n')
% ylabel('y_n')
% title(['Parameter estimation with errorvar=1, N=', num2str(N)])