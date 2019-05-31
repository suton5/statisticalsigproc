N=100;
error_var=0.5;
pole1=0.5+0.1j;
pole2=0.5-0.1j;
a1=pole1+pole2;
a2=-pole1*pole2;

signal_AR=[0.01 0.01];
for i=2:N-1
    next=a1*signal_AR(i) + a2*signal_AR(i-1) + sqrt(error_var)*randn(1,1);
    signal_AR=[signal_AR next];
end

% spec = abs(fft(signal_AR));
% figure(2)
% plot(spec)
% xlim([-N/2, 3*N/2])

% P=2
g1=signal_AR(2:N-1)';
g2=signal_AR(1:N-2)';
y=signal_AR(3:N)';
G=[g1 g2];

% ML estimate for theta
theta_ML = inv(transpose(G)*G)*transpose(G)*y

% Prior distribution parameters

% Setting prior as 0s and being very confident about it basically means 
% that MAP tries to balance between the clearly non-zero data and the zero
% prior estimation.
theta_prior = [0;0];
prior_var = eye(2);

% Likelihood parameters
likelihood_var = error_var*inv((transpose(G)*G));

% Posterior parameters
[theta_MAP,phi,big_theta,post_var] = post_param(G,error_var,theta_prior,prior_var,y);
theta_MAP

y_real = G*[a1;a2];
y_ML = G*theta_ML;
y_MAP = G*theta_MAP;

figure(1)
plot(y_real)
hold on
plot(y_ML)
plot(y_MAP)
hold off
title(["P=2, (ML MSE, MAP MSE)=" num2str(immse(y_real, y_ML)) ...
    num2str(immse(y_real, y_MAP))])
legend('Original', 'ML', 'MAP')

% P=10
g1=signal_AR(10:N-1)';
g2=signal_AR(9:N-2)';
g3=signal_AR(8:N-3)';
g4=signal_AR(7:N-4)';
g5=signal_AR(6:N-5)';
g6=signal_AR(5:N-6)';
g7=signal_AR(4:N-7)';
g8=signal_AR(3:N-8)';
g9=signal_AR(2:N-9)';
g10=signal_AR(1:N-10)';

y=signal_AR(11:N)';
G=[g1 g2 g3 g4 g5 g6 g7 g8 g9 g10];

% ML estimate for theta
theta_ML = inv(transpose(G)*G)*transpose(G)*y

% Prior distribution parameters
theta_prior = zeros(10,1);
prior_var = eye(10);

% Likelihood parameters
likelihood_var = error_var*inv((transpose(G)*G));

% Posterior parameters
[theta_MAP,phi,big_theta,post_var] = post_param(G,error_var,theta_prior,prior_var,y);
theta_MAP

y_real = G*[a1;a2;zeros(8,1)];
y_ML = G*theta_ML;
y_MAP = G*theta_MAP;

figure(2)
plot(y_real)
hold on
plot(y_ML)
plot(y_MAP)
hold off
title(["P=10, (ML MSE, MAP MSE)=" num2str(immse(y_real, y_ML)) ...
    num2str(immse(y_real, y_MAP))])
legend('Original', 'ML', 'MAP')