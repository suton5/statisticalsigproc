N=100;
P=2;
error_var=2;
pole1=0.5+0.1j;
pole2=0.5-0.1j;
a1=pole1+pole2;
a2=-pole1*pole2;

signal_AR=[1 1];
for i=2:N-1
    next=a1*signal_AR(i) + a2*signal_AR(i-1) + sqrt(error_var)*randn(1,1);
    signal_AR=[signal_AR next];
end

figure(1)
plot(signal_AR)

spec = abs(fft(signal_AR));
figure(2)
plot(spec)
xlim([-N/2, 3*N/2])

g1=signal_AR(2:N-1)';
g2=signal_AR(1:N-2)';
y=signal_AR(3:N)';
G=[g1 g2];

% ML estimate for theta
theta_ML = inv(transpose(G)*G)*transpose(G)*y

% Prior distribution parameters
theta_prior = [1;-0.26];
prior_var = [0.1 0 ; 0 0.5];

% Likelihood parameters
likelihood_var = error_var*inv((transpose(G)*G));

% Posterior parameters
[theta_MAP,phi,big_theta,post_var] = post_param(G,error_var,theta_prior,prior_var,y);
theta_MAP