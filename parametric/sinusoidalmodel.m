clear all;
rng(3,'twister');

% Truncated organ extract
[x, Fs]=audioread('../audio/organ.wav');
xt=x(18000:23000,1);

J = 5;
N = size(xt, 1);
w_vect = [0.25*pi, 0.5*pi, 0.75*pi, pi, 1.25*pi];
G = zeros(N, 2*J);

% Populate with cos
for i=1:N
    for j=1:J
        G(i, 2*j-1) = cos(i*w_vect(j));
    end
end

% Populate with sin
for i=1:N
    for j=1:J
        G(i, 2*j) = sin(i*w_vect(j));
    end
end

% Generate data
error_var = 1;
% theta=-0.5:0.1:0.4;
theta=[0.1, 0.2];
y = G*theta' + sqrt(error_var)*randn(N,1);

% ML estimate for theta
theta_ML = (G'*G)\G'*xt;

% Prior distribution parameters
theta_prior = zeros(2*J, 1);
prior_var = eye(2*J);

% Posterior parameters
[theta_MAP,phi,big_theta,post_var] = post_param(G,error_var,...
    theta_prior,prior_var,xt);

% plot(y)
% hold on
% plot(G*theta_ML)
% plot(G*theta_MAP)
% hold off

% plot(abs(fft([xt; zeros(100000, 1)])))

fvtool(G*theta_MAP)