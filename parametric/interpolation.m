clear all

% 5th note in piano extract (AR 2)
[x, Fs]=audioread('../audio/piano_clean.wav');
xt=x(11800:12900,1)';
error_var=1.6105e-05;

% % The consonant 'J' (AR 2)
% [x, Fs]=audioread('../audio/f1lcapae.wav');
% xt=x(5380:6001,1)';
% error_var=9.7499e-05;

% % The vowel 'E' (AR 7)
% [x, Fs]=audioread('../audio/alphabet/E.wav');
% xt=x(1000:4300,1)';
% error_var=3.2507e-04;

P=2;

N=size(xt,2);
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

% % Forward prediction mode (w adding noise)
% for packet = p:p+L-1
%     x_f(packet) = x_f(packet-1:-1:packet-P)*theta_MAP + ...
%         sqrt(error_var)*randn(1,1);
% end
% 
% % Backward prediction mode (w adding noise)
% for packet = p+L-1:-1:p
%     x_b(packet) = x_b(packet+1:packet+P)*theta_MAP + ...
%         sqrt(error_var)*randn(1,1);
% end

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

[immse(xt, x_predict1), immse(xt, x_predict2)]
(immse(xt, x_predict1) - immse(xt, x_predict2))/immse(xt, x_predict1)

plot(xt)
hold on
plot(x_predict1)
plot(x_predict2)
hold off
legend('Real','Weighted Interpolation','Bayesian Interpolation')
xlabel('n')
ylabel('x_n')
title('Interpolation of piano musical note')

% plot(xt_lossy)
% xlabel('n')
% ylabel('x_n')
% title('Piano musical note, with 100 packets zeroed out')