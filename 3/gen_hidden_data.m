clear all

N=100;
sigma_n=2;
sigma_theta=1;

signal = [-3 5 -2 4 1 3 5 -1 2 4 6 5 -2 -2 1];

noise = sigma_n*randn(N,1);
theta = sigma_theta*randn(1);
offset = round(90*rand(1));

y=noise;

signal_offset=0*noise;
signal_offset(offset:offset+14)=signal'*theta;

y=y+signal_offset;
plot(y)
hold on
plot(signal_offset)
hold off