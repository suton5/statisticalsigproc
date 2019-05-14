N=1000;
n=1:N;
signal=real(exp(j*0.1*pi*n));

%linear modulation
A=0.05;
B=1;
a_lin=A*n+B;
linsig=a_lin.*signal;
% plot(linsig)

%periodic modulation
a_per=1+0.2*sin(0.01*pi*n);
persig=a_per.*signal;
% plot(persig)

%random AR1 modulation
a_ran=[1];
for i=1:N-1
    a_ran=[a_ran; a_ran(i)+0.05*randn(1,1)];
end
ransig=a_ran'.*signal;
plot(ransig)