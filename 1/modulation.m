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
a_per1=1+0.2*sin(0.01*pi*n);
a_per2=1+0.2*sin(0.05*pi*n);
persig1=a_per1.*signal;
persig2=a_per2.*signal;
% plot(persig)

%random AR1 modulation
a_ran=[1];
for i=1:N-1
    a_ran=[a_ran; a_ran(i)+0.05*randn(1,1)];
end
ransig=a_ran'.*signal;
%plot(ransig)

fvtool(persig2)