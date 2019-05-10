% Plot sound files
%x=audioread('audio/alphabet/A.wav');
%xt=x(1200:3500,1);
%x=audioread('audio/piano_clean.wav');
%xt=x(12100:12300,1);
[x, Fs]=audioread('audio/piano_clean.wav');
xt=x(12100:12900,1);
N=size(xt,1);
% xhamm=xt.*hamming(N);
% xhann=xt.*hanning(N);
% xcheb=xt.*chebwin(N);
%fvtool(xt,1,xhamm,2,xhann,3,xcheb,4)

% Plot spectrum
%audiofft(xt, 'rect', N);
%hold on
xw=xt.* hamming(N);
Y = fft(xw);
P2 = abs(Y/N);
P1 = P2(1:N/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = Fs*(0:(N/2))/N;
plot(f,P1) 
title('Single-Sided Frequency Spectrum')
xlabel('f (Hz)')
ylabel('Amplitude')

% audiofft(xt, 'hamm', N);
% %audiofft(xt, 'hamm', N);
%audiofft(xt, 'cheb', N);
%legend('rect','hamm','hann','cheb')
%hold off
