% Full 5th note in piano extract
% x=audioread('audio/piano_clean.wav');
% xt=x(11800:12900,1);
% plot(xt)

% Full 6th note in piano extract
% x=audioread('audio/piano_clean.wav');
% xt=x(14800:end,1);
% plot(xt)

% The consonant 'J' in the word 'John'
% [x, Fs]=audioread('audio/f1lcapae.wav');
% xt=x(4500:6001,1);
% plot(xt)

% The vowel 'I' in the word 'resigned'
% [x, Fs]=audioread('audio/f1lcapae.wav');
% xt=x(11000:13501,1);
% plot(xt) 

% Full organ extract
% [x, Fs]=audioread('audio/organ.wav');
% xt=x(:,1);
% plot(xt)

% Truncated organ extract
[x, Fs]=audioread('audio/organ.wav');
xt=x(18000:23000,1);
% plot(xt)

% Compute length of truncated sequence
N=size(xt,1);

% Plot spectrum using Hamming window
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