N=100;
freq_norm=0.2;
n=1:N;
y=cos(freq_norm*pi*n); %essentially a rect window
y_hamm=y' .* hamming(N);
y_hann=y' .* hanning(N);
%fft implementation
% yfreq=20*log10(abs(fft(y)));
% yfreq_hamm=20*log10(abs(fft(y_hamm)));
% yfreq_hann=20*log10(abs(fft(y_hann)));
% plot(yfreq)
% hold on
% plot(yfreq_hamm)
% plot(yfreq_hann)
% hold off

%using built-in tool
%fvtool(y,1,y_hamm,2,y_hann,3)
% 
noise=0.1*randn(N,1)';
yn=y+noise; %essentially a rect window
yn_hamm=yn' .* hamming(N);
yn_hann=yn' .* hanning(N);
%fvtool(yn,2,y_hamm,3,yn_hamm,4,y_hann,5,yn_hann,6)
fvtool(y_hamm,1,yn_hamm,2)