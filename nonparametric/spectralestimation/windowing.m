N=100;
freq_norm=0.5;
n=1:N;

%essentially a rect window
y=cos(freq_norm*pi*n);
%implement other windows
y_hamm=y' .* hamming(N);
y_hann=y' .* hanning(N);
y_cheb=y' .* chebwin(N);

%using built-in tool
fvtool(y,1,y_hamm,2,y_hann,3,y_cheb,4)

%add different noise levels
noise1=0.01*randn(N,1)';
noise2=0.1*randn(N,1)';
yn1=y+noise1;
yn1_hamm=yn1' .* hamming(N);
yn2=y+noise2;
yn2_hamm=yn2' .* hamming(N);

%using built-in tool
fvtool(y_hamm,1,yn1_hamm,2,yn2_hamm,3)

%fft implementation
% yfreq=20*log10(abs(fft(y)));
% yfreq_hamm=20*log10(abs(fft(y_hamm)));
% yfreq_hann=20*log10(abs(fft(y_hann)));
% plot(yfreq)
% hold on
% plot(yfreq_hamm)
% plot(yfreq_hann)
% hold off