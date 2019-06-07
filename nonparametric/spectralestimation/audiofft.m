function [f,P1] = audiofft(data, window, window_length)
% Generates one-sided frequency spectrum of a .wav file, with optional
% windows
if window == 'rect'
    data=data.* rectwin(window_length);
elseif window == 'hamm'
    data=data.* hamming(window_length);
elseif window == 'hann'
    data=data.* hanning(window_length);
elseif window == 'cheb'
    data=data.* chebwin(window_length);
end
        
Y = fft(data);
L = length(Y);
P2 = abs(Y(:,1)/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);

plot(P1)
grid
title('Single-Sided Frequency Spectrum')
ylabel('Amplitude')
end

