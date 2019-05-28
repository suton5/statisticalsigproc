clear all;
rng(1, 'twister');

% parameters: N=1024, threshold=1000, factor=1
% clean=0.5*cos([1:10000]*pi/4)+sin([1:10000]*pi/100);

% parameters: N=2048, threshold=0.8(best)/1.5(lowest MSE), factor=0.01
% [x, Fs]=audioread('../audio/piano_clean.wav');
% clean=x(:,1)';

[x, Fs]=audioread('../audio/grosse_original.wav');
clean=x(:,1)';

% parameters: N=2048, threshold=1.8(best)/40(lowest MSE), factor=0.01 
% [x, Fs]=audioread('../audio/f1lcapae.wav');
% clean=x(:,1)';

%data length
L=size(clean,2);
%add noise
factor=0.01;
noise=factor*randn(1,L);
x_in=clean+noise;

%estimate psd of noise
threshold=1.5;
spec=abs(fft(x_in));
spec(spec>threshold) = [];
Sn = mean(spec.^2);


%generate zero array of same length as x_in
y_out=0*x_in;
y_noisy=0*x_in;

%select length of each data partition
N=2048;
%select overlap
overlap=N/2;

%create a matrix of the data, with length-N partitions in each column
%with an overlap between successive partitions.
x=buffer(x_in,N,overlap);

%[length of each partition, number of partitions (zero padded if needed)]
[N_samps,N_frames]=size(x);

%create N_frames columns of repeated length-N rect windows as each
%column. perform windowing with the partitioned data matrix.
x_w=repmat(rectwin(N),1,N_frames).*x;

%for each data partition column, excluding the final 2
for frame_no=1:N_frames-2
    
    %perform FFT on each data partition
    X_w(:,frame_no)=fft(x_w(:,frame_no));
    
    %copy spectra into another matrix
    Y_w(:,frame_no)=X_w(:,frame_no);
    
%     %%%Filter implementations
%     
%     %%Bandpass filter
%     %attenuate by 0.1 for 2:64
%     Y_w(2:N/8,frame_no)=0.1*X_w(2:N/8,frame_no);
%     %attenuate by 0.2 for 129:256
%     Y_w(N/4+1:N/2,frame_no)=0.2*X_w(N/4+1:N/2,frame_no);
%     %fill up 258:512 with conjugate of 2:256, in reverse
%     Y_w(N:-1:N/2+2,frame_no)=conj(Y_w(2:N/2,frame_no));
    
    %%Wiener filter
    a = Y_w(:,frame_no);
    a(abs(a).^2 < Sn) = 0;
    a(abs(a).^2 > Sn) = (((abs(a(abs(a).^2 > Sn)).^2) - Sn) ./ ...
        (abs(a(abs(a).^2 > Sn)).^3)) .* a(abs(a).^2 > Sn);
    Y_w(:,frame_no) = a;

%     %%Wiener filter (modified exponent)
%     a = Y_w(:,frame_no);
%     a(abs(a).^2 < Sn) = 0;
%     a(abs(a).^2 > Sn) = (((abs(a(abs(a).^2 > Sn)).^2) - Sn) ./ ...
%         (abs(a(abs(a).^2 > Sn)).^1.5)) .* a(abs(a).^2 > Sn);
%     Y_w(:,frame_no) = a;

%     %%Spectral subtraction filter
%     a = Y_w(:,frame_no);
%     a(abs(a).^2 < Sn) = 0;
%     a(abs(a).^2 > Sn) = (1 - sqrt(Sn) ./ (abs(a(abs(a).^2 > Sn)))) ...
%         .* a(abs(a).^2 > Sn);
%     Y_w(:,frame_no) = a;

%     %%Power subtraction filter
%     a = Y_w(:,frame_no);
%     a(abs(a).^2 < Sn) = 0;
%     a(abs(a).^2 > Sn) = sqrt(1 - Sn ./ (abs(a(abs(a).^2 > Sn)).^2)) ...
%         .* a(abs(a).^2 > Sn);
%     Y_w(:,frame_no) = a;
    
    %transform processed spectra into time-domain
    y_w(:,frame_no)=ifft(Y_w(:,frame_no));
    
    %reconstruct proper signal, removing the added overlap
    y_out((frame_no-1)*overlap+1:(frame_no-1)*overlap+N)=...
     y_out((frame_no-1)*overlap+1:(frame_no-1)*overlap+N)+y_w(:,frame_no)';
    y_noisy((frame_no-1)*overlap+1:(frame_no-1)*overlap+N)=...
     y_noisy((frame_no-1)*overlap+1:(frame_no-1)*overlap+N)+x_w(:,frame_no)';
    
end    

% Normalise the outputs
clean=(clean./max(clean));
y_noisy=(y_noisy./max(y_noisy));
y_out=(y_out./max(y_out));

% %observed that shifting signal results in almost perfect overlap
% %(for sum of sinusoids signal)
% clean=clean(1000:9000);
% x_in=x_in(1000:9000);
% y_out=y_out(1112:9112);


initialerror=immse(clean, y_noisy);
finalerror=immse(clean, y_out);

% figure(1)
% plot(clean)
% ylim([-1.5 1.5])
% title("Clean")
% xlabel("Time")
% ylabel("Amplitude")
% 
% figure(2)
% plot(y_noisy)
% ylim([-1.5 1.5])
% title(["Noisy, MSE=" num2str(initialerror)])
% xlabel("Time")
% ylabel("Amplitude")
% 
% figure(3)
% plot(y_out)
% ylim([-1.5 1.5])
% title(["Processed, MSE=" num2str(finalerror)])
% xlabel("Time")
% ylabel("Amplitude")
% 
% figure(4)
% plot(abs(fft(clean)))
% title("Clean")
% xlabel("Frequency")
% ylabel("Amplitude")
% 
% figure(5)
% plot(abs(fft(y_noisy)))
% title(["Noisy, MSE=" num2str(initialerror)])
% xlabel("Frequency")
% ylabel("Amplitude")
% 
% figure(6)
% plot(abs(fft(y_out)))
% title(["Processed, MSE=" num2str(finalerror)])
% xlabel("Frequency")
% ylabel("Amplitude")
% sound(y_out, Fs)