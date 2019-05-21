clear all;
rng(1,'twister');

%create clean signal
%N=1024, wiener, thres=1000
% clean=0.5*cos([1:10000]*pi/4)+sin([1:10000]*pi/100);

[x, Fs]=audioread('../audio/piano_clean.wav');
clean=x(:,1)';

% [x, Fs]=audioread('../audio/armst_37_orig.wav');
% clean=x(:,1)';

%data length
L=size(clean,2);
%add noise
factor=0.01;
x_in=clean+factor*randn(1,L);

%estimate psd of noise, method 1
threshold=20;
spec=abs(fft(x_in));
spec(spec>threshold) = [];
noisepsd = mean(spec.^2);

% %estimate psd of noise, method 2
% noisepsd = factor * L;

%generate zero array of same length as x_in
y_out=0*x_in;

%select length of each data partition
N=1024;
%select overlap
overlap=N/2;

%create a matrix of the data, with length-N partitions in each column
%with an overlap between successive partitions.
x=buffer(x_in,N,overlap);

%[length of each partition, number of partitions (zero padded if needed)]
[N_samps,N_frames]=size(x);

%create N_frames columns of repeated length-N hanning windows as each
%column. perform windowing with the partitioned data matrix.
x_w=repmat(rectwin(N),1,N_frames).*x;

%for each data partition column, excluding the final 2
for frame_no=1:N_frames-2
    
    %perform FFT on each data partition
    X_w(:,frame_no)=fft(x_w(:,frame_no));
    
    %copy spectra into another matrix
    Y_w(:,frame_no)=X_w(:,frame_no);
    
    %%%Filter implementations
    
%     %%Bandpass filter
%     %attenuate by 0.1 for 2:64
%     Y_w(2:N/8,frame_no)=0.1*X_w(2:N/8,frame_no);
%     %attenuate by 0.2 for 129:256
%     Y_w(N/4+1:N/2,frame_no)=0.2*X_w(N/4+1:N/2,frame_no);
%     %fill up 258:512 with conjugate of 2:256, in reverse
%     Y_w(N:-1:N/2+2,frame_no)=conj(Y_w(2:N/2,frame_no));
    
    %%Wiener filter
    a = Y_w(:,frame_no);
    a(abs(a).^2 < noisepsd) = 0;
    a(abs(a).^2 > noisepsd) = (1 - noisepsd ./ (abs(a(abs(a).^2 > noisepsd)).^2)) .* a(abs(a).^2 > noisepsd);
    Y_w(:,frame_no) = a;

%     %%Spectral subtraction filter
%     a = Y_w(:,frame_no);
%     a(abs(a).^2 < noisepsd) = 0;
%     a(abs(a).^2 > noisepsd) = (1 - sqrt(noisepsd) ./ (abs(a(abs(a).^2 > noisepsd)))) .* a(abs(a).^2 > noisepsd);
%     Y_w(:,frame_no) = a;
    
%     %%Power subtraction filter
%     a = Y_w(:,frame_no);
%     a(abs(a).^2 < noisepsd) = 0;
%     a(abs(a).^2 > noisepsd) = sqrt(1 - noisepsd ./ (abs(a(abs(a).^2 > noisepsd)).^2)) .* a(abs(a).^2 > noisepsd);
%     Y_w(:,frame_no) = a;
    
    %transform processed spectra into time-domain
    y_w(:,frame_no)=ifft(Y_w(:,frame_no));
    
    %reconstruct proper signal, removing the added overlap
    y_out((frame_no-1)*overlap+1:(frame_no-1)*overlap+N)=y_out((frame_no-1)*overlap+1:(frame_no-1)*overlap+N)+y_w(:,frame_no)';
    
end    

clean=(clean./max(clean));
x_in=(x_in./max(x_in));
y_out=(y_out./max(y_out));

% %observed that shifting signal results in almost perfect overlap
% clean=clean(1000:9000);
% x_in=x_in(1000:9000);
% y_out=y_out(1112:9112);


initialerror=immse(clean, x_in);
finalerror=immse(clean, y_out);
%percentage reduction
per=(initialerror-finalerror)/initialerror*100;

figure(1)
plot(clean)
ylim([-1.5 1.5])
title('Clean')

figure(2)
plot(x_in)
ylim([-1.5 1.5])
title(['Noisy, MSE=' num2str(initialerror)])

figure(3)
plot(y_out)
ylim([-1.5 1.5])
title(['Processed, MSE=' num2str(finalerror)])% 'Percentage red=' num2str(per)])

figure(4)
plot(abs(fft(clean)))
title('Clean')

figure(5)
plot(abs(fft(x_in)))
title(['Noisy, MSE=' num2str(initialerror)])

figure(6)
plot(abs(fft(y_out)))
title(['Processed, MSE=' num2str(finalerror)])