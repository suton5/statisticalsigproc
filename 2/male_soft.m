clear all;

[x, Fs]=audioread('../audio/tests/male_speech_noisy_soft.wav');
x_in=x(:,1)';

%estimate psd of noise
threshold=0.25;
spec=abs(fft(x_in));
spec(spec>threshold) = [];
Sn = mean(spec.^2);

%generate zero array of same length as x_in
y_out=0*x_in;

%select length of each data partition
N=2048;
%select overlap
overlap=N/2;

%create a matrix of the data, with length-N partitions in each column
%with an overlap between successive partitions.
x=buffer(x_in,N,overlap);

%[length of each partition, number of partitions (zero padded if needed)]
[N_samps,N_frames]=size(x);

%create N_frames columns of repeated length-N hanning windows as each
%column. perform windowing with the partitioned data matrix.
x_w=repmat(hanning(N),1,N_frames).*x;

%for each data partition column, excluding the final 2
for frame_no=1:N_frames-2
    
    %perform FFT on each data partition
    X_w(:,frame_no)=fft(x_w(:,frame_no));
    
    %copy spectra into another matrix
    Y_w(:,frame_no)=X_w(:,frame_no);
    
    %%Wiener filter
    a = Y_w(:,frame_no);
    a(abs(a).^2 < Sn) = 0;
    a(abs(a).^2 > Sn) = (((abs(a(abs(a).^2 > Sn)).^2) - Sn) ./ ...
        (abs(a(abs(a).^2 > Sn)).^2)) .* a(abs(a).^2 > Sn);
    Y_w(:,frame_no) = a;

    %transform processed spectra into time-domain
    y_w(:,frame_no)=ifft(Y_w(:,frame_no));
    
    %reconstruct proper signal, removing the added overlap
    y_out((frame_no-1)*overlap+1:(frame_no-1)*overlap+N)=...
     y_out((frame_no-1)*overlap+1:(frame_no-1)*overlap+N)+y_w(:,frame_no)';
 
end    

sound(y_out, Fs)