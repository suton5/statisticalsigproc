%create noisy signal
x_in=0.5*cos([1:10000]*pi/4)+sin([1:10000]*pi/100)+randn(1,10000);

%generate zero array of same length as x_in
y_out=0*x_in;

N=512;
%use overlap M=N/2
overlap=256;

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
    
    %process 2:64, 129:256, 258:512 of the spectra
    
    %attenuate by 0.1 for 2:64
    Y_w(2:N/8,frame_no)=0.1*X_w(2:N/8,frame_no);
    %attenuate by 0.2 for 129:256
    Y_w(N/4+1:N/2,frame_no)=0.2*X_w(N/4+1:N/2,frame_no);
    %fill up 258:512 with conjugate of 2:256, in reverse
    Y_w(N:-1:N/2+2,frame_no)=conj(Y_w(2:N/2,frame_no));
    
    %transform processed spectra into time-domain
    y_w(:,frame_no)=ifft(Y_w(:,frame_no));
    
    %reconstruct original signal, removing the added overlap
    y_out((frame_no-1)*overlap+1:(frame_no-1)*overlap+N)=y_out((frame_no-1)*overlap+1:(frame_no-1)*overlap+N)+y_w(:,frame_no)';
    
end    
    


