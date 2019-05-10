tdft=[];
tfft=[];
N=floor(logspace(1.4,3,500))

for i = N
    tempdft=[];
    tempfft=[];
    for j=1:1
        x = randn(i,1);

        tic;
        X1 = DFT(x);
        toc;
        t1 = toc;

        tic;
        X2 = fft(x);
        toc;
        t2 = toc;

        tempdft = [tempdft;t1];
        tempfft = [tempfft;t2];
    end
    tdft = [tdft; mean(tempdft)];
    tfft = [tfft; mean(tempfft)];
end
n2=0.0000001*N.^2;
nlogn=10^(-7)*N.*log2(N);
figure
loglog(N,tdft)
hold on
loglog(N,tfft)
loglog(N,n2)
loglog(N,nlogn)
hold off