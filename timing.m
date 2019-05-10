tdft=[];
tfft=[];

for i=0:4

    x = randn(10^i,1);

    tic
    X1 = DFT(x);
    toc
    t1 = toc;

    tic
    X2 = fft(x);
    toc
    t2 = toc;
    
    tdft = [tdft;t1];
    tfft = [tfft;t2];
    
    
end

N = 0:4;
figure
semilogy(N,tdft)
hold on
semilogy(N,tfft)
hold off