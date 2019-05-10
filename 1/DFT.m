function X = DFT(x)
%DFT.m Compute DFT of data sequence

%Calculate length of sequence
length=size(x,1);

%Construct frequency matrix (vectorised)
temp=0:length-1;
freq=exp(-j*2*pi*temp' * temp/length);

%Vectorised DFT
X = freq*x;
end