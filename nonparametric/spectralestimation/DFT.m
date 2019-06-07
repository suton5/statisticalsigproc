function xfreq = DFT(xdata)
%DFT.m Computes DFT of data sequence

%Calculate length of sequence
length=size(xdata,1);

%Construct frequency matrix (vectorised)
temp=0:length-1;
freq=exp(-j*2*pi*temp' * temp/length);

%Vectorised DFT
xfreq = freq*xdata;
end