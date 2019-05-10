function X = DFT(x)
%DFT.m Compute DFT of data sequence

%Calculate length of sequence
length=size(x,1);

%Construct frequency matrix (vectorised)
temp=0:length-1;
freq=exp(-j*2*pi*temp' * temp/length);

%Construct frequency matrix (for loops)
% freq=zeros(length);
% for a=1:length
%     for b=1:length
%         freq(a,b)=exp(-j*2*pi*(a-1)*(b-1)/length);
%     end
% end

%Vectorised DFT
X = freq*x;

end