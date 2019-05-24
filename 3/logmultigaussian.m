function logmultivargaussian = logmultigaussian(N,y,mean,var)
%   Compute log probability of observed data under noise model
quadterm = transpose(y-mean)*inv(var)*(y-mean);
logmultivargaussian = -log(det(2*pi*var))/2 - quadterm/2;
end

