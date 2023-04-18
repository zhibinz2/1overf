function [fcoef1,fcoef2,cprod, freqs] = spectra(x1,x2,maxfreq,sr);
% x1 and x2 should be of the same length in the first dimension
% This function can be used for Powerspectra or Cross spectra

SamplesLength=size(x2,1);
T=SamplesLength/sr; % in second 
df=1/T;
nbins = ceil(maxfreq/df) + 1;
freqs = [0:(nbins-1)]*df;

fcoef1 = fft(ndetrend(x1,1),[],1)/SamplesLength;
fcoef2 = fft(ndetrend(x2,1),[],1)/SamplesLength;

cprod=[];
for i=1:nbins
    cprod(i,:)=fcoef1(i,:).*conj(fcoef2(i,:))/SamplesLength;
end

end

