function [freqs,fcoef] = oneoverf(y,Fs);
% this function compute do fft, plot power spectrum
%   Dand fit the 1overf 
    fcoef=fft(y);
    N=length(y);
    fcoef=fcoef/N;
    halfN=floor(N/2);
    Fs=1/2;% assuming a sampling frequency of 0.5Hz (2 second per sample)
    freqs=linspace(0,Fs/2,halfN+1); % same df
    fcoef=2*fcoef(1:halfN);
    amplitude = abs(fcoef);
    p=amplitude.^2;
    
    %xx=log10(freqs(2:end));yy=log10(p(2:length(freqs)));
    xx=log10(freqs(2:end));yy=log10(p(1:length(p)));
    plot(xx,yy);hold on;
    tbl=table(xx',yy');
    mdl=fitlm(tbl,'linear');
    plotAdded(mdl);
    % xlabel('Log(f)');ylabel('Log(power)');title('C. Spectrum');
    %ylim([-3 3]);
    %ylim([-7 0]);
end

