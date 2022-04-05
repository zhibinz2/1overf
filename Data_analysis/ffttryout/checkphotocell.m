%% Bottom two analog singnal
% plot photocell analog signal
analog1=133; 
analog1data=dataEEG(133,:);
plot(datatimes,analog1data);


analog2=134;%trigger=139;
analog2data=dataEEG(134,:);
plot(datatimes,analog2data);

% plot fft series
F133=fft(analog1data);
A133=abs(F133);
N=length(analog1data);
fV=linspace(0,Fs/2,N/2+1);
plot(fV,A133(1:length(fV)));
xlim([0 200]);
xlabel('frequency');
ylabel('amplitude (uV)');
title('Spetra');

F134=fft(analog2data);
A134=abs(F134);
N=length(analog1data);
fV=linspace(0,Fs/2,N/2+1);
plot(fV,A134(1:length(fV)));
xlim([0 200]);
xlabel('frequency');
ylabel('amplitude (uV)');
title('Spetra');