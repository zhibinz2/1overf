% https://www.mathworks.com/help/signal/ref/filtfilt.html

wform = ecg(500);

plot(wform)
axis([0 500 -1.25 1.25])
text(155,-0.4,"Q")
text(180,1.1,"R")
text(205,-1,"S")

%%
% rng default
rng(0)

x = wform' + 0.25*randn(500,1);
d = designfilt("lowpassfir", ...
    'PassbandFrequency',0.15,'StopbandFrequency', 0.2, ...
    'PassbandRipple',1,'StopbandAttenuation',60, ...
    'DesignMethod',"equiripple");
y = filtfilt(d,x);
y1 = filter(d,x);

subplot(2,1,1)
plot([y y1])
title("Filtered Waveforms")
legend("Zero-phase Filtering","Conventional Filtering")

subplot(2,1,2)
plot(wform)
title("Original Waveform")

%%
d1 = designfilt("lowpassiir",'FilterOrder',12, ...
    'HalfPowerFrequency',0.15,'DesignMethod',"butter");
y = filtfilt(d1,x);

subplot(1,1,1)
plot(x)
hold on
plot(y,'LineWidth',3)
legend("Noisy ECG","Zero-Phase Filtering")