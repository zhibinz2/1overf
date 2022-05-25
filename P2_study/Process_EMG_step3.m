%% plot raw EMG
timeL;samplesL;TRIGGERindL;srL;channels_infoL;
timeR;samplesR;TRIGGERindR;srR;channels_infoR;

figure;
subplot(2,1,1);plot(samplesL(33,:)');
subplot(2,1,2);plot(samplesR(33,:)');

%% Extract EMG
EMGL=samplesL(33,:)';
EMGR=samplesR(33,:)';

%% EMG preprocessing 
% detrend
detrend_EMGL=detrend(EMGL,2); % even better
detrend_EMGR=detrend(EMGR,2); 
figure('units','normalized','outerposition',[0 0 1 0.6]);
subplot(2,1,1);plot(timeL,detrend_EMGL,'r'); title('detrend_EMGL');% move to zero level
subplot(2,1,2);plot(timeR,detrend_EMGR,'b'); title('detrend_EMGR');
% high pass filter % a little bit edge effect at the end
Hd = makefilter(srL,0.2,0.15,6,20,0); 
filtered_EMGL=filtfilthd(Hd,detrend_EMGL);
Hd = makefilter(srR,0.2,0.15,6,20,0); 
filtered_EMGR=filtfilthd(Hd,detrend_EMGR);
figure('units','normalized','outerposition',[0 0 1 0.6]);
subplot(2,1,1);plot(timeL,filtered_EMGL,'r'); title('filtered_EMGL');% move to zero level
subplot(2,1,2);plot(timeR,filtered_EMGR,'b'); title('filtered_EMGR');

%% Hilbert transform
% method 4 *******  spline interpolation over local maxima (best)
np = 50; % N-tap Hilbert filter
[upL,loL] = envelope(filtered_EMGL,np,'peak');
[upR,loR] = envelope(filtered_EMGR,np,'peak');
close;figure('units','normalized','outerposition',[0 0 1 0.6]);
subplot(2,1,1);
plot(timeL,filtered_EMGL,'g');
hold on; plot(timeL,upL,'r');
ylim([-1000 1000]);title(['np = ' num2str(np)]); % xlim([1000 1100]);
hold off;
subplot(2,1,2);
plot(timeR,filtered_EMGR,'g');
hold on; plot(timeR,upR,'r');
ylim([-1000 1000]);title(['np = ' num2str(np)]); % xlim([1000 1100]);
hold off;