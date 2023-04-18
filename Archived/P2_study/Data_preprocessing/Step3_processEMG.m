%% plot raw EMG
timeL;samplesL;TRIGGERindL;srL;channels_infoL;labelsL;
timeR;samplesR;TRIGGERindR;srR;channels_infoR;labelsR;
% find EMG chan (BIP 01, BIP 02)
EMGindL=find(labelsL=='BIP 01');
EMGindR=find(labelsR=='BIP 01');
% EMGindR=find(labelsR=='BIP 02');
figure;
subplot(2,1,1);plot(samplesL(EMGindL,:)');
subplot(2,1,2);plot(samplesR(EMGindR,:)');


%% Extract EMG
EMGL=samplesL(EMGindL,:)';
EMGR=samplesR(EMGindR,:)';

%% EMG preprocessing 
% detrend
detrend_EMGL=detrend(EMGL,2); % even better
detrend_EMGR=detrend(EMGR,2); 
figure('units','normalized','outerposition',[0 0 1 0.6]);
subplot(2,1,1);plot(timeL,detrend_EMGL,'r'); title('detrend-EMGL');% move to zero level
subplot(2,1,2);plot(timeR,detrend_EMGR,'b'); title('detrend-EMGR');
% high pass filter % a little bit edge effect at the end
Hd = makefilter(srL,0.2,0.15,6,20,0); 
filtered_EMGL=filtfilthd(Hd,detrend_EMGL);
Hd = makefilter(srR,0.2,0.15,6,20,0); 
filtered_EMGR=filtfilthd(Hd,detrend_EMGR);
figure('units','normalized','outerposition',[0 0 1 0.6]);
subplot(2,1,1);plot(timeL,filtered_EMGL,'r'); title('filtered-EMGL');% move to zero level
subplot(2,1,2);plot(timeR,filtered_EMGR,'b'); title('filtered-EMGR');

%% Hilbert transform - synchronization
% method 4 *******  spline interpolation over local maxima (best)
np = 50; % N-tap Hilbert filter
[upL,loL] = envelope(filtered_EMGL,np,'peak');
[upR,loR] = envelope(filtered_EMGR,np,'peak');

close;figure('units','normalized','outerposition',[0 0 1 0.6]);
subplot(2,1,1);
plot(filtered_EMGL,'g');
hold on; plot(upL,'r');
ylim([-1000 1000]);title(['np = ' num2str(np)]); % xlim([1000 1100]);
hold off;
subplot(2,1,2);
plot(filtered_EMGR,'g');
hold on; plot(upR,'r');
ylim([-1000 1000]);title(['np = ' num2str(np)]); % xlim([1000 1100]);
hold off;
% selet two point on x axis to zoom in
[x, ~] = ginput(2); % read two mouse clicks on the plot % x were index, y were real values
% get the proximate index
string(x)
ind1=round(x(1))
ind2=round(x(2))

% zoom in and plot again
figure('units','normalized','outerposition',[0 0 1 0.6]);
subplot(2,1,1);
plot(filtered_EMGL(ind1:ind2),'g');
hold on; plot(upL(ind1:ind2),'r');
ylim([-1000 1000]);title(['np = ' num2str(np)]); % xlim([1000 1100]);
hold off;
subplot(2,1,2);
plot(filtered_EMGR(ind1:ind2),'g');
hold on; plot(upR(ind1:ind2),'r');
ylim([-1000 1000]);title(['np = ' num2str(np)]); % xlim([1000 1100]);
hold off;

%% extract EMG variables
filtered_EMGL;filtered_EMGR;
loL;loR; % Hilbert envelop


