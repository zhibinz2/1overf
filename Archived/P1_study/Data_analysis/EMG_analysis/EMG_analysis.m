% View raw EMG
close;figure('units','normalized','outerposition',[0 0 1 0.3]);
plot(time,EMG,'k');
%% EMG preprocessing 
% detrend
detrend_data=detrend(EMG,2); % even better
hold on;
plot(time,detrend_data,'r'); % move to zero level
legend('raw-EMG','detrend-EMG');
% high pass filter % a little bit edge effect at the end
Hd = makefilter(sr,0.2,0.15,6,20,0); 
filtered_data=filtfilthd(Hd,detrend_data);
plot(time,filtered_data,'b');
legend('raw-EMG','detrend-EMG','high-pass-EMG');

%% Hilbert transform
% method 1 ******
EnveEMG=abs(hilbert(filtered_data));
close;figure('units','normalized','outerposition',[0 0 1 0.3]);
plot(time,filtered_data,'b');
hold on; plot(time',EnveEMG,'r');

% method 2 ****** adjust Hilbert filter length
fl = 50*sr;
[up1,lo1] = envelope(filtered_data,fl,'analytic');
close;figure('units','normalized','outerposition',[0 0 1 0.3]);
plot(time,filtered_data,'b');
hold on; plot(time,up1,'r');

% method 3 ******** adjust length of window
wl = 100;
[up1,lo1] = envelope(filtered_data,wl,'rms');
close;figure('units','normalized','outerposition',[0 0 1 0.3]);
plot(time,filtered_data,'b');
hold on; plot(time,up1,'r');

% method 4 *******  spline interpolation over local maxima (taking longer)
np = 500;
[up1,lo1] = envelope(filtered_data,np,'peak');
close;figure('units','normalized','outerposition',[0 0 1 0.3]);
plot(time,filtered_data,'b');
hold on; plot(time,up1,'r');
ylim([-1000 1000]);title(['np = ' num2str(np)]);

%%
