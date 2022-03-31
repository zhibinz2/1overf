% View raw EEG
% plot(time(Startpoint:Endpoint)',samples(2:33,Startpoint:Endpoint)');
figure('units','normalized','outerposition',[0 0 1 0.3]);
plot(samples(2:33,:)');
[x, y] = ginput(2); % read two mouse clicks on the plot % x were index, y were real values
% get the proximate index
string(x)
ind1=round(x(1))
ind2=round(x(2))
% plot the section
plot(samples(2:33,ind1:ind2)');

time';
datatimes; % 
% EEGdata=samples(2:33,Startpoint:Endpoint)';


%% detrend the data
detrend_data=ndetrend(samples',1); 
figure('units','normalized','outerposition',[0 0 1 0.3]);
plot(detrend_data(ind1:ind2,2:33));

%% Broadband filter method 1 (run very slow and might crash)
% https://mikexcohen.com/lectures.html

% specify Nyquist freuqency
nyquist = sr/2;
% filter frequency band
filtbound = [0.2 50]; % Hz
% transition width
trans_width = 0.2; % fraction of 1, thus 20%
% filter order
filt_order = round(3*(sr/filtbound(1)));
% frequency vector (as fraction of Nyquist
ffrequencies  = [ 0 (1-trans_width)*filtbound(1) filtbound (1+trans_width)*filtbound(2) nyquist ]/nyquist;
% shape of filter (must be the same number of elements as frequency vector
idealresponse = [ 0 0 1 1 0 0 ];
% get filter weights
filterweights = firls(filt_order,ffrequencies,idealresponse);
% plot for visual inspection
figure(1), clf
subplot(211)
plot(ffrequencies*nyquist,idealresponse,'k--o','markerface','m')
set(gca,'ylim',[-.1 1.1],'xlim',[-2 nyquist+2])
xlabel('Frequencies (Hz)'), ylabel('Response amplitude')
subplot(212)
plot((0:filt_order)*(1000/sr),filterweights)
xlabel('Time (ms)'), ylabel('Amplitude')

% apply filter to data
nbchan=32;
detrend_data=detrend_data';
filtered_data = zeros(size(detrend_data));
for chani=1:nbchan
    filtered_data(chani,:) = filtfilt(filterweights,1,detrend_data(chani,:));
end

plot(datatimes,filtered_data);

figure(2), clf
plot(datatimes,squeeze(detrend_data(chani,:)))
hold on
plot(datatimes,squeeze(filtered_data(chani,:)),'r','linew',2)
xlabel('Time (ms)'), ylabel('Voltage (\muV)')
legend({'raw data';'filtered'})

%% filtfilt method 2 (matlab. did not specify Fs, reliable?)
% https://www.mathworks.com/help/signal/ref/filtfilt.html
% low pass filtering
d = designfilt("lowpassfir", ...
    'PassbandFrequency',0.15,'StopbandFrequency', 0.2, ...
    'PassbandRipple',1,'StopbandAttenuation',60, ...
    'DesignMethod',"equiripple");
filtered_data = filtfilt(d,detrend_data);
plot(datatimes,filtered_data);

%% filtfilthd method (hnl)
% this method seem to be dependent on EEGLAB. because EELAB reset the path for filtfilt
% eeglab; 
% which -all filtfilt.m
rmpath /home/zhibin/Documents/GitHub/matlab-archive/hnlcode/external/fieldtrip-20141007

cd /home/zhibin/Documents/GitHub/Motor_cordination/EEGanalysis/20211102
open EEG_filter_step2b.m

which filtfilthd
open filtfilthd
which makefilter
open makefilter
% Hd = makefilter(sr,Fpass,Fstop,Apass,Astop,doplot)

% low pass
Hd = makefilter(sr,50,51,6,20,1); 
filtered_data=filtfilthd(Hd,detrend_data);% after calling EEGLAB filtfilt reset to matlab default
% open /home/zhibin/Documents/GitHub/matlab-archive/hnlcode/external/fieldtrip-20141007/external/signal/filtfilt.m % this not work
% /usr/local/MATLAB/R2019a/toolbox/signal/signal/filtfilt.m % this works
figure('units','normalized','outerposition',[0 0 1 0.3]);
plot(filtered_data(ind1:ind2,2:33));

% high pass
Hd = makefilter(sr,0.2,0.15,6,20,1); % Fpass/Fstop does not seem to work under 1 Hz
filtered_data2=filtfilthd(Hd,double(filtered_data));
figure('units','normalized','outerposition',[0 0 1 0.3]);
plot(filtered_data2(ind1:ind2,2:33));


%% examine EEG quality 
cd /home/zhibin/Documents/GitHub/Motor_cordination/Data_processing_streamline

% find index of each photocell for stimulus
for i=1:length(locs)
    PhotocellInd(i)=find(datatimes==locs(i)); % locs are values in datatimes 
end
figure('units','normalized','outerposition',[0 0 1 0.3]);
plot(datatimes(ind1:ind2),filtered_data2(ind1:ind2,2:33));
hold on;
xline(datatimes(PhotocellInd(1)),'r','first tap');
xline(datatimes(PhotocellInd(end)),'r','last tap');
xline(datatimes(PhotocellInd(241)),'r','start checking power spectrum');
xline(datatimes(PhotocellInd(480)),'r','end checking power spectrum');

%% look at power spectrum
cd /home/zhibin/Documents/GitHub/Motor_cordination/EEGanalysis
open Powerspetra.m

% select a section for power spectrum
figure('units','normalized','outerposition',[0 0 1 0.3]);
plot(datatimes(PhotocellInd(241):PhotocellInd(480)),filtered_data2(PhotocellInd(241):PhotocellInd(480),2:33));

% some settings for plotting
Timeselected=datatimes(PhotocellInd(241):PhotocellInd(480));
EEGselected=filtered_data2(PhotocellInd(241):PhotocellInd(480),2:33);
figure('units','normalized','outerposition',[0 0 1 0.3]);
plot(Timeselected,EEGselected);

% fft
fcoef=fft(EEGselected);
N=size(EEGselected,1); % N=length(Timeselected);
fcoef=fcoef/N;
halfN=floor(N/2);
% df=1/T; 
% fV=[0:df:(halfN-1)*df]; 
fV=linspace(0,sr/2,halfN+1); % same df
fcoef=2*fcoef(1:halfN,:);
amplitude = abs(fcoef);
% plot one channel 
plot(fV(1:length(amplitude)),amplitude(:,1));
xlim([0 25]);
xlabel('frequency');
ylabel('amplitude (uV)');
title('Spetrum');
% plot all channels
plot(fV(1:length(amplitude)),amplitude(:,:));
xlim([0 25]);
xlabel('frequency');
ylabel('amplitude (uV)');
title('Spetra of 32 channels');
legend(labels);
% set(gca,'yscale','log')
% ylabel('log of amplitude (uV)');
% set(gca,'xscale','log')

% plot all channels separately
figure('units','normalized','outerposition',[0 0 1 1]);
for chan=1:32
    subplot(4,8,chan);
    plot(fV(1:size(amplitude,1)),amplitude(:,chan));
    xlabel('frequency');ylabel('amplitude (uV)');
    xlim([0 25]);title([labels{chan}]);ylim([0 50]);
end
suptitle('spectra of all channels')

%% Plot all erp of planning on scalp map 
figure('units','normalized','outerposition',[0 0 1 1]);
for chan=1:32
    subplot('Position',[XXPLOT(chan) YYPLOT(chan) 0.05 0.05]); % not showing, why
    % plot([1 1],[1 1],'ro');
    plot(fV(1:size(amplitude,1)),amplitude(:,chan));
    if ~isempty(find([1:30 32]==chan))
    set(gca,'XTick',[]); set(gca,'YTick',[]); 
    end
    if chan==31
        xlabel('frequency');
        ylabel('amplitude (uV)');
    end
    xlim([0 25]);title([labels{chan}]);ylim([0 50]);
end
suptitle('spectra of all channels on scalp map')
%% baseline normalization 




