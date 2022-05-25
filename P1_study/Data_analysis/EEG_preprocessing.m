%% View raw EEG (skip, let's look at Photocell,reponses, EMG and EEG all together)
% plot(time(Startpoint:Endpoint)',samples(2:33,Startpoint:Endpoint)');
figure('units','normalized','outerposition',[0 0 1 0.3]);
% plot(samples(2:33,:)');
plot(samples(1:32,:)');% hold on; xline(x(1),'r');xline(x(2),'r');
[x, y] = ginput(2); % read two mouse clicks on the plot % x were index, y were real values
% get the proximate index
string(x)
ind1=round(x(1))
ind2=round(x(2))
xline(ind1,'b');xline(ind2,'b');
hold off;
% ind1=1;
% ind2=length(time);

% plot the section
% plot(samples(2:33,ind1:ind2)');
% figure('units','normalized','outerposition',[0 0 1 0.3]);
% plot(samples(1:32,ind1:ind2)');

% time';
% datatimes; % 
% EEGdata=samples(2:33,Startpoint:Endpoint)';

%% View Photocell, button presses, EMG, EEG together
figure('units','normalized','outerposition',[0 0 1 0.6]);
subplot(3,1,1);
plot(Photocell);
ylabel('photocell signal');title('photocell signal');
hold on;
plot(BottonPresTimeInd,ones(1,length(BottonPresTimeInd)),'r.'); title('botton press');
legend({'photocell signal','botton press'});
subplot(3,1,2);
plot(EMG,'k'); title('EMG'); % EMG channel
subplot(3,1,3);
plot(EEG);
xlabel('time (sample)');title('EEG');

[x, ~] = ginput(2); % read two mouse clicks on the plot % x were index, y were real values
% get the proximate index
string(x)
ind1=round(x(1))
ind2=round(x(2))

close;
figure('units','normalized','outerposition',[0 0 1 0.6]);
subplot(3,1,1);
plot(Photocell);
ylabel('photocell signal');title('photocell signal');
hold on;
plot(BottonPresTimeInd,ones(1,length(BottonPresTimeInd)),'r.'); title('botton press');
legend({'photocell signal','botton press'});
hold on; xline(ind1,'r');xline(ind2,'r');hold off;

subplot(3,1,2);
plot(EMG,'k'); title('EMG'); % EMG channel
hold on; xline(ind1,'r');xline(ind2,'r');hold off;

subplot(3,1,3);
plotx(EEG);
xlabel('time (sample)');title('EEG');
hold on; xline(ind1,'r');xline(ind2,'r');hold off;


%% cut off artifact at the beginning and end, then extended the two sides to the original length
EEG2=[repmat(EEG(ind1,:),ind1-1,1); EEG(ind1:ind2,:); repmat(EEG(ind2,:),length(time)-ind2,1)];
figure('units','normalized','outerposition',[0 0 1 0.3]);
plot(EEG2);

EMG2=[repmat(EMG(ind1),ind1-1,1); EMG(ind1:ind2); repmat(EMG(ind2),length(time)-ind2,1)];
figure('units','normalized','outerposition',[0 0 1 0.3]);
plot(EMG2);

%% %% View all together again
figure('units','normalized','outerposition',[0 0 1 0.6]);
subplot(3,1,1);
plot(Photocell);
ylabel('photocell signal');title('photocell signal');
hold on;
plot(BottonPresTimeInd,ones(1,length(BottonPresTimeInd)),'r.'); title('botton press');
legend({'photocell signal','botton press'});
hold on; xline(ind1,'r');xline(ind2,'r');hold off;

subplot(3,1,2);
plot(EMG2,'k'); title('EMG'); % EMG channel
hold on; xline(ind1,'r');xline(ind2,'r');hold off;

subplot(3,1,3);
plotx(EEG2);
xlabel('time (sample)');title('EEG');
hold on; xline(ind1,'r');xline(ind2,'r');hold off;

%% %% select a section and Zoom in all together again
[x, y] = ginput(2); % read two mouse clicks on the plot % x were index, y were real values
% get the proximate index
ind1=round(x(1))
ind2=round(x(2))

close;
figure('units','normalized','outerposition',[0 0 1 0.6]);
subplot(3,1,1);
plot(Photocell);
ylabel('photocell signal');title('photocell signal');
hold on;
plot(BottonPresTimeInd,ones(1,length(BottonPresTimeInd)),'r.'); title('botton press');
legend({'photocell signal','botton press'});
hold on; xline(ind1,'r');xline(ind2,'r');hold off;
xlim([ind1 ind2]);

subplot(3,1,2);
plot(EMG2,'k'); title('EMG'); % EMG channel
hold on; xline(ind1,'r');xline(ind2,'r');hold off;
xlim([ind1 ind2]);

subplot(3,1,3);
plotx(EEG2);
xlabel('time (sample)');title('EEG');
hold on; xline(ind1,'r');xline(ind2,'r');hold off;
xlim([ind1 ind2]);

%% detrend the data (no padding needed)
% lab's detrend function
% detrend_data=ndetrend(samples(2:33,:)',1); 
% detrend_data=ndetrend(samples(1:32,:)',1); 
% detrend_data=ndetrend(EEG,1); 
% % mathwork's detrend function
% detrend_data=detrend(EEG,0); % not good
% detrend_data=detrend(EEG,1); % similar to ndetrend
detrend_data=detrend(EEG2,2); % even better
figure('units','normalized','outerposition',[0 0 1 0.3]);
% plot(detrend_data(ind1:ind2,2:33));
plot(detrend_data);
% plot(detrend_data(ind1:ind2,:));

% [x, y] = ginput(2); % read two mouse clicks on the plot % x were index, y were real values
% get the proximate index
% string(x)
% ind1=round(x(1))
% ind2=round(x(2))
% hold on;
% xline(ind1,'b');xline(ind2,'b');
% xlim(round([x(1) x(2)]));ylim(round([y(1) y(2)]));
% hold off;
%% Broadband filter method 1 (run very slow and might crash) (skip)
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

%% filtfilt method 2 (matlab. did not specify Fs, reliable?) (skip)
% https://www.mathworks.com/help/signal/ref/filtfilt.html
% low pass filtering
d = designfilt("lowpassfir", ...
    'PassbandFrequency',0.15,'StopbandFrequency', 0.2, ...
    'PassbandRipple',1,'StopbandAttenuation',60, ...
    'DesignMethod',"equiripple");
filtered_data = filtfilt(d,detrend_data);
plot(datatimes,filtered_data);

%% filtfilthd method (hnl) low pass first (skip)
% this method seem to be dependent on EEGLAB. because EELAB reset the path for filtfilt
% eeglab; 
% which -all filtfilt.m
rmpath /home/zhibin/Documents/GitHub/matlab-archive/hnlcode/external/fieldtrip-20141007
cd  /home/zhibin/Documents/GitHub/matlab-archive
cd /usr/local/MATLAB/R2019a/toolbox/signal/signal
which filtfilt.m

% cd /home/zhibin/Documents/GitHub/Motor_cordination/EEGanalysis/20211102
% open EEG_filter_step2b.m

% which filtfilthd
% open filtfilthd
% which makefilter
% open makefilter
% Hd = makefilter(sr,Fpass,Fstop,Apass,Astop,doplot)

% low pass 
% (this will create short edge artifact )
Hd = makefilter(sr,50,51,6,20,0);  % Hd = makefilter(sr,50,51,6,20,1); 
filtered_data=filtfilthd(Hd,detrend_data);% after calling EEGLAB filtfilt reset to matlab default
% open /home/zhibin/Documents/GitHub/matlab-archive/hnlcode/external/fieldtrip-20141007/external/signal/filtfilt.m % this not work
% /usr/local/MATLAB/R2019a/toolbox/signal/signal/filtfilt.m % this works
figure('units','normalized','outerposition',[0 0 1 0.3]);
% plot(filtered_data(ind1:ind2,2:33));
% plot(filtered_data(ind1:ind2,:));
plot(filtered_data(:,:));
ylim([-500 500]);

% cut off artifact at the beginning and end
[x, y] = ginput(2); % read two mouse clicks on the plot % x were index, y were real values
% get the proximate index
ind1=round(x(1));
ind2=round(x(2));
filtered_data=filtered_data(ind1:ind2,:);
time=time(ind1:ind2);
figure('units','normalized','outerposition',[0 0 1 0.3]);
plot(filtered_data(:,:));
ylim([y(1) y(2)]);;

% add padding
% padding=zeros(size(filtered_data));
% filtered_data=cat(1,padding,filtered_data,padding);

% high pass 
% (very large artifact if not cut off edges after low pass filter, even padding won't slove it)
% (no edge artifact if cut off edges after low pass filter, no padding needed)
Hd = makefilter(sr,0.2,0.15,6,20,0); % Hd = makefilter(sr,0.2,0.15,6,20,1); xlim([0 2]);% Fpass/Fstop does not seem to work under 1 Hz
filtered_data2=filtfilthd(Hd,double(filtered_data));
figure('units','normalized','outerposition',[0 0 1 0.3]);
% plot(filtered_data2(ind1:ind2,2:33));
% plot(filtered_data2(ind1:ind2,:));
plot(filtered_data2(:,:));ylim([-5000 5000]);
hold on;xline(x(1),'r');xline(x(2),'r');hold off;

%% filtfilthd method (hnl) high pass first
% high pass (no paddings needed)
Hd = makefilter(sr,0.2,0.15,6,20,0); 
filtered_data1=filtfilthd(Hd,detrend_data);
figure('units','normalized','outerposition',[0 0 1 0.3]);
plot(filtered_data1);

% add padding
padding=zeros(round(size(filtered_data1,1)/10), size(filtered_data1,2));
filtered_data2=cat(1,padding,filtered_data1,padding);
figure('units','normalized','outerposition',[0 0 1 0.3]);
plot(filtered_data2);

% low pass 
% (this will create short edge artifact)
% (if added zero paddings, edge artifact disappear)
% (remove existing filtered_data variable from workspace might fasten)
tic
filtered_data3=[];
Hd = makefilter(sr,50,51,6,20,0);  
filtered_data3=filtfilthd(Hd,filtered_data2);
figure('units','normalized','outerposition',[0 0 1 0.3]);
plot(filtered_data3);
% ylim([-100 100]);
toc

% remove padding
filtered_data4=filtered_data3((size(padding,1)+1):(size(padding,1)+size(detrend_data,1)),:);
figure('units','normalized','outerposition',[0 0 1 0.3]);
plotx(filtered_data4);
% ylim([-100 100]);

 
clearvars samples filtered_data1 filtered_data2 filtered_data3
% clear samples filtered_data1 filtered_data2 filtered_data3
% clear samples
%% select a section to examine EEG quality 
cd /home/zhibin/Documents/GitHub/Motor_cordination/Data_processing_streamline

% find index of each photocell for stimulus
for i=1:length(locs)
    PhotocellTimeInd(i)=find(time==locs(i)); % locs are values in time 
end

figure('units','normalized','outerposition',[0 0 1 0.3]);
plot(time,Photocell);
hold on; xline(time(PhotocellTimeInd(5)),'m');
xline(time(PhotocellTimeInd(245)),'m');
xline(time(PhotocellTimeInd(485)),'m');
% ylim([-1000 1000]);
% xlim([ind1 ind2]);
% hold on;xline(x(1),'r');xline(x(2),'r');
% xline(time(PhotocellInd(1)),'r','first photocell');
% xline(time(PhotocellInd(end)),'r','last photocell');
% xline(datatimes(PhotocellInd(241)),'r','start checking power spectrum');
% xline(datatimes(PhotocellInd(480)),'r','end checking power spectrum');

%% plot Photocell, button presses, EMG, EEG together and zoom in together
figure('units','normalized','outerposition',[0 0 1 0.6]);
subplot(3,1,1);
plot(Photocell);
ylabel('photocell signal');title('photocell signal');
hold on;
plot(BottonPresTimeInd,ones(1,length(BottonPresTimeInd)),'r.'); title('botton press');
legend({'photocell signal','botton press'});
subplot(3,1,2);
plot(EMG2,'k'); title('EMG'); % ylim([-5000 5000]);% EMG channel
subplot(3,1,3);
plot(filtered_data4);
xlabel('time (sample)');title('EEG');% ylim([-500 500]);

[x, y] = ginput(2); % read two mouse clicks on the plot % x were index, y were real values
% get the proximate index
string(x)
ind1=round(x(1))
ind2=round(x(2))
hold on; xline(ind1,'r');xline(ind2,'r');hold off;

% close;
% figure('units','normalized','outerposition',[0 0 1 0.6]);
subplot(3,1,1);
plot(Photocell');
ylabel('photocell signal');title('photocell signal');
hold on;
plot(BottonPresTimeInd,ones(1,length(BottonPresTimeInd)),'r.'); title('botton press');
legend({'photocell signal','botton press'});
xlim([ind1 ind2]);

subplot(3,1,2);
plot(EMG2,'k'); title('EMG'); % EMG channel
xlim([ind1 ind2]);

subplot(3,1,3);
plotx(filtered_data4);
xlabel('time (sample)');title('EEG');
ylim([-500 500]);
xlim([ind1 ind2]);


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
    xlim([0 25]);title([labels{chan}]);ylim([0 15]);
end
suptitle('spectra of all channels')

%% Plot on scalp map for the spectra (for my own examing)
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
    xlim([0 25]);title([labels{chan}]);ylim([0 15]);
end
suptitle('spectra of all channels on scalp map')
%% run ICA to remove eye blinds
% look at data before ICA
plot(filtered_data4);
% previous code for ICA
cd /home/zhibin/Documents/GitHub/Motor_cordination/Data_processing_streamline
open combine_sessions.m
cd /home/zhibin/Documents/GitHub/Motor_cordination/EEGanalysis/20211102
open RunICA_step6b.m

% run ICA
tic
[icasig, A, W] = fastica(filtered_data4');
toc

% Plot ICA component
for i=1:size(A,2)
    SqA(i)=sumsqr(A(:,i));
end
figure;
plot(1:size(A,2),SqA,'ro');ylabel('sum of square of column in A');xlabel('ICs');
[B,I]=sort(SqA,'descend');

% topoplot to examine them
cd /home/zhibin/Documents/GitHub/Motor_cordination/1_over_f/data_analysis/channels_info
load('chaninfo.mat')
ComponentsExam=I(1:10);
figure;
for i=1:length(ComponentsExam)
    subplot(5,2,i);
    topoplot(A(:,ComponentsExam(i)),chaninfo,'nosedir','+X');title(['component' num2str(ComponentsExam(i))]);colorbar;
end

% Display signal before ICA 
figure('units','normalized','outerposition',[0 0 1 0.5]);
plot(time,filtered_data4);
hold on; hold off;
title('EEG Singal Before ICA');

% Display ICA signals
figure('units','normalized','outerposition',[0 0 1 1]);
subplot(2,1,1);
plot(time,icasig');
hold on;hold off;
title('all IC time series');ylabel('uV');
subplot(2,1,2);
plot(time,icasig(ComponentsExam,:)');
hold on;hold off;
legend(strsplit(num2str(ComponentsExam)));
title('first 10 IC time series');xlabel('time samples');ylabel('uV');

% Calculate Correlation
% FP1 and FP2 are channel 1 and 3;
% compute correlation between FP1 and FP2;
[RHO1,PVAL1] = corr(filtered_data4(:,1),icasig');
figure
subplot(1,2,1);
plot(1:length(RHO1),RHO1,'bo');hold on;
yyaxis left;ylabel('correlation coefficient','color','b');set(gca,'ycolor','b');
plot(1:length(RHO1),PVAL1,'ro');
yyaxis right; ylabel('p value','color','r');set(gca,'ycolor','r'); xlabel('ICs');
legend({'correlation coefficient','p-values '});title('correlation with FP1');hold off;
subplot(1,2,2);
[RHO3,PVAL3] = corr(filtered_data4(:,3),icasig');
plot(1:length(RHO1),RHO3,'bo');hold on;
yyaxis left;ylabel('correlation coefficient','color','b');set(gca,'ycolor','b');
plot(1:length(RHO1),PVAL3,'ro');
yyaxis right; ylabel('p value','color','r');set(gca,'ycolor','r'); xlabel('ICs');
legend({'correlation coefficient','p-values '});title('correlation with FP2');hold off;

[B1,I1]=sort(abs(RHO1),'descend');[B3,I3]=sort(abs(RHO3),'descend');
ComponentsExam=unique([I1(1) I3(1)]);

figure;
for i=1:length(ComponentsExam)
    subplot(length(ComponentsExam),1,i);
    topoplot(A(:,ComponentsExam(i)),chaninfo,'nosedir','+X');title(['component' num2str(ComponentsExam(i))]);colorbar;
end

%% Plot ICs Topoplot Spetrogram ERP AND Power spectra,  (Similar to EEGLAB) (skip)
ComponentsExam=I(1);

for i=1:length(ComponentsExam)
    ComponentPick=ComponentsExam(i);
    
    % Put together trial matrics
    icasig_trials=[];
    for i=1:length(goodepochs)
        icasig_trials(:,1,i)=icasig(ComponentPick,IndStart(i):IndEnds(i))';
    end

    % Just to check
    goodchans;

    figure('units','normalized','outerposition',[0 0 0.38 0.85]);
    % Topoplot
    subplot('Position',[0.05 0.55 0.4 0.4]);
    topoplot(A(:,ComponentPick),test,'nosedir','+Y');
    title(['component' num2str(ComponentPick)]);
    colorbar;

    % Plot Spectrogram
    subplot('Position',[0.55 0.65 0.4 0.3]);
    ColorLim=5;
    imagesc(squeeze(icasig_trials)');colormap jet; colorbar ; ylabel('Shorted Trials'); 
    caxis([-1*ColorLim ColorLim]);
    xline(500,'k');xline(1000,'k');
    title(['Component ' num2str(ComponentPick)]);

    % Plot IC ERP
    subplot('Position',[0.55 0.55 0.3 0.07]);
    plot(1:size(icasig_trials,1),mean(icasig_trials,3)); 
    xline(500,'k');xline(1000,'k');xlabel('time');

    % Plot Powerspetra
    subplot('Position',[0.1 0.05 0.8 0.4]);
    rate=Fs;maxfreq=50;
    [pow,freqs,df,eppow,corr,cprod,fcoef] = allspectra(icasig_trials,rate,maxfreq);
    plot(freqs,pow);xlabel('frequecy (Hz)');ylabel('Magnitude');
end

%% Deside which components to remove and mix back the signal and display
ComponentRemove=I(1);

A(:,ComponentRemove)=0; icasig(ComponentRemove,:)=0;

mixedsig=A*icasig;

% Plot before and after
figure('units','normalized','outerposition',[0 0 1 1]);
subplot(2,1,1);
% before ICA
plotx(time,filtered_data4);
hold on;hold off;
title('EEG Singal Before ICA');
% after ICA
subplot(2,1,2);
plotx(time,mixedsig);
hold on;hold off;
title('Mixed Signal with ICs removed');

%% segment EEG according to conditions
separationsTimeInd; % indices of frist press in each condition in time

PhotocellTimeInd;
figure('units','normalized','outerposition',[0 0 1 0.3]);
plot(time,Photocell);
hold on; xline(time(PhotocellTimeInd(5)),'m');
xline(time(PhotocellTimeInd(245)),'m');
xline(time(PhotocellTimeInd(485)),'m');
xline(time(PhotocellTimeInd(724)),'m');

% segment EEG into trials
EEG_eye_open_resting=mixedsig(:,PhotocellTimeInd(1):PhotocellTimeInd(2))';
EEG_eye_close_resting=mixedsig(:,PhotocellTimeInd(3):PhotocellTimeInd(4))';
allPerm
EEG_syncopation=mixedsig(:,PhotocellTimeInd(5):PhotocellTimeInd(244))';
EEG_randomization=mixedsig(:,PhotocellTimeInd(245):PhotocellTimeInd(484))';
EEG_synchronization=mixedsig(:,PhotocellTimeInd(485):PhotocellTimeInd(724))';

cd /ssd/zhibin/1overf/20220331
save('segmented_clean_Data.mat','EEG_eye_open_resting','EEG_eye_close_resting',...
    'EEG_syncopation','EEG_randomization','EEG_synchronization');



EEG_trials=zeros(size(mixedsig'));
events=[];
for i=1:3
    eventstoadd=[PhotocellTimeInd(5+240*(i-1)) PhotocellTimeInd(5+240*i-1)];
    events=[events;eventstoadd];
end
events(:,2)-events(:,1)
EEG_trials=zeros(size(mixedsig'));
EEG

%% baseline normalization 




