%% load TMSi data
% open /home/zhibin/Documents/GitHub/TMSi_SAGA_matlab/SAGA_interface/testdatatryout2022.m

clear;close all;

% Step 1: Open Poly5 file.
addpath /home/zhibin/Documents/GitHub/TMSi_SAGA_matlab % hnlb
cd /ssd/zhibin/1overf/20220324
d = TMSiSAGA.Poly5.read('/ssd/zhibin/1overf/20220324/20220324.Poly5');

samples=d.samples;
sr=d.sample_rate;
channels=d.channels;

% Organize the channels

% unit_names=cell(length(channels),1);
% altertive_names=cell(length(channels),1);
% for i=1:length(channels)
%     altertive_names(i)={channels{i}.alternative_name};
%     unit_names(i)={channels{i}.unit_name};
% end

numbers=num2str([1:length(channels)]');
labels=strings(length(channels),1);
units=strings(length(channels),1);
for i=1:length(channels)
    labels(i)=channels{i}.alternative_name;
    units(i)=channels{i}.unit_name;
end
channels_info=table(numbers,labels,units)

% Create time stamps
num2str(d.time)
time=[1/sr:1/sr:d.time];

% Plot channels of Key presses, photocells, EMG
plot(samples(54,:),'ro'); % trigger = key presses = 223 (255-2^5)
unique(samples(54,:))
hold on; 
plot(samples(38,:),'b'); %  ISO aux = analog
plot(samples(34,:),'k'); % EMG channel

% Save Channels of presses, photocells, EMG
BottonPres=samples(54,:);
Photocell=samples(38,:);
EMG=samples(34,:);

% View EEG
plot(time',samples(2:33,:)');

%% load TMSi for EEGLAB
% open /home/zhibin/Documents/GitHub/TMSi_SAGA_matlab/SAGA_interface/testdatatryout2022.m

cd /ssd/zhibin/1overf/20220324
[file, pathname] = uigetfile({'*.Poly5';'*.S00';'*.TMS32'},'Pick your file');
[path,filename,extension] = fileparts(file);
open Poly5toEEGlab

% Step 2: Plot a single channel.
plot((0:(d.num_samples - 1)) / d.sample_rate, d.samples(2, :));

% Step 3: Save dataset in the same directory as the *.Poly5 file.
% Open EEGlab
eeglab

% Load TMSi's channel location file
load('EEGChannels64TMSi.mat', 'ChanLocs');

% Transform data to eeglab. 
eegdataset = toEEGLab(d, ChanLocs);

% Save dataset in the same directory as the *.Poly5 file.
pop_saveset(eegdataset,'filename',filename,'filepath',pathname)
disp(['Data saved as EEGlab dataset (.set) in this folder: ',pathname])

%% load stimulus data
cd /home/zhibin/Documents/GitHub/1overf/stimulus_data_storage
load('20220324.mat')

%% Covert BottonPres and Photocell signals to time points
Photocell=samples(38,:);
BottonPres=samples(54,:);
unique(BottonPres) 

time;

open /home/zhibin/Documents/GitHub/Motor_cordination/EEGanalysis/20211102/organize_photocells_step2.m

% for photocell
% view the time course of photocell signals
figure('units','normalized','outerposition',[0 0 1 0.3]);
plot(Photocell');xlabel('time');ylabel('photocell signal');
% plot EEG on top
% hold on; plot(time',samples(2:33,:)'); % it zoom out the phtocell amplitude, too small to see
% click and select the start and end point for peak extraction
[x, y] = ginput(2); % read two mouse clicks on the plot % x were index, y were real values
Startpoint=round(x(1));Endpoint=round(x(2)); % Startpoint and Endpoint are sample number or index in time

% replace the beginning and end with baseline value
Photocell(1:Startpoint)=mean(y);Photocell(Endpoint:end)=mean(y); % plot(Photocell');

analog1data=Photocell; 
datatimes=time;
plot(datatimes,analog1data,'b'); 

% Examine peaks detection in analog1
Halfhigh1=3/4*(max(analog1data)-min(analog1data)); % value for 'MinPeakProminence'
% Check if need to adjust the Halfhigh cutoff
close;figure;
findpeaks(analog1data,datatimes,'MinPeakProminence',Halfhigh1,'Annotate','extents');
yline(Halfhigh1,'m','MinPeakProminence');

% locate the trial sessions % pks=value of the peak % locs=time of the peak
[pks,locs] = findpeaks(analog1data,datatimes,'MinPeakProminence',Halfhigh1,'Annotate','extents');

% examine pks and locs (both are values of analog1data and datatimes, not indices)
% i=1;
% find(analog1data==pks(i)) % return the index of peak in analog1data
% find(datatimes==locs(i)) % return the same "index" in datatimes (ie, i=1, index=4226)
% so, infer from above, beacause the above "find values" for the same "index" works
% pks are values in analog1data
% locs are values in datatimes ((ie, i=1, value=4225))

% figure;plot(locs,'bo');
% figure;bar(locs);
% figure;plot(pks,'bo');
% figure;bar(pks);
% figure;plot(locs,pks,'bo');ylim([min(analog1data) max(analog1data)]);

% locsDurations=diff(locs);% time durations between peaks
% close;figure;
% plot(locsDurations,'ro');% look at distribution of these durations   
% xlabel('each peak');ylabel('locsDurations = Between-peaks time duration (ms)');
% Lowcutoff=5*mean(locsDurations);% cutoff standard of the between-peak duarations to separate between trials
% Highcutoff=180*mean(locsDurations);
% hold on; yline(Lowcutoff,'m--','lowcut');hold off; % examine the cutoff line
% hold on; yline(Highcutoff,'m--','highcut');hold off;

% for botton presses
botton1data=BottonPres(Startpoint:Endpoint);
% look for values other than 0 and 255
plot(datatimes,botton1data,'bo'); % view time course of botton press values
PresInd=find(botton1data~=255); % extract Index of real key presses in the values
plot(PresInd,ones(1,length(PresInd)),'ro'); % look at the above Index
threshold = NumFramesInterval*ifi*sr/4; % determine a threshold of key press interval
BottonPresTimeInd=PresInd(find([1 diff(PresInd)>threshold])); % exact index of key press onset in datatimes
% examine key press interval
figure;
plot(1:length(diff(BottonPresTimeInd)), diff(BottonPresTimeInd),'ro'); 
xlabel('key press');ylabel('Intervals (ms)');

%% Compare photocell amd botton presses timing
PhotocellTime=locs; % locs are values in datatimes
BottonPressTime=datatimes(BottonPresTimeInd);

% view and compare
close;figure;
plot(PhotocellTime',ones(length(PhotocellTime),1),'bo');
hold on;
plot(BottonPressTime',ones(length(BottonPressTime),1),'ro');
xlabel('time');title('Photocell / Botton press indicator');
hold off;
legend('Photocell','Botton press');
% I ended up doing synchronization all the time


%% Compute errors
cd /ssd/zhibin/1overf/20220324

% compute the error
SynError=[];
% find the time difference with the closest botton press 
for i=1:480
    [minValue,closetIndex]=min(abs(BottonPressTime-PhotocellTime(i))); % closetIndex in BottonPressTime
    SynError(i)=PhotocellTime(i)-BottonPressTime(closetIndex)
end
plot(SynError,'ro');title('synchronization error');xlabel('taps');ylabel('timing error (s)');

RanError=[];
% find the time difference with the next botton press
for i=491:length(PhotocellTime)
    diffInd=(find((BottonPressTime-PhotocellTime(i))>0));
    RanError(i)=BottonPressTime(diffInd(1))-PhotocellTime(i);
end
RanError=RanError(491:end);
plot(RanError,'ro'); title('randomization error');xlabel('taps');ylabel('timing error (s)');

%% auto correlation and spectrum
open /ssd/zhibin/1overf/20220324/explore_1_over_f.m

SupTitles={'Synchronization' 'Randomization'};

for i=1:2 % i=2;
    if i==1
        y=SynError;
    else
        y=RanError;
    end
    
    F=fft(y);
    p=abs(F).^2;
    Fs=1/2;% assuming a sampling frequency of 0.5Hz (2 second per sample)
    freqs=linspace(0,Fs/2,length(y)/2+1);

    % plot(log10(freqs),log10(p(1:length(freqs))));
    % xlabel('log of frequency(Hz)');
    % ylabel('log of power');
    % title('{\bf LogLogPlot}');

    % loglog(x,y)
    % grid on


    % 3 plots
    figure('units','normalized','outerposition',[0 0 0.4 0.6]);
    
    subplot(1,3,1);
    plot(y);
    xlabel('Trial');ylabel('Error(sec.)');title('A. Time Series');
    
    subplot(1,3,2);
    autocorr(y,'NumLags',200,'NumSTD',2);
    % [acf,lags,bounds] = autocorr(RTs(ConditionExplore),'NumLags',25,'NumSTD',2);
    xlabel('Lag');ylabel('Corr.');title('B. Autocorr.');
    
    subplot(1,3,3);
    xx=log10(freqs(2:end));yy=log10(p(2:length(freqs)));
    plot(xx,yy);hold on;
    tbl=table(xx',yy');
    mdl=fitlm(tbl,'linear');
    plotAdded(mdl);
    xlabel('Log(f)');ylabel('Log(power)');title('C. Spectrum');
    
    suptitle(char(SupTitles(i)));
end


