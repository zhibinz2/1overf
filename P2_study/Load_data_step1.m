%% load TMSi data
clear;close all;
addpath /home/zhibin/Documents/GitHub/Motor_cordination/1_over_f/data_analysis/2P_testing;
cd /ssd/zhibin/1overf/20220515_2P
cd /ssd/zhibin/1overf/20220517_2P
cd /ssd/zhibin/1overf/20220518_2P
% add "https://github.com/zhibinz2/TMSi_SAGA_matlab/" to your path

[EEGfileNameL]=uigetfile('*.Poly5');% select the Left player EEG
Path_filenameL=[pwd '/' EEGfileNameL];
[timeL,samplesL,TRIGGERindL,srL,channels_infoL] = LoadTMSi(Path_filenameL);

[EEGfileNameR]=uigetfile('*.Poly5');% select the Right player EEG
Path_filenameR=[pwd '/' EEGfileNameR];
[timeR,samplesR,TRIGGERindR,srR,channels_infoR] = LoadTMSi(Path_filenameR);

%% Extract Trigger channel info
TriggersL=samplesL(TRIGGERindL,:)';unique(TriggersL)
TriggersR=samplesR(TRIGGERindR,:)';unique(TriggersR)

%% extract time points for button presses 
% Left player
PresIndL=unique([find(TriggersL == 239);find(TriggersL == 235);find(TriggersL == 111 )]); % extract Index of real key presses in the values
plot(PresIndL,ones(1,length(PresIndL)),'ro'); % look at the above Index (one press produced several indices)
% ifi=1/100;waitframes=2;sr=2000;
% determine a threshold of numbers of frames in the button press interval
threshold = 1*1.5; % interval larger than 1 consecutive samples
% button presses recorded via DP port into stimtracker is sampled at 2000 Hz
% feedback via photocell to stimtracker has minial interval every other waitframes (minmal of 1000 Hz)
BottonPresTimeIndL=PresIndL(find([0 diff(PresIndL')]>threshold)); % exact index of key press onset in datatimes (reduce several indices into one)
% create a time series that assign botton presses as 1, all other as 0
BottonPresTimeL01=zeros(size(timeL));
BottonPresTimeL01(BottonPresTimeIndL)=1;
plot(BottonPresTimeL01,'r.')
sum(BottonPresTimeL01) % num of button presses


% Right player
PresIndR=unique([find(TriggersR == 239);find(TriggersR == 235);find(TriggersR == 111 )]); % extract Index of real key presses in the values
plot(PresIndR,ones(1,length(PresIndR)),'bo'); % look at the above Index (one press produced several indices)
% ifi=1/100;waitframes=2;sr=2000;
% determine a threshold of numbers of frames in the button press interval
threshold = 1*1.5; % interval larger than 1 consecutive samples
% button presses recorded via DP port into stimtracker is sampled at 2000 Hz
% feedback via photocell to stimtracker has minial interval every other waitframes (minmal of 1000 Hz)
BottonPresTimeInd=PresIndR(find([0 diff(PresIndR')]>threshold)); % exact index of key press onset in datatimes (reduce several indices into one)
% create a time series that assign botton presses as 1, all other as 0
BottonPresTimeR01=zeros(size(timeR));
BottonPresTimeR01(BottonPresTimeInd)=1;
plot(BottonPresTimeR01,'b.')
numPresR=sum(BottonPresTimeR01) % num of button presses

%% extract time points for feedbacks from the other player
% Left Player
FeedbIndL=unique([find(TriggersL == 251);find(TriggersL == 235)]); % extract Index of real key presses in the values
plot(FeedbIndL,ones(1,length(FeedbIndL)),'bo'); % look at the above Index (one press produced several indices)
% ifi=1/100;waitframes=2;sr=2000;
% determine a threshold of numbers of frames in the button press interval
threshold = ifi*waitframes*srL*1.5; % number of samples within a interval larger than the waitframes
% button presses recorded via DP port into stimtracker is sampled at 2000 Hz
% feedback via photocell to stimtracker has minial interval every other waitframes (minmal of 1000 Hz)
FeedbTimeIndL=FeedbIndL(find([0 diff(FeedbIndL')]>threshold)); % exact index of key press onset in datatimes (reduce several indices into one)
% create a time series that assign botton presses as 1, all other as 0
FeedbTimeL01=zeros(size(timeL));
FeedbTimeL01(FeedbTimeIndL)=1;
plot(FeedbTimeL01,'b.')
numFbL=sum(FeedbTimeL01) % num of feedbacks from the other plus 4+2x5 resting photocells

% Right Player
FeedbIndR=unique([find(TriggersR == 251);find(TriggersR == 235)]); % extract Index of real key presses in the values
plot(FeedbIndR,ones(1,length(FeedbIndR)),'ro'); % look at the above Index (one press produced several indices)
% ifi=1/100;waitframes=2;sr=2000;
% determine a threshold of numbers of frames in the button press interval
threshold = ifi*waitframes*srR*1.5; % number of samples within a interval larger than the waitframes
% button presses recorded via DP port into stimtracker is sampled at 2000 Hz
% feedback via photocell to stimtracker has minial interval every other waitframes (minmal of 1000 Hz)
FeedbTimeIndR=FeedbIndR(find([0 diff(FeedbIndR')]>threshold)); % exact index of key press onset in datatimes (reduce several indices into one)
% create a time series that assign botton presses as 1, all other as 0
FeedbTimeR01=zeros(size(timeL));
FeedbTimeR01(FeedbTimeIndR)=1;
plot(FeedbTimeR01,'r.')
numFbR=sum(FeedbTimeR01) % num of feedbacks from the other plus 4+2x5 resting photocells

%% extract time points for pacers (using this as index to segment data into conditions)
% Left player
% PacerIndL=unique([find(TriggersL == 127)]); % extract Index of real key presses in the values
% PacerIndL=unique([find(TriggersL == 127);find(TriggersL == 111)); % extract Index of real key presses in the values
PacerIndL=unique([find(TriggersL == 127);find(TriggersL == 111);find(TriggersL == 123)]); % extract Index of real key presses in the values
plot(PacerIndL,ones(1,length(PacerIndL)),'ro'); % look at the above Index (one press produced several indices)
% Pacer at 2Hz has intervals of 500ms (1000 samples)
% Pacer at 3Hz has interval of 333ms (666.6 samples)
% determine a threshold of numbers of samples in pacer interval
threshold = ifi*srL*0.75; % number of samples within a interval larger than the screen refresh rate
% button presses recorded via DP port into stimtracker is sampled at 2000 Hz
% feedback via photocell to stimtracker has minial interval every other waitframes (minmal of 1000 Hz)
PacerTimeIndL=PacerIndL(find([diff([0 PacerIndL'])]>threshold)); % exact index of key press onset in datatimes (reduce several indices into one)
% 2+2+30*6+2*5 = 194 % 2+30*6
% create a time series that assign botton presses as 1, all other as 0
PacerTimeL01=zeros(size(timeL));
PacerTimeL01(PacerTimeIndL)=1;
plot(PacerTimeL01,'r');
sum(PacerTimeL01)
% xlim([-1 length(PacerTimeL01)]);

% Right player
PacerIndR=unique([find(TriggersR == 127);find(TriggersR == 111)]); % extract Index of real key presses in the values
plot(PacerIndR,ones(1,length(PacerIndR)),'bo'); % look at the above Index (one press produced several indices)
% Pacer at 2Hz has intervals of 500ms (1000 samples)
% Pacer at 3Hz has interval of 333ms (666.6 samples)
% determine a threshold of numbers of samples in pacer interval
threshold = ifi*srR*0.75; % number of samples within a interval larger than the screen refresh rate
% button presses recorded via DP port into stimtracker is sampled at 2000 Hz
% feedback via photocell to stimtracker has minial interval every other waitframes (minmal of 1000 Hz)
PacerTimeIndR=PacerIndR(find([diff([0 PacerIndR'])]>threshold)); % exact index of key press onset in datatimes (reduce several indices into one)
% 2+2+30*6+2*5 = 194
% create a time series that assign botton presses as 1, all other as 0
PacerTimeR01=zeros(size(timeR));
PacerTimeR01(PacerTimeIndR)=1;
plot(PacerTimeR01,'b');
sum(PacerTimeR01)

% Plot the not aligned time points
figure;
PacersL=PacerTimeIndL([1 2 3 4 5 34 35 36 37 66 67 68 69 98 99 100 101 130 131 132 133 162 163 164 165 194]);
plot(PacersL, ones(length(PacersL)),'r.');
hold on;
PacersR=PacerTimeIndR([1 2 3 4 5 34 35 36 37 66 67 68 69 98 99 100 101 130 131 132 133 162 163 164 165 194]);
plot(PacersR, ones(length(PacersR)),'b.');

% make shift for missming pacers in one player if that happened
ShiftPacer=PacerTimeIndL(end)-PacerTimeIndR(end);
PacersL=PacersR+ShiftPacer;
figure;
plot(PacersL,ones(1,length(PacersL)),'r.');hold on;ylim([0 3]);
plot(PacerTimeIndL,2*ones(1,length(PacerTimeIndL)),'r.');

% Plot the aligned time points
figure;
plot(PacersL-PacersL(1), ones(length(PacersL)),'r.');
hold on;
plot(PacersR-PacersR(1), ones(length(PacersR)),'b.');

