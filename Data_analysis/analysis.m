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
altertive_names
unit_names
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
Plot a single channel.

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

%% load behaviral data
cd /home/zhibin/Documents/GitHub/1overf/stimulus_data_storage
load('20220324.mat')

%% Covert BottonPres and Photocell signals to time points
BottonPres=samples(54,:);unique(BottonPres) 
Photocell=samples(38,:);plot(Photocell)
time;

open /home/zhibin/Documents/GitHub/Motor_cordination/EEGanalysis/20211102/organize_photocells_step2.m






%% auto correlation and spectrum
open /ssd/zhibin/1overf/20220324/explore_1_over_f.m

cd /ssd/zhibin/1overf/20220324

