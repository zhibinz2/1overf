cd /home/zhibin/Downloads/HNL/FFT 2019/Demo code/code20191216.m
cd /home/zhibin/Downloads/HNL/FFT 2019/IN-CLASS+FFT+LAB
open IN-CLASS+FFT+LAB.pdf

Fs=1000
f=5
theta = 0; %phase in degrees
T = 1; %duration of 2 secs
A = 1; %amplitude of 1
N=Fs*T;
dt = 1/Fs; % time steps
time = [dt:dt:T]; %time at which each sample is taken
theta = theta*pi/180; %convert theta to radians
sinewave = A*cos(2*pi*f*time+theta); %build the sine wave
figure;plot(time,sinewave,'b');

fcoef = fft(sinewave);
fcoef=fcoef/N;
amplitude = abs(fcoef);
plot(amplitude)
amplitude = 2*amplitude;

% spectra
df=1/T;
halfN=floor(N/2);
fV=[0:df:(halfN-1)*df]; 
fV=linspace(0,Fs/2,halfN+1); % same df
fcoef=2*fcoef(1:halfN)
amplitude = abs(fcoef);
plot(fV(1:length(amplitude)),amplitude);
xlim([0 25]);
xlabel('frequency');
ylabel('amplitude (uV)');
title('Spetrum');

%% theta =0
clc; clear;
close all;

Fs=1000
f=5
theta = 0; %phase in degrees
T = 2; %duration of 2 secs
A = 2; %amplitude of 1
dt = 1/Fs; % time steps
time = [dt:dt:T]; %time at which each sample is taken
theta = theta*pi/180; %convert theta to radians
sinewave = A*cos(2*pi*f*time+theta); %build the sine wave
plot(time,sinewave,'b')

clf
fcoef = fft(sinewave);
amplitude = abs(fcoef);
plot(amplitude)

angle(fcoef(11))
angle(fcoef(11))*180/pi   % 1.8

%% theta = 45
clc; clear;
close all;

Fs=1000
f=5
theta = 45; %phase in degrees
T = 2; %duration of 2 secs
A = 2; %amplitude of 1
dt = 1/Fs; % time steps
time = [dt:dt:T]; %time at which each sample is taken
theta = theta*pi/180; %convert theta to radians
sinewave = A*cos(2*pi*f*time+theta); %build the sine wave
plot(time,sinewave,'b')

clf
fcoef = fft(sinewave);
amplitude = abs(fcoef);
plot(amplitude)
amplitude = abs(fcoef)/(T*Fs);
plot(amplitude)
angle(fcoef(11))*180/pi  %46.8


%% Fs increase to 10000
Fs=10000
f=5
theta = 45; %phase in degrees
T = 2; %duration of 2 secs
A = 2; %amplitude of 1
dt = 1/Fs; % time steps
time = [dt:dt:T]; %time at which each sample is taken
theta = theta*pi/180; %convert theta to radians
sinewave = A*cos(2*pi*f*time+theta); %build the sine wave
plot(time,sinewave,'b')

clf
fcoef = fft(sinewave);
amplitude = abs(fcoef)/(T*Fs);
plot(amplitude);

angle(fcoef(11))*180/pi  %45.18


%%
Fs=1000
f=5
theta = 45; %phase in degrees
T = 2; %duration of 2 secs
A = 2; %amplitude of 1
dt = 1/Fs; % time steps
time = [dt:dt:T]; %time at which each sample is taken
theta = theta*pi/180; %convert theta to radians
sinewave = A*cos(2*pi*f*time+theta); %build the sine wave
plot(time,sinewave,'b')

clf
fcoef = fft(sinewave);
amplitude = abs(fcoef)/(T*Fs);
plot(amplitude)

%% 

Fs=1000
f=5
theta = 45; %phase in degrees
T = 2; %duration of 2 secs
A = 2; %amplitude of 1
dt = 1/Fs; % time steps
time = [dt:dt:T]; %time at which each sample is taken
theta = theta*pi/180; %convert theta to radians
sinewave = A*cos(2*pi*f*time+theta); %build the sine wave
plot(time,sinewave,'b')

clf
fcoef = fft(sinewave);
amplitude = abs(fcoef)/(T*Fs);
plot(amplitude)

fcoef(11)

%% theta to 90
Fs=1000
f=5
theta = 90; %phase in degrees
T = 2; %duration of 2 secs
A = 2; %amplitude of 1
dt = 1/Fs; % time steps
time = [dt:dt:T]; %time at which each sample is taken
theta = theta*pi/180; %convert theta to radians
sinewave = A*cos(2*pi*f*time+theta); %build the sine wave
plot(time,sinewave,'b')


fcoef = fft(sinewave);
amplitude = abs(fcoef)/(T*Fs);
plot(amplitude)
angle(fcoef(11))*180/pi  %91.8




hold on 
Fs=1000
f=5
theta = 45; %phase in degrees
T = 2; %duration of 2 secs
A = 2; %amplitude of 1
dt = 1/Fs; % time steps
time = [dt:dt:T]; %time at which each sample is taken
theta = theta*pi/180; %convert theta to radians
sinewave = A*cos(2*pi*f*time+theta); %build the sine wave
plot(time,sinewave,'g')



