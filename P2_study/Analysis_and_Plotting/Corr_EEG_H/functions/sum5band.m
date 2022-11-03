function [delta, theta, alpha, beta, gamma] = sum5band(EEG,sr,win,varargin)
%This funciton get the valuses of power sum from the 5 EEG band width

% fft to get power sum of delta, theta, alpha, beta, gamma 
% define delta, theta, alpha, beta, gamma band width as:
% https://www.sinhaclinic.com/what-are-brainwaves/
% delta 1-3, theta 4-7, alpha 8-12, beta 13-38, gamma 39-42
% https://nhahealth.com/brainwaves-the-language/
% delta 1-3, theta 4-7, alpha 8-12, beta 13-30, gamma 30-44
% my band class
% delta 1-3, theta 4-7, alpha 8-12, beta 13-30, gamma 31-50
%%
% Input: EEG in time x chan
% Output: xx is EEG in time x chan x chunks

%%
if nargin < 2 | isempty(sr)
  sr = 2000;
end
if nargin < 3 | isempty(win)
  win = 0.5;
end

SamplesLength=sr*win; % number of samples in the window
nTrials=floor(size(EEG,1)/SamplesLength);

% reorganize the time series into trials (time x chunks)
xx=[];
for i=1:nTrials
    xx(:,:,nTrials)=EEG((SamplesLength*(i-1)+1):(SamplesLength*i),:);
end


%%

df=1/win;
maxfreq=50; % Hz
nbins = ceil(maxfreq/df) + 1;
freqs=(0:1:nbins-1)*df;
fcoef = fft(xx,[],1);
% pow= squeeze(var(fcoef(1:nbins,:,:),[],3)); % why all zeros?
eppow = abs(squeeze(mean(fcoef(1:nbins,:,:),3)));
% my band class
% delta 1-3, theta 4-7, alpha 8-12, beta 13-30, gamma 30-50
deltainds=find(freqs>0 & freqs <3.1);
thetainds=find(freqs>3 & freqs <7.9);
alphainds=find(freqs>7.9 & freqs <12.1);
betainds=find(freqs>12.1 & freqs <30.1);
gammainds=find(freqs>30.1 & freqs <50.1);
delta = sum(eppow(deltainds,:),1); % sum_of_pow x chan
theta = sum(eppow(thetainds,:),1); 
alpha = sum(eppow(alphainds,:),1); 
beta  = sum(eppow(betainds,:),1); 
gamma = sum(eppow(gammainds,:),1); 
end