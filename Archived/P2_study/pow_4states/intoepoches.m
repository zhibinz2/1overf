function [epoches] = intoepoches(data,sr,win)
% data in timexchans
% this function organize data into timexchansxepoches

if win==[]; win=1; end % chunck of 1 seconds
SamplesLength=sr*win; % number of samples in the window

% maximum number of 1 second windows (# of chuncks)
nTrials=floor(size(data,1)/SamplesLength);

% zcore transform
zdata=zscore(data);

% reorganize the time series into trials (time x chans x chunks)
epoches=[];
for i=1:nTrials
    epoches(:,:,nTrials)=zdata((SamplesLength*(i-1)+1):(SamplesLength*i),:);
end

end