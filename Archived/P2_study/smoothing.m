function [smoothx] = smoothing(x,win)
% x could be a time series of errors or button presses intervals
% win is the number of samples in the smoothing window
smoothx=[];
for i=1:length(x)-win
    smoothx(i)=mean(x(i:i+win-1));
end
end

