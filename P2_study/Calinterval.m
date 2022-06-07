function [pressIntervals] = Calinterval(A)
% A is a time series of 0 and 1 (should be in a row)
% pressIntervals is the number of samples (zeros) between each button press (one)
% A could be the button press or the feedback
pressInd=find([0 diff(A)]==1);
pressIntervals=[pressInd(1) diff(pressInd)];
end

