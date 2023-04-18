% Set screen flash 100HZ
% 0.01s per flash, 0.02 per frame (flip every other flash / every 20 ms)
% Set screen flash 144HZ
ifi=0.0069; % 0.0069 s (refresh every 6.9ms, flip every other fresh / every 13.8ms)

% Mean stimulus interval
MeanTapInterval=0.5; % second , or 500 ms
% NumFramesInterval=MeanTapInterval/ifi*2; % let's say NumFramesInterval=25 frames
NumFramesInterval=round(MeanTapInterval/ifi);  % on average 72 frames per stimulus 

% 3 trials/conditions 5 blocks
% estimated about 2 min (240 taps) per trial/condition, total 30 min of tapping
numTaps=240;

%% Total time Condition 1:
TotalTimeCondition1=MeanTapInterval*240; % in seconds

Showframes1=[1:NumFramesInterval:NumFramesInterval*numTaps];
% total time
Showframes1(end)*ifi % in seconds

%% Total time Condition 2:
% Same as condition 1 except the instruction
Showframes2=Showframes1;

%% Total time Condition 3:
% TotalTimeCondition2=0.5*240;
% RandomIntervals = round(normrnd(NumFramesInterval,NumFramesInterval/5,[1,numTaps]));% require statistic and machine learning toolbox
% RandomIntervals = round(NumFramesInterval + (NumFramesInterval/5).*randn(1,numTaps)); % normal distribution
RandomIntervals = round(NumFramesInterval + NumFramesInterval.*(rand(1,numTaps)-0.5)); % uniform distribution
mean(RandomIntervals) % mean interval , should be about equal to NumFramesInterval, which is 72 (Screen Refresh 144Hz)
plot(RandomIntervals,'ro');
% in case some interval are negative (for normal distribution)
% if any(RandomIntervals<0)
%     RandomIntervals(find(RandomIntervals<0))=1;
% end
% RandomIntervals
Showframes3=cumsum(RandomIntervals);
plot(Showframes3,ones(1,length(Showframes3)),'ro');hold on;xline(Showframes3)
Showframes3(end)
% sum(RandomIntervals)*20/1000/60 % total time
% Showframes(end)*20/1000/60 % total time
Showframes3(end)*ifi % in seconds


%% Flip loop each trial

% select a condition

Showframes=[Showframes1;Showframes2;Showframes3];

for n=1:max(Showframes)+NumFramesInterval
    if any(Showframes(:) == n);
        Screen('flip'...
        
    end
end

