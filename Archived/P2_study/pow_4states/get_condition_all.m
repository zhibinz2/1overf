% get the sequence of conditions
% organize conditions from all sessions into a matrix of (12 sessions x 12 trials)
condition_all=[];
tic
for s=1:numSes
    clear conditions
    runid=num2str(seeds(s,:));
    load(['/ssd/zhibin/1overf/' runid '_2P/Cleaned_data/clean_' runid '.mat'],'conditions');
    condition_all(s,:)=conditions; % 12 session x 12 trials
end
toc % 0.019422 seconds.
% reshape into a vector in time sequence
% condition_all=reshape(condition_all',[],1); % 144 x 1 