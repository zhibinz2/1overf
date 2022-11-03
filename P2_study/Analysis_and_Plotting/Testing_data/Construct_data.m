seeds=[20220713;20220721;20220804;20220808;20220810;20220811;20220815;20220816;20221003;2022100401;
        2022100402;20221005];
numSes=size(seeds,1);

% organize conditions in all sessions into a vector of 288
condition_all=[]; % 12 session x 12 trials
% initialize the testing dataset
testing_data=[]; % 12 session x 12 trials
tic
for s=1:numSes
    clear conditions
    runid=num2str(seeds(s,:));
    load(['/ssd/zhibin/1overf/' runid '_2P/Cleaned_data/clean_' runid '.mat' ]);
    condition_all(s,:)=conditions; 
    if mod(s,2)==0;
        % assigned negative to synco condi
        testing_data(s,:)=-1*conditions-s/100; 
    else mod(s,2)==1;
        % append sequential integer as two decimals for each pair of subject in each session
        testing_data(s,:)=conditions+s/100; 
    end
end
toc % took 78 sec

% reshape into a vector in time sequence
condition_all=reshape(condition_all',[],1); % 144 x 1 
testing_data=reshape(testing_data',[],1); % 144 x 1 

