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
    % tesing_data has the following structure: L&R assinged same value
    % synco are negaitive and synch are positive
    % sessions marked with 0.01 to 0.12
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
testing_data_L=reshape(testing_data',[],1); % 144 x 1 
testing_data_R=testing_data_L; % 144 x 1 
testing_data_LR=[testing_data_L; testing_data_R];
% testing EEG dataset in 32channels
% same structure as testing data but mark each channel with 0.0001-0.0032
channel_mark1=nan(12,32);
for s=1:numSes
    if mod(s,2)==0;
        % assigned negative to synco condi
        channel_mark1(s,:)=-[1:32]/10000; 
    else mod(s,2)==1;
        % append sequential integer as two decimals for each pair of subject in each session
        channel_mark1(s,:)=[1:32]/10000; 
    end
end
channel_mark2=nan(144,32);
for s=1:numSes
    channel_mark2([1:12]+12*(s-1),:)=repmat(channel_mark1(s,:),12,1);
end
testing_EEG_L=repmat(testing_data_L,1,32) + channel_mark2;
testing_EEG_R=repmat(testing_data_R,1,32) + channel_mark2;
testing_EEG_LR=[testing_EEG_L;testing_EEG_R];
