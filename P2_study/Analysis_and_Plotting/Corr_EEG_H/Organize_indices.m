seeds=[20220713;20220721;20220804;20220808;20220810;20220811;20220815;20220816;20221003;2022100401;
        2022100402;20221005];
numSes=size(seeds,1);

% organize conditions in all sessions into a vector of 288
condition_all=[];
tic
for s=1:numSes
    clear conditions
    runid=num2str(seeds(s,:));
    load(['/ssd/zhibin/1overf/' runid '_2P/Cleaned_data/clean_' runid '.mat' ]);
    condition_all(s,:)=conditions; % 12 session x 12 trials
end
toc % took 78 sec

% reshape into a vector in time sequence
condition_all=reshape(condition_all',[],1); % 144 x 1 

% indices for 4 states
states4names={'Uncouple','Leading','Following','Mutual'};

% find the indices for each state in the L&R subject conbined sequence (4 states)
% organize the indicies for PLS and correlation
uncoupleInd_LR=[find(condition_all==1);12*numSes+find(condition_all==1)];
leadingInd_LR=[find(condition_all==2);12*numSes+find(condition_all==3)];
followingInd_LR=[find(condition_all==3);12*numSes+find(condition_all==2)];
mutualInd_LR=[find(condition_all==4);12*numSes+find(condition_all==4)];
Inds4_LR=[uncoupleInd_LR leadingInd_LR followingInd_LR mutualInd_LR];