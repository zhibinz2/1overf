seeds=[20220713;20220721;20220804;20220808;20220810;20220811;20220815;20220816;20221003;2022100401;
        2022100402;20221005];
numSes=size(seeds,1);

% Compute H and save them in each session's directory
cd /ssd/zhibin/1overf/
H=[];intervals_H_removed=[];
tic
for s=1:numSes
    clear intervals H intervals_H_removed
    runid=num2str(seeds(s,:));
    load(['/ssd/zhibin/1overf/' runid '_2P/Cleaned_data/clean_' runid '.mat' ])   
    for b=1:12
        [~,H(1,b)]=DFA_main(intervals{b}(:,1));
        [~,H(2,b)]=DFA_main(intervals{b}(:,2));
        % d removal
        [intervals_H_removed{b}(:,1)]=remove_d(intervals{b}(:,1),H(1,b)-0.5);
        [intervals_H_removed{b}(:,2)]=remove_d(intervals{b}(:,2),H(2,b)-0.5);
        % save H and intervals_H_removed
        eval(['save ' '/ssd/zhibin/1overf/' runid '_2P/Cleaned_data/d_removal' runid '.mat H intervals_H_removed'])
    end
end
toc


% Extract H from all sessions
cd /ssd/zhibin/1overf/
H_all=[];
for s=1:numSes
    runid=num2str(seeds(s,:));
    clear H
    % load the previously computed and saved H in each session
    load(['/ssd/zhibin/1overf/' runid '_2P/Cleaned_data/d_removal' runid '.mat' ])  
    for b=1:12
        H_all(:,s,b)=H(:,b); % 2 subject x 12 session x 12 trials
    end
end

% Organize H_all for corr
H_all; % (2xnumSesx12) for all sessions from SECT 10-1 (matched int)
H_all_L=squeeze(H_all(1,:,:));
H_all_R=squeeze(H_all(2,:,:));
% squeeze into 1 vector from the 96 blocks for each subject, for corr with pow in each chan
H_all_L=reshape(H_all_L',[],1);% 144x1 (each element from one block in time sequence) 
H_all_R=reshape(H_all_R',[],1);
% Combine L&R in a sequence to be used for computing correlation and PLS prediction
H_all_LR=[H_all_L;H_all_R]; 