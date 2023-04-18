% Compute sum of EEG power
seeds=[20220713;20220721;20220804;20220808;20220810;20220811;20220815;20220816;20221003;2022100401;
        2022100402;20221005];
numSes=size(seeds,1);
EEGwin=0.5; % second

% EEG -500ms before the matched tap
zEEG500_L=cell(numSes,12);zEEG500_R=cell(numSes,12);
delta_L=cell(numSes,12);theta_L=cell(numSes,12);alpha_L=cell(numSes,12);beta_L=cell(numSes,12);gamma_L=cell(numSes,12);
delta_R=cell(numSes,12);theta_R=cell(numSes,12);alpha_R=cell(numSes,12);beta_R=cell(numSes,12);gamma_R=cell(numSes,12);
cd /ssd/zhibin/1overf/
tic
for s=1:numSes % each session
    runid=num2str(seeds(s,:));
    clear dataL dataR
    load(['/ssd/zhibin/1overf/' runid '_2P/Cleaned_data/clean_' runid '.mat' ])
    zEEG_L={};zEEG_R={};
    for b=1:12 % each block  
        zEEG_L{b}=zscore(dataL{b}(:,1:32),[],1);
        zEEG_R{b}=zscore(dataR{b}(:,1:32),[],1);
        for i=1:size(samples{b},1) % each matched interval           
            zEEG500_L{s,b}{i}=zEEG_L{b}(samples{b}(i,1)-999:samples{b}(i,1),:); % start from first matched tap (if not enough samples before, then consider starting from 2nd tap)
            zEEG500_R{s,b}{i}=zEEG_R{b}(samples{b}(i,2)-999:samples{b}(i,2),:);
            [delta_L{s,b}(i,:), theta_L{s,b}(i,:), alpha_L{s,b}(i,:), beta_L{s,b}(i,:), gamma_L{s,b}(i,:)]...
                =sum5band(zEEG500_L{s,b}{i},sr,EEGwin);
            [delta_R{s,b}(i,:), theta_R{s,b}(i,:), alpha_R{s,b}(i,:), beta_R{s,b}(i,:), gamma_R{s,b}(i,:)]...
                =sum5band(zEEG500_R{s,b}{i},sr,EEGwin);
        end
    end
end
toc % about 2 min

% organize EEG power
zEEG500_L; zEEG500_R; % for all sessions from SECT 12
delta_L;theta_L;alpha_L;beta_L;gamma_L;
delta_R;theta_R;alpha_R;beta_R;gamma_R;

% save these EEG power from each tap
save('20220713_1005_EEG_5band_pow_taps.mat','delta_L','theta_L','alpha_L','beta_L','gamma_L',...
    'delta_R','theta_R','alpha_R','beta_R','gamma_R');

% sum the band power in each channel in each block (12session x 12blocks = 144)
delta_L_sum=[];theta_L_sum=[];alpha_L_sum=[];beta_L_sum=[];gamma_L_sum=[];
delta_R_sum=[];theta_R_sum=[];alpha_R_sum=[];beta_R_sum=[];gamma_R_sum=[];
for s=1:numSes
    for b=1:12
        delta_L_sum(s,b,:)=sum([delta_L{s,b}]);
        theta_L_sum(s,b,:)=sum([theta_L{s,b}]);
        alpha_L_sum(s,b,:)=sum([alpha_L{s,b}]);
        beta_L_sum(s,b,:)=sum([beta_L{s,b}]);
        gamma_L_sum(s,b,:)=sum([gamma_L{s,b}]);
        delta_R_sum(s,b,:)=sum([delta_R{s,b}]);
        theta_R_sum(s,b,:)=sum([theta_R{s,b}]);
        alpha_R_sum(s,b,:)=sum([alpha_R{s,b}]);
        beta_R_sum(s,b,:)=sum([beta_R{s,b}]);
        gamma_R_sum(s,b,:)=sum([gamma_R{s,b}]);
    end
end

% squeeze into 32 vectors for corr with H
delta_L_chan=[];theta_L_chan=[];alpha_L_chan=[];beta_L_chan=[];gamma_L_chan=[];
delta_R_chan=[];theta_R_chan=[];alpha_R_chan=[];beta_R_chan=[];gamma_R_chan=[];
for c=1:32
    delta_L_chan(:,c)=reshape(delta_L_sum(:,:,c)',[],1); % 144x32 (each element from one block in time sequence) 
    theta_L_chan(:,c)=reshape(theta_L_sum(:,:,c)',[],1);
    alpha_L_chan(:,c)=reshape(alpha_L_sum(:,:,c)',[],1);
    beta_L_chan(:,c)=reshape(beta_L_sum(:,:,c)',[],1);
    gamma_L_chan(:,c)=reshape(gamma_L_sum(:,:,c)',[],1);
    delta_R_chan(:,c)=reshape(delta_R_sum(:,:,c)',[],1);
    theta_R_chan(:,c)=reshape(theta_R_sum(:,:,c)',[],1);
    alpha_R_chan(:,c)=reshape(alpha_R_sum(:,:,c)',[],1);
    beta_R_chan(:,c)=reshape(beta_R_sum(:,:,c)',[],1);
    gamma_R_chan(:,c)=reshape(gamma_R_sum(:,:,c)',[],1);
end

% Combine L and R for 288 predictors in PLS
delta_LR_chan=[delta_L_chan;delta_R_chan]; % (288 x 32)
theta_LR_chan=[theta_L_chan;theta_R_chan];
alpha_LR_chan=[alpha_L_chan;alpha_R_chan];
beta_LR_chan=[beta_L_chan;beta_R_chan];
gamma_LR_chan=[gamma_L_chan;gamma_R_chan];

% fix the scale in the data 
delta_LR_chan = delta_LR_chan./(ones(288,1)*std(delta_LR_chan)); 
theta_LR_chan = theta_LR_chan./(ones(288,1)*std(theta_LR_chan));
alpha_LR_chan = alpha_LR_chan./(ones(288,1)*std(alpha_LR_chan));
beta_LR_chan = beta_LR_chan./(ones(288,1)*std(beta_LR_chan));
gamma_LR_chan = gamma_LR_chan./(ones(288,1)*std(gamma_LR_chan));
