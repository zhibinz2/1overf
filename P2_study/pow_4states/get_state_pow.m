% reorganzie pow_all according to the 4 states
condi_pow_all=cell(1,4); % 4 state {51 freq x 32 chan x 72 trials}
tic
for s=1:numSes
        for t=1:12 % 12 trials
            if condition_all(s,t)==1; % Uncouple
                condi_pow_all{1}=cat(3,condi_pow_all{1}, squeeze(pow_all(s,1,t,:,:)));
                condi_pow_all{1}=cat(3,condi_pow_all{1}, squeeze(pow_all(s,2,t,:,:)));
            elseif condition_all(s,t)==2; % L-lead
                condi_pow_all{2}=cat(3,condi_pow_all{2}, squeeze(pow_all(s,1,t,:,:)));% L leading
                condi_pow_all{3}=cat(3,condi_pow_all{3}, squeeze(pow_all(s,2,t,:,:)));% R following
            elseif condition_all(s,t)==3; % R-lead
                condi_pow_all{2}=cat(3,condi_pow_all{2}, squeeze(pow_all(s,2,t,:,:)));% R leading
                condi_pow_all{3}=cat(3,condi_pow_all{3}, squeeze(pow_all(s,1,t,:,:)));% L following
            else
                condi_pow_all{4}=cat(3,condi_pow_all{4}, squeeze(pow_all(s,1,t,:,:)));
                condi_pow_all{4}=cat(3,condi_pow_all{4}, squeeze(pow_all(s,2,t,:,:)));
            end
        end
end
toc % Elapsed time is 0.055227 seconds.
% get average pow for each state
avg_state_pow=nan(4,51,32); % 4 state x 51 freq x 32 chan
for st=1:4
    avg_state_pow(st,:,:)=mean(condi_pow_all{st},3);
end
