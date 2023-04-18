%% compute pow in all sessions 
% original order of the experiment sessions
seeds=[20220713;20220721;20220804;20220808;20220810;20220811;20220815;20220816;20221003;2022100401;
        2022100402;20221005];
numSes=size(seeds,1); 
% get the sequence of conditions in all sessions
run get_condition_all.m 
% get the pow matrix in all essions
run get_pow_all.m 

%% compute average pow for each condition
% name of the 4 states
states4names={'Uncoupled','Leading','Following','Mutual'};
% get the averge pow for each state
run get_state_pow.m

%% Plot the average pow for the 4 states
% set plot colors
run plot_settings.m
% draw the plot
run draw_plots.m
