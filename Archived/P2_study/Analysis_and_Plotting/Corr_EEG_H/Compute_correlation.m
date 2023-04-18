% Organize the data needed for the analysis
addpath functions/
run Organize_EEG.m
run Organize_H.m
run Organize_indices.m

% Compute correlation in 4 states
delta_LR_H_LR_4corr=zeros(4,32);
theta_LR_H_LR_4corr=zeros(4,32);
alpha_LR_H_LR_4corr=zeros(4,32);
beta_LR_H_LR_4corr=zeros(4,32);
gamma_LR_H_LR_4corr=zeros(4,32);
for s=1:4
    for c=1:32
        delta_LR_H_LR_4corr(s,c)=corr(delta_LR_chan(Inds4_LR(:,s),c),H_all_LR(Inds4_LR(:,s)));
        theta_LR_H_LR_4corr(s,c)=corr(theta_LR_chan(Inds4_LR(:,s),c),H_all_LR(Inds4_LR(:,s)));
        alpha_LR_H_LR_4corr(s,c)=corr(alpha_LR_chan(Inds4_LR(:,s),c),H_all_LR(Inds4_LR(:,s)));
        beta_LR_H_LR_4corr(s,c)=corr(beta_LR_chan(Inds4_LR(:,s),c),H_all_LR(Inds4_LR(:,s)));
        gamma_LR_H_LR_4corr(s,c)=corr(gamma_LR_chan(Inds4_LR(:,s),c),H_all_LR(Inds4_LR(:,s)));
    end
end

% Design color scheme for the plotting
red   = [1 0 0];
pink  = [1 0.65 0.75];
blue  = [0 0 1];
mediumblue = [0 0.4 0.7];
green = [0 1 0];
darkgreen = [0 0.5 0];
grey  = [0.5 0.5 0.5];
yellow  = [1 1 0];
deepyellow  = [1 0.8 0.2];
megenta = [1 0 1];% fill([0 1 1 0],[0 0 1 1],megenta)
cyan = [0 1 1]; % fill([0 1 1 0],[0 0 1 1],cc)
purple = [0.6 0.1 0.9];
condicolors=[darkgreen;red;blue;megenta;purple;purple];
HNLcolors = [darkgreen; deepyellow; pink];
addpath /home/zhibin/Documents/GitHub/matlab-archive/hnlcode/common/gen_code/color
hnc = hotncold(100);

% Plot
% Combine L and R in 4 states(4x5) 
canvas(0.3,0.5);
cmin=-0.7;cmax=0.7;
for s=1:4
    subplot(4,5,5*(s-1)+1);
    topoplot(delta_LR_H_LR_4corr(s,:),channels,'nosedir','+X','style','map');
    clim([cmin cmax]);
    subplot(4,5,5*(s-1)+2);
    topoplot(theta_LR_H_LR_4corr(s,:),channels,'nosedir','+X','style','map');
    clim([cmin cmax]);
    subplot(4,5,5*(s-1)+3);
    topoplot(alpha_LR_H_LR_4corr(s,:),channels,'nosedir','+X','style','map');
    clim([cmin cmax]);
    subplot(4,5,5*(s-1)+4);
    topoplot(beta_LR_H_LR_4corr(s,:),channels,'nosedir','+X','style','map');
    clim([cmin cmax]);
    subplot(4,5,5*(s-1)+5);
    topoplot(gamma_LR_H_LR_4corr(s,:),channels,'nosedir','+X','style','map');
    clim([cmin cmax]);
end
colormap(hnc)
cb=colorbar;
cb.AxisLocation = 'out';
cb.Position = [0.92 0.15 0.01 0.75];
delete(findall(gcf,'type','annotation'))
h0=annotation('textbox',[0.17 0.95 0.05 0.03],'string','Delta','color',[0 0 0])
h1=annotation('textbox',[0.33 0.95 0.05 0.03],'string','Theta','color',[0 0 0])
h2=annotation('textbox',[0.5 0.95 0.05 0.03],'string','Alpha','color',[0 0 0])
h3=annotation('textbox',[0.66 0.95 0.05 0.03],'string','Beta','color',[0 0 0])
h4=annotation('textbox',[0.81 0.95 0.05 0.03],'string','Gamma','color',[0 0 0])
v0=annotation('textbox',[0.14 0.15 0.05 0.03],'string','Mutual','color',condicolors(4,:))
v1=annotation('textbox',[0.14 0.37 0.05 0.03],'string','Following','color',condicolors(3,:))
v2=annotation('textbox',[0.14 0.59 0.05 0.03],'string','Leading','color',condicolors(2,:))
v3=annotation('textbox',[0.14 0.81 0.05 0.03],'string','Uncouple','color',condicolors(1,:))
set(v0,'Rotation',90);set(v1,'Rotation',90);set(v2,'Rotation',90);set(v3,'Rotation',90);
set([h0 h1 h2 h3 h4 v0 v1 v2 v3], 'fitboxtotext','on',...
    'edgecolor','none')
sg=annotation('textbox',[0.3 0.01 0.4 0.07],'string',...
    {['Correlation of sum-EEG (-500ms) and H-int ^{* PLOT 13}'],char(datetime('now'))});
set(gcf,'color','w'); % set background white for copying in ubuntu
