
% Correlation in 4 states
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
% Combine L and R in 4 states(4x5)
canvas(0.3,0.5);
cmin=-0.7;cmax=0.7;
for s=1:4
    subplot(4,5,5*(s-1)+1);
    topoplot(delta_LR_H_LR_4corr(s,:),channels,'nosedir','+X','style','map');
    % title([states4names{s} ': delta & H'],'Color',condicolors(s,:));
    % colorbar;colormap('jet');
    clim([cmin cmax]);
    subplot(4,5,5*(s-1)+2);
    topoplot(theta_LR_H_LR_4corr(s,:),channels,'nosedir','+X','style','map');
    % title([states4names{s} ': theta & H'],'Color',condicolors(s,:));
    % colorbar;colormap('jet');
    clim([cmin cmax]);
    subplot(4,5,5*(s-1)+3);
    topoplot(alpha_LR_H_LR_4corr(s,:),channels,'nosedir','+X','style','map');
    % title([states4names{s} ': alpha & H'],'Color',condicolors(s,:));
    % colorbar;colormap('jet');
    clim([cmin cmax]);
    subplot(4,5,5*(s-1)+4);
    topoplot(beta_LR_H_LR_4corr(s,:),channels,'nosedir','+X','style','map');
    % title([states4names{s} ': beta & H'],'Color',condicolors(s,:));
    % colorbar;colormap('jet');
    clim([cmin cmax]);
    subplot(4,5,5*(s-1)+5);
    topoplot(gamma_LR_H_LR_4corr(s,:),channels,'nosedir','+X','style','map');
    % title([states4names{s} ': gamma & H'],'Color',condicolors(s,:));
    % colorbar;colormap('jet');
    clim([cmin cmax]);
end
% sgtitle('4states: Corr of sum-EEG (-500ms) and H-int ^{* PLOT 13}')
colormap(hnc)
% https://www.mathworks.com/help/matlab/ref/matlab.graphics.illustration.colorbar-properties.html
cb=colorbar;
cb.AxisLocation = 'out';
cb.Position = [0.92 0.15 0.01 0.75];
% dim is [x position, y position of starting point, width, height]
% for i=1:5
% h(i)=annotation('textbox',[(i-1)/5+15/100 0/5+19/20 2/100 1/100],'string','A','color',[0 0 0]);
% end
% for i=1:4
% v(i)=annotation('textbox',[0/5+2/20 (i-1)/5+1/20 2/100 1/100],'string','H','color',...
%     condicolors(i,:));
% end
% set([h v], 'fitboxtotext','on',...
%     'edgecolor','none')
% clf h(1)
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
delete(sg)
sg=annotation('textbox',[0.3 0.01 0.4 0.07],'string',...
    {['Correlation of sum-EEG (-500ms) and H-int ^{* PLOT 13}'],char(datetime('now'))});
set([h0 h1 h2 h3 h4 v0 v1 v2 v3], 'fitboxtotext','on',...
    'edgecolor','none')
set(gcf,'color','w'); % set background white for copying in ubuntu
