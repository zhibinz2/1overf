figure;
tic
cmin=0;cmax=1e-4;
for st=1:4
    subplot(1,4,st);
    imagesc(squeeze(avg_state_pow(st,:,:)));colorbar;colormap('jet');
    clim([cmin cmax]);
    xlabel('chans');ylabel('freq');set(gca, 'YDir','normal');
    title(states4names{st},'color',condicolors(st,:));
end
toc
