%% Viwe raw EEG
timeL;samplesL;TRIGGERindL;srL;channels_infoL;
timeR;samplesR;TRIGGERindR;srR;channels_infoR;

figure;
subplot(2,1,1);plot(samplesL(1:32,:)');
subplot(2,1,2);plot(samplesR(1:32,:)');

figure; % look at the bad chan P8
subplot(2,1,1);plot(samplesL(28,:)');
subplot(2,1,2);plot(samplesR(28,:)');

%% Extract EEG
EEGL=samplesL(1:32,:)';
EEGR=samplesR(1:32,:)';

%% detrend the EEG data (no padding needed)
detrend_dataL=detrend(EEGL,2);
detrend_dataR=detrend(EEGR,2);

figure('units','normalized','outerposition',[0 0 1 0.6]);
subplot(2,1,1);plot(detrend_dataL); ylim([-1000 1000]);
subplot(2,1,2);plot(detrend_dataR); ylim([-1000 1000]);

%% Filter the data
% filtfilthd method (hnl) high pass first then low pass filter
cd /usr/local/MATLAB/R2019a/toolbox/signal/signal
which filtfilt
which filtfilt -all
cd D:\Program Files\MATLAB\R2019b\toolbox\signal\signal\filtfilt.m   
% high pass (no paddings needed)
Hd = makefilter(srL,0.2,0.15,6,20,0); 
filtered_dataL1=filtfilthd(Hd,detrend_dataL);
filtered_dataR1=filtfilthd(Hd,detrend_dataR);

figure('units','normalized','outerposition',[0 0 1 0.6]);
subplot(2,1,1);plot(filtered_dataL1); ylim([-300 300]);
subplot(2,1,2);plot(filtered_dataR1); ylim([-300 300]);

% add padding (this step takes several minutes)
paddingL=zeros(round(size(filtered_dataL1,1)/10), size(filtered_dataL1,2));
filtered_dataL2=cat(1,paddingL,filtered_dataL1,paddingL);
paddingR=zeros(round(size(filtered_dataR1,1)/10), size(filtered_dataR1,2));
filtered_dataR2=cat(1,paddingR,filtered_dataR1,paddingR);
figure('units','normalized','outerposition',[0 0 1 0.6]);
subplot(2,1,1);;plot(filtered_dataL2);ylim([-300 300]);title('filtered-dataL2');
subplot(2,1,2);plot(filtered_dataR2);ylim([-300 300]);title('filtered-dataR2');

% low pass 
% (this will create short edge artifact)
% (if added zero paddings, edge artifact disappear)
% (remove existing filtered_data variable from workspace might fasten)
% filtered_dataL3=[];
Hd = makefilter(srL,50,51,6,20,0);  
filtered_dataL3=filtfilthd(Hd,filtered_dataL2);
filtered_dataR3=filtfilthd(Hd,filtered_dataR2);
figure('units','normalized','outerposition',[0 0 1 0.6]);
subplot(2,1,1);plot(filtered_dataL3);ylim([-300 300]);
subplot(2,1,2);plot(filtered_dataR3);ylim([-300 300]);
% plotx(filtered_dataR3(ind1:ind2,:));ylim([-300 300]); % channel 30 bad

% remove padding
filtered_dataL4=filtered_dataL3((size(paddingL,1)+1):(size(paddingL,1)+size(detrend_dataL,1)),:);
filtered_dataR4=filtered_dataR3((size(paddingR,1)+1):(size(paddingR,1)+size(detrend_dataR,1)),:);
figure('units','normalized','outerposition',[0 0 1 0.6]);
subplot(2,1,1);plotx(filtered_dataL4);ylim([-300 300]);
subplot(2,1,2);plotx(filtered_dataR4);ylim([-300 300]);
% ylim([-100 100]);

clearvars filtered_dataL1 filtered_dataL2 filtered_dataL3 paddingL
clearvars filtered_dataR1 filtered_dataR2 filtered_dataR3 paddingR

%% ICA
filtered_dataL4; filtered_dataR4;
% Left player
% run ICA
[icasigL, AL, WL] = fastica(filtered_dataL4');
% Plot ICA component
for i=1:size(AL,2)
    SqAL(i)=sumsqr(AL(:,i));
end
figure;
plot(1:size(AL,2),SqAL,'ro');ylabel('sum of square of column in A');xlabel('ICs'); title('L');
[BL,IL]=sort(SqAL,'descend');
ComponentsExam=IL(1:10);
% topoplot to examine them
cd /home/zhibin/Documents/GitHub/Motor_cordination/1_over_f/data_analysis/channels_info
load('chaninfo.mat')
figure;
for i=1:length(ComponentsExam)
    subplot(5,2,i);
    topoplot(AL(:,ComponentsExam(i)),chaninfo,'nosedir','+X');title(['component' num2str(ComponentsExam(i))]);colorbar;
end
suptitle('L');
% Calculate Correlation
% FP1 and FP2 are channel 1 and 3;
% compute correlation between FP1 and FP2;
[RHO1,PVAL1] = corr(filtered_dataL4(:,1),icasigL');
[RHO3,PVAL3] = corr(filtered_dataL4(:,3),icasigL');
[B1,I1]=sort(abs(RHO1),'descend');[B3,I3]=sort(abs(RHO3),'descend');
ComponentsExamL=unique([I1(1:2) I3(1:2)])
figure;
for i=1:length(ComponentsExamL)
    subplot(length(ComponentsExamL),1,i);
    topoplot(AL(:,ComponentsExamL(i)),chaninfo,'nosedir','+X');title(['component' num2str(ComponentsExamL(i))]);colorbar;
end
suptitle('L');
% Decide component to remove
ComponentRemoveL=2; % ComponentRemoveL=ComponentsExamL;

% Right player
% run ICA
[icasigR, AR, WR] = fastica(filtered_dataR4');
% Plot ICA component
for i=1:size(AR,2)
    SqAR(i)=sumsqr(AR(:,i));
end
figure;
plot(1:size(AR,2),SqAR,'ro');ylabel('sum of square of column in A');xlabel('ICs');title('R');
[BR,IR]=sort(SqAR,'descend');
ComponentsExam=IR(1:10);
figure;
for i=1:length(ComponentsExam)
    subplot(5,2,i);
    topoplot(AR(:,ComponentsExam(i)),chaninfo,'nosedir','+X');title(['component' num2str(ComponentsExam(i))]);colorbar;
end
suptitle('R');
% Calculate Correlation
% FP1 and FP2 are channel 1 and 3;
% compute correlation between FP1 and FP2;
[RHO1,PVAL1] = corr(filtered_dataR4(:,1),icasigR');
[RHO3,PVAL3] = corr(filtered_dataR4(:,3),icasigR');
[B1,I1]=sort(abs(RHO1),'descend');[B3,I3]=sort(abs(RHO3),'descend');
ComponentsExamR=unique([I1(1:2) I3(1:2)])
% topoplot to examine them
figure;
for i=1:length(ComponentsExamR)
    subplot(length(ComponentsExamR),1,i);
    topoplot(AR(:,ComponentsExamR(i)),chaninfo,'nosedir','+X');title(['component' num2str(ComponentsExamR(i))]);colorbar;
end
% Decide component to remove
ComponentRemoveR=[27 29]; % ComponentRemoveL=ComponentsExamL;

%% Remove ICA components, mix back the signal and display for comparison
% Left player
ALrm=AL;icasigLrm=icasigL; % make new, backup AL, icasigL
ALrm(:,ComponentRemoveL)=0; icasigLrm(ComponentRemoveL,:)=0;
mixedsigL=ALrm*icasigLrm;
mixedsigL=mixedsigL';
% Plot before and after
figure('units','normalized','outerposition',[0 0 1 1]);
% before ICA
subplot(2,1,1);
plotx(timeL,filtered_dataL4);
hold on;hold off; ylim([-200 200]);
title('EEG Singal Before ICA');
% after ICA
subplot(2,1,2);
plotx(timeL,mixedsigL);
hold on;hold off; ylim([-200 200]);
title('Mixed Signal with ICs removed');
% selet two point on x axis to zoom in
ind1=round(length(timeL)*0.55);
ind2=round(length(timeL)*0.6);
figure('units','normalized','outerposition',[0 0 1 1]);
% before ICA
subplot(2,1,1);
plotx(timeL(ind1:ind2),filtered_dataL4(ind1:ind2,:));
hold on;hold off; ylim([-200 200]);
title('EEG Singal Before ICA');
% after ICA
subplot(2,1,2);
plotx(timeL(ind1:ind2),mixedsigL(ind1:ind2,:));
hold on;hold off; ylim([-200 200]);
title('Mixed Signal with ICs removed');

% Right player
ARrm=AR;icasigRrm=icasigR; % make new, backup AR, icasigR
ARrm(:,ComponentRemoveR)=0; icasigRrm(ComponentRemoveR,:)=0;
mixedsigR=ARrm*icasigRrm;
mixedsigR=mixedsigR';
% Plot before and after
figure('units','normalized','outerposition',[0 0 1 1]);
% before ICA
subplot(2,1,1);
plotx(timeR,filtered_dataR4);
hold on;hold off;ylim([-100 100]);
title('EEG Singal Before ICA');
% after ICA
subplot(2,1,2);
plotx(timeR,mixedsigR);
hold on;hold off;ylim([-100 100]);
title('Mixed Signal with ICs removed');
% selet two point on x axis to zoom in
ind1=round(length(timeL)*0.55);
ind2=round(length(timeL)*0.6);
% zoom in and plot again
figure('units','normalized','outerposition',[0 0 1 1]);
% before ICA
subplot(2,1,1);
plotx(timeR(ind1:ind2),filtered_dataR4(ind1:ind2,:));
% plotx(timeR(ind1:ind2),filtered_dataR4(ind1:ind2,1:29)); 
hold on;hold off;ylim([-100 100]);
title('EEG Singal Before ICA');
% after ICA
subplot(2,1,2);
plotx(timeR(ind1:ind2),mixedsigR(ind1:ind2,:));
% plotx(timeR(ind1:ind2),mixedsigR(ind1:ind2,1:29));
hold on;hold off;ylim([-100 100]);
title('Mixed Signal with ICs removed');

%% examine spectra on scalp map 
run /home/zhibin/Documents/GitHub/Motor_cordination/1_over_f/data_analysis/channels_info.m
Timeselected=timeL(ind1:ind2);% Timeselected=timeR(ind1:ind2); 
EEGselected=filtered_dataL4(ind1:ind2,:);% EEGselected=filtered_dataR4(ind1:ind2,:);
% EEGselected=mixedsigL(ind1:ind2,:);% EEGselected=mixedsigR(ind1:ind2,:);
% fft
fcoef=fft(EEGselected);
N=size(EEGselected,1); % N=length(Timeselected);
fcoef=fcoef/N;
halfN=floor(N/2);
% df=1/T; 
% fV=[0:df:(halfN-1)*df]; 
fV=linspace(0,srR/2,halfN+1); % same df
fcoef=2*fcoef(1:halfN,:);
amplitude = abs(fcoef);
% Plot on scalp map for the spectra (for my own examing)
figure('units','normalized','outerposition',[0 0 1 1]);
for chan=1:32
    subplot('Position',[XXPLOT(chan) YYPLOT(chan) 0.05 0.05]); % not showing, why
    % plot([1 1],[1 1],'ro');
    plot(fV(1:size(amplitude,1)),amplitude(:,chan));
    if ~isempty(find([1:30 32]==chan))
    set(gca,'XTick',[]); 
    % set(gca,'YTick',[]); 
    end
    if chan==31
        xlabel('frequency');
        ylabel('amplitude (uV)');
    end
    xlim([0 25]);title([labels{chan}]);%ylim([0 1]);
end
suptitle('spectra of all channels on scalp map')