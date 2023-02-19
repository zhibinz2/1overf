%% Load data in mat file
load('example_data.mat');

%% detrend the EEG data 
detrend_dataL=detrend(EEGL,2);
detrend_dataR=detrend(EEGR,2);

%% Filter the data
% high pass 
Hd = makefilter(srL,1.5,1,6,20,0); 
filtered_dataL1=filtfilthd(Hd,detrend_dataL);
filtered_dataR1=filtfilthd(Hd,detrend_dataR);
% add padding 
paddingL=zeros(round(size(filtered_dataL1,1)/10), size(filtered_dataL1,2));
filtered_dataL2=cat(1,paddingL,filtered_dataL1,paddingL);
paddingR=zeros(round(size(filtered_dataR1,1)/10), size(filtered_dataR1,2));
filtered_dataR2=cat(1,paddingR,filtered_dataR1,paddingR);
% low pass
Hd = makefilter(srL,50,51,6,20,0);  
filtered_dataL3=filtfilthd(Hd,filtered_dataL2);
filtered_dataR3=filtfilthd(Hd,filtered_dataR2);
% remove padding
filtered_dataL5=filtered_dataL3((size(paddingL,1)+1):(size(paddingL,1)+size(detrend_dataL,1)),:);
filtered_dataR5=filtered_dataR3((size(paddingR,1)+1):(size(paddingR,1)+size(detrend_dataR,1)),:);
% clear workspace
clearvars filtered_dataL1 filtered_dataL2 filtered_dataL3 paddingL
clearvars filtered_dataR1 filtered_dataR2 filtered_dataR3 paddingR
%% add FASTICA repo
addpath(genpath('/home/zhibinz2/Documents/GitHub/matlab/external/FastICA_25'))
%% ICA for L Left player
% run ICA
[icasigL, AL, WL] = fastica(filtered_dataL5');
% compute correlation with FP1 and FP2;
[RHO1,PVAL1] = corr(filtered_dataL5(:,1),icasigL');
[RHO3,PVAL3] = corr(filtered_dataL5(:,3),icasigL');
[B1,I1]=sort(abs(RHO1),'descend');[B3,I3]=sort(abs(RHO3),'descend');
ComponentsExamL=unique([I1(1:5) I3(1:5)])
% Choose component to remove based on my removal criteria
TentativeRemoveL=[];
for i=1:length(ComponentsExamL)
    ComponentExam=AL(:,ComponentsExamL(i));
    if (abs(ComponentExam(1))+abs(ComponentExam(3)))>(((2/32)*sum(abs(ComponentExam)))*3); % removal criteria
        TentativeRemoveL=[TentativeRemoveL i];
    end
end
ComponentRemoveL=ComponentsExamL(TentativeRemoveL);
%% ICA for Right player
% run ICA
[icasigR, AR, WR] = fastica(filtered_dataR5');
% compute correlation with FP1 and FP2;
[RHO1,PVAL1] = corr(filtered_dataR5(:,1),icasigR');
[RHO3,PVAL3] = corr(filtered_dataR5(:,3),icasigR');
[B1,I1]=sort(abs(RHO1),'descend');[B3,I3]=sort(abs(RHO3),'descend');
ComponentsExamR=unique([I1(1:5) I3(1:5)])
% Choose component to remove based on my removal criteria
TentativeRemoveR=[];
for i=1:length(ComponentsExamR)
    ComponentExam=AR(:,ComponentsExamR(i));
    if (abs(ComponentExam(1))+abs(ComponentExam(3)))>(((2/32)*sum(abs(ComponentExam)))*3); % removal criteria
        TentativeRemoveR=[TentativeRemoveR i];
    end
end
ComponentRemoveR=ComponentsExamR(TentativeRemoveR);
%% Remove ICA components for Left player
ALrm=AL;icasigLrm=icasigL; % make new, backup AL, icasigL
ALrm(:,ComponentRemoveL)=0; icasigLrm(ComponentRemoveL,:)=0;
mixedsigL=ALrm*icasigLrm;
mixedsigL=mixedsigL';
%% Remove ICA components for Right player
ARrm=AR;icasigRrm=icasigR; % make new, backup AR, icasigR
ARrm(:,ComponentRemoveR)=0; icasigRrm(ComponentRemoveR,:)=0;
mixedsigR=ARrm*icasigRrm;
mixedsigR=mixedsigR';
