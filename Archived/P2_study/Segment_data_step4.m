
%% Segment EEG
% mixedsigL=mixedsigL';
% mixedsigR=mixedsigR';
% PacersL=PacerTimeIndL([1 2 3 4 5 34 35 36 37 66 67 68 69 98 99 100 101 130 131 132 133 162 163 164 165 194]);
% PacersR=PacerTimeIndR([1 2 3 4 5 34 35 36 37 66 67 68 69 98 99 100 101 130 131 132 133 162 163 164 165 194]);

EEGOpenEyeRestingL=mixedsigL(PacersL(1):PacersL(2),:);
EEGCloseEyeRestingL=mixedsigL(PacersL(3):PacersL(4),:);
EEGCondi1L=mixedsigL(PacersL(6):PacersL(7),:);
EEGResting2L=mixedsigL(PacersL(7):PacersL(8),:);
EEGCondi2L=mixedsigL(PacersL(10):PacersL(11),:);
EEGResting3L=mixedsigL(PacersL(11):PacersL(12),:);
EEGCondi3L=mixedsigL(PacersL(14):PacersL(15),:);
EEGResting4L=mixedsigL(PacersL(15):PacersL(16),:);
EEGCondi4L=mixedsigL(PacersL(18):PacersL(19),:);
EEGResting5L=mixedsigL(PacersL(19):PacersL(20),:);
EEGCondi5L=mixedsigL(PacersL(22):PacersL(23),:);
EEGResting6L=mixedsigL(PacersL(23):PacersL(24),:);
EEGCondi6L=mixedsigL(PacersL(26):end,:);

EEGOpenEyeRestingR=mixedsigR(PacersR(1):PacersR(2),:);
EEGCloseEyeRestingR=mixedsigR(PacersR(3):PacersR(4),:);
EEGCondi1R=mixedsigR(PacersR(6):PacersR(7),:);
EEGResting2R=mixedsigR(PacersR(7):PacersR(8),:);
EEGCondi2R=mixedsigR(PacersR(10):PacersR(11),:);
EEGResting3R=mixedsigR(PacersR(11):PacersR(12),:);
EEGCondi3R=mixedsigR(PacersR(14):PacersR(15),:);
EEGResting4R=mixedsigR(PacersR(15):PacersR(16),:);
EEGCondi4R=mixedsigR(PacersR(18):PacersR(19),:);
EEGResting5R=mixedsigR(PacersR(19):PacersR(20),:);
EEGCondi5R=mixedsigR(PacersR(22):PacersR(23),:);
EEGResting6R=mixedsigR(PacersR(23):PacersR(24),:);
EEGCondi6R=mixedsigR(PacersR(26):end,:);

save(['EEG' num2str(seed) '.mat'],'EEGCondi1L','EEGCondi2L','EEGCondi3L',...
    'EEGCondi4L','EEGCondi5L','EEGCondi6L',...
    'EEGCondi1R','EEGCondi2R','EEGCondi3R',...
    'EEGCondi4R','EEGCondi5R','EEGCondi6R');


%% segment button presses
BottonPresTimeL01; BottonPresTimeR01;

BPOpenEyeRestingL=BottonPresTimeL01(PacersL(1):PacersL(2),:);
BPCloseEyeRestingL=BottonPresTimeL01(PacersL(3):PacersL(4),:);
BPCondi1L=BottonPresTimeL01(PacersL(6):PacersL(7),:);
BPResting2L=BottonPresTimeL01(PacersL(7):PacersL(8),:);
BPCondi2L=BottonPresTimeL01(PacersL(10):PacersL(11),:);
BPResting3L=BottonPresTimeL01(PacersL(11):PacersL(12),:);
BPCondi3L=BottonPresTimeL01(PacersL(14):PacersL(15),:);
BPResting4L=BottonPresTimeL01(PacersL(15):PacersL(16),:);
BPCondi4L=BottonPresTimeL01(PacersL(18):PacersL(19),:);
BPResting5L=BottonPresTimeL01(PacersL(19):PacersL(20),:);
BPCondi5L=BottonPresTimeL01(PacersL(22):PacersL(23),:);
BPResting6L=BottonPresTimeL01(PacersL(23):PacersL(24),:);
BPCondi6L=BottonPresTimeL01(PacersL(26):end,:);

BPOpenEyeRestingR=BottonPresTimeR01(PacersR(1):PacersR(2),:);
BPCloseEyeRestingR=BottonPresTimeR01(PacersR(3):PacersR(4),:);
BPCondi1R=BottonPresTimeR01(PacersR(6):PacersR(7),:);
BPResting2R=BottonPresTimeR01(PacersR(7):PacersR(8),:);
BPCondi2R=BottonPresTimeR01(PacersR(10):PacersR(11),:);
BPResting3R=BottonPresTimeR01(PacersR(11):PacersR(12),:);
BPCondi3R=BottonPresTimeR01(PacersR(14):PacersR(15),:);
BPResting4R=BottonPresTimeR01(PacersR(15):PacersR(16),:);
BPCondi4R=BottonPresTimeR01(PacersR(18):PacersR(19),:);
BPResting5R=BottonPresTimeR01(PacersR(19):PacersR(20),:);
BPCondi5R=BottonPresTimeR01(PacersR(22):PacersR(23),:);
BPResting6R=BottonPresTimeR01(PacersR(23):PacersR(24),:);
BPCondi6R=BottonPresTimeR01(PacersR(26):end,:);

save(['BP' num2str(seed) '.mat'],'BPCondi1L','BPCondi2L','BPCondi3L',...
    'BPCondi4L','BPCondi5L','BPCondi6L',...
    'BPCondi1R','BPCondi2R','BPCondi3R',...
    'BPCondi4R','BPCondi5R','BPCondi6R');


%% segment feedbacks from the other 
FeedbTimeL01; FeedbTimeR01;

FBOpenEyeRestingL=FeedbTimeL01(PacersL(1):PacersL(2),:);
FBCloseEyeRestingL=FeedbTimeL01(PacersL(3):PacersL(4),:);
FBCondi1L=FeedbTimeL01(PacersL(6):PacersL(7),:);
FBResting2L=FeedbTimeL01(PacersL(7):PacersL(8),:);
FBCondi2L=FeedbTimeL01(PacersL(10):PacersL(11),:);
FBResting3L=FeedbTimeL01(PacersL(11):PacersL(12),:);
FBCondi3L=FeedbTimeL01(PacersL(14):PacersL(15),:);
FBResting4L=FeedbTimeL01(PacersL(15):PacersL(16),:);
FBCondi4L=FeedbTimeL01(PacersL(18):PacersL(19),:);
FBResting5L=FeedbTimeL01(PacersL(19):PacersL(20),:);
FBCondi5L=FeedbTimeL01(PacersL(22):PacersL(23),:);
FBResting6L=FeedbTimeL01(PacersL(23):PacersL(24),:);
FBCondi6L=FeedbTimeL01(PacersL(26):end,:);

FBOpenEyeRestingR=FeedbTimeR01(PacersR(1):PacersR(2),:);
FBCloseEyeRestingR=FeedbTimeR01(PacersR(3):PacersR(4),:);
FBCondi1R=FeedbTimeR01(PacersR(6):PacersR(7),:);
FBResting2R=FeedbTimeR01(PacersR(7):PacersR(8),:);
FBCondi2R=FeedbTimeR01(PacersR(10):PacersR(11),:);
FBResting3R=FeedbTimeR01(PacersR(11):PacersR(12),:);
FBCondi3R=FeedbTimeR01(PacersR(14):PacersR(15),:);
FBResting4R=FeedbTimeR01(PacersR(15):PacersR(16),:);
FBCondi4R=FeedbTimeR01(PacersR(18):PacersR(19),:);
FBResting5R=FeedbTimeR01(PacersR(19):PacersR(20),:);
FBCondi5R=FeedbTimeR01(PacersR(22):PacersR(23),:);
FBResting6R=FeedbTimeR01(PacersR(23):PacersR(24),:);
FBCondi6R=FeedbTimeR01(PacersR(26):end,:);

save(['FB' num2str(seed) '.mat'],'FBCondi1L','FBCondi2L','FBCondi3L',...
    'FBCondi4L','FBCondi5L','FBCondi6L',...
    'FBCondi1R','FBCondi2R','FBCondi3R',...
    'FBCondi4R','FBCondi5R','FBCondi6R');

%% %% segment EMG
filtered_EMGL; filtered_EMGR; % before hilbert transform
upL;upR; % after hilbert transform

EMGOpenEyeRestingL=upL(PacersL(1):PacersL(2),:);
EMGCloseEyeRestingL=upL(PacersL(3):PacersL(4),:);
EMGCondi1L=upL(PacersL(6):PacersL(7),:);
EMGResting2L=upL(PacersL(7):PacersL(8),:);
EMGCondi2L=upL(PacersL(10):PacersL(11),:);
EMGResting3L=upL(PacersL(11):PacersL(12),:);
EMGCondi3L=upL(PacersL(14):PacersL(15),:);
EMGResting4L=upL(PacersL(15):PacersL(16),:);
EMGCondi4L=upL(PacersL(18):PacersL(19),:);
EMGResting5L=upL(PacersL(19):PacersL(20),:);
EMGCondi5L=upL(PacersL(22):PacersL(23),:);
EMGResting6L=upL(PacersL(23):PacersL(24),:);
EMGCondi6L=upL(PacersL(26):end,:);

EMGOpenEyeRestingR=upR(PacersR(1):PacersR(2),:);
EMGCloseEyeRestingR=upR(PacersR(3):PacersR(4),:);
EMGCondi1R=upR(PacersR(6):PacersR(7),:);
EMGResting2R=upR(PacersR(7):PacersR(8),:);
EMGCondi2R=upR(PacersR(10):PacersR(11),:);
EMGResting3R=upR(PacersR(11):PacersR(12),:);
EMGCondi3R=upR(PacersR(14):PacersR(15),:);
EMGResting4R=upR(PacersR(15):PacersR(16),:);
EMGCondi4R=upR(PacersR(18):PacersR(19),:);
EMGResting5R=upR(PacersR(19):PacersR(20),:);
EMGCondi5R=upR(PacersR(22):PacersR(23),:);
EMGResting6R=upR(PacersR(23):PacersR(24),:);
EMGCondi6R=upR(PacersR(26):end,:);

save(['EMG' num2str(seed) '.mat'],'EMGCondi1L','EMGCondi2L','EMGCondi3L',...
    'EMGCondi4L','EMGCondi5L','EMGCondi6L',...
    'EMGCondi1R','EMGCondi2R','EMGCondi3R',...
    'EMGCondi4R','EMGCondi5R','EMGCondi6R');