%% TMSi (32 channels)
AllchanNames={'Fp1','Fpz','Fp2','F7','F3','Fz','F4','F8','FC5','FC1','FC2','FC6','M1','T7','C3','Cz','C4','T8','M2','CP5','CP1','CP2','CP6','P7','P3','Pz','P4','P8','POz','O1','Oz','O2'};
labels = {' 1-Fp1',' 2-Fpz',' 3-Fp2',' 4-F7',' 5-F3',' 6-Fz',' 7-F4',' 8-F8','9-FC5','10-FC1','11-FC2','12-FC6','13-M1','14-T7','15-C3','16-Cz','17-C4','18-T8','19-M2','20-CP5','21-CP1','22-CP2','23-CP6','24-P7','25-P3','26-Pz','27-P4','28-P8','29-POz','30-O1','31-Oz','32-O2'};
%% include channel 1 (33 channels)
% AllchanNames={'CREF','Fp1','Fpz','Fp2','F7','F3','Fz','F4','F8','FC5','FC1','FC2','FC6','M1','T7','C3','Cz','C4','T8','M2','CP5','CP1','CP2','CP6','P7','P3','Pz','P4','P8','POz','O1','Oz','O2'};
% labels = {'1-CREF',' 2-Fp1',' 3-Fpz',' 4-Fp2',' 5-F7',' 6-F3',' 7-Fz',' 8-F4',' 9-F8','10-FC5','11-FC1','12-FC2','13-FC6','14-M1','15-T7','16-C3','17-Cz','18-C4','19-T8','20-M2','21-CP5','22-CP1','23-CP2','24-CP6','25-P7','26-P3','27-Pz','28-P4','29-P8','30-POz','31-O1','32-Oz','33-O2'};
%% 3D coordinates of the 33 channels
X=[0.8333,0.94899,0.99585,0.94899,0.58755,0.65354,0.67069,0.65354,0.58755,0.33402,0.35419,0.35419,0.33402,-0.17365,0,0,0,0,0,-0.17365,-0.33402,-0.36228,-0.36228,-0.33402,-0.58755,-0.65747,-0.68225,-0.65747,-0.58755,-0.89101,-0.98006,-0.99613,-0.98006];
Y=[0,0.30834,0,-0.30834,0.80869,0.52922,0,-0.52922,-0.80869,0.87016,0.35419,-0.35419,-0.87016,0.98481,1,0.70041,0,-0.70041,-1,-0.98481,0.87016,0.34985,-0.34985,-0.87016,0.80869,0.53241,0,-0.53241,-0.80869,0,0.1905,0,-0.1905];
Z=[0.4164,0.065926,0.09098,0.065926,0.028271,0.54112,0.74174,0.54112,0.028271,0.36228,0.8655,0.8655,0.36228,0,0,0.71374,1,0.71374,0,0,0.36228,0.86392,0.86392,0.36228,0.028271,0.53317,0.73112,0.53317,0.028271,0.45399,0.056519,0.087851,0.056519];
%% Reduced to 32 channels
X=X(2:33);Y=Y(2:33);Z=Z(2:33);
%% Calculate subplot coordinates for scalp map
% cd /home/zhibin/Documents/GitHub/Motor_cordination/EEGanalysis
% load('coordinates.mat');
% X = [chanlocs128.X];
% Y = [chanlocs128.Y];
% Z = [chanlocs128.Z];
% labels = {'FP1','FPZ','FP2','AF3','AF4','F11','F7','F5','F3','F1','FZ','F2','F4','F6','F8','F12','FT11','FC5','FC3','FC1','FCZ','FC2','FC4','FC6','FT12','T7','C5','C3','C1','CZ','C2','C4','C6','T8','TP7','CP5','CP3','CP1','CPZ','CP2','CP4','CP6','TP8','M1','M2','P7','P5','P3','P1','PZ','P2','P4','P6','P8','PO7','PO3','POZ','PO4','PO8','O1','OZ','O2','CB1','CB2','AFP1','AFP2','AF7','AF5','AFZ','AF6','AF8','AFF5H','AFF3H','AFF1H','AFF2H','AFF4H','AFF6H','F9','F10','FFT7H','FFC5H','FFC3H','FFC1H','FFC2H','FFC4H','FFC6H','FFT8H','FT9','FT7','FT8','FT10','FTT7H','FCC5H','FCC3H','FCC1H','FCC2H','FCC4H','FCC6H','FTT8H','TTP7H','CCP5H','CCP3H','CCP1H','CCP2H','CCP4H','CCP6H','TTP8H','TPP7H','CPP5H','CPP3H','CPP1H','CPP2H','CPP4H','CPP6H','TPP8H','P9','P10','PPO3H','PPO1H','PPO2H','PPO4H','PO9','PO5','PO1','PO2','PO6','PO10','CBZ'};
% labels = {'1-FP1', '2-FPZ', '3-FP2', '4-AF3', '5-AF4', '6-F11', '7-F7', '8-F5', '9-F3', '10-F1', '11-FZ', '12-F2', '13-F4', '14-F6', '15-F8', '16-F12', '17-FT11', '18-FC5', '19-FC3', '20-FC1', '21-FCZ', '22-FC2', '23-FC4', '24-FC6', '25-FT12', '26-T7', '27-C5', '28-C3', '29-C1', '30-CZ', '31-C2', '32-C4', '33-C6', '34-T8', '35-TP7', '36-CP5', '37-CP3', '38-CP1', '39-CPZ', '40-CP2', '41-CP4', '42-CP6', '43-TP8', '44-M1', '45-M2', '46-P7', '47-P5', '48-P3', '49-P1', '50-PZ', '51-P2', '52-P4', '53-P6', '54-P8', '55-PO7', '56-PO3', '57-POZ', '58-PO4', '59-PO8', '60-O1', '61-OZ', '62-O2', '63-CB1', '64-CB2', '65-AFP1', '66-AFP2', '67-AF7', '68-AF5', '69-AFZ', '70-AF6', '71-AF8', '72-AFF5H', '73-AFF3H', '74-AFF1H', '75-AFF2H', '76-AFF4H', '77-AFF6H', '78-F9', '79-F10', '80-FFT7H', '81-FFC5H', '82-FFC3H', '83-FFC1H', '84-FFC2H', '85-FFC4H', '86-FFC6H', '87-FFT8H', '88-FT9', '89-FT7', '90-FT8', '91-FT10', '92-FTT7H', '93-FCC5H', '94-FCC3H', '95-FCC1H', '96-FCC2H', '97-FCC4H', '98-FCC6H', '99-FTT8H', '100-TTP7H', '101-CCP5H', '102-CCP3H', '103-CCP1H', '104-CCP2H', '105-CCP4H', '106-CCP6H', '107-TTP8H', '108-TPP7H', '109-CPP5H', '110-CPP3H', '111-CPP1H', '112-CPP2H', '113-CPP4H', '114-CPP6H', '115-TPP8H', '116-P9', '117-P10', '118-PPO3H', '119-PPO1H', '120-PPO2H', '121-PPO4H', '122-PO9', '123-PO5', '124-PO1', '125-PO2', '126-PO6', '127-PO10', '128-CBZ'};

sort(Z) %just to examine
ZZ=((max(Z)-min(Z))*1.5-(Z-min(Z)))./((max(Z)-min(Z))*1.5); %Make higher Z smaller fraction
sort(ZZ) %just to examine

XXPLOT=(10*Y)*0.95;
YYPLOT=(10*X)*0.8;

XXPLOT=XXPLOT.*ZZ;
YYPLOT=YYPLOT.*ZZ;

XXPLOT=-1*XXPLOT;
YYPLOT=YYPLOT;

XXPLOT=XXPLOT*0.6;
YYPLOT=YYPLOT*0.7;

% adjust to [0 1]
XXPLOT=XXPLOT/13+0.5;
YYPLOT=YYPLOT/12+0.5;

figure('units','normalized','outerposition',[0 0 1 1]);
plot(XXPLOT,YYPLOT,'r.');

% labels = {'FP1','FPZ','FP2','AF3','AF4','F11','F7','F5','F3','F1','FZ','F2','F4','F6','F8','F12','FT11','FC5','FC3','FC1','FCZ','FC2','FC4','FC6','FT12','T7','C5','C3','C1','CZ','C2','C4','C6','T8','TP7','CP5','CP3','CP1','CPZ','CP2','CP4','CP6','TP8','M1','M2','P7','P5','P3','P1','PZ','P2','P4','P6','P8','PO7','PO3','POZ','PO4','PO8','O1','OZ','O2','CB1','CB2','AFP1','AFP2','AF7','AF5','AFZ','AF6','AF8','AFF5H','AFF3H','AFF1H','AFF2H','AFF4H','AFF6H','F9','F10','FFT7H','FFC5H','FFC3H','FFC1H','FFC2H','FFC4H','FFC6H','FFT8H','FT9','FT7','FT8','FT10','FTT7H','FCC5H','FCC3H','FCC1H','FCC2H','FCC4H','FCC6H','FTT8H','TTP7H','CCP5H','CCP3H','CCP1H','CCP2H','CCP4H','CCP6H','TTP8H','TPP7H','CPP5H','CPP3H','CPP1H','CPP2H','CPP4H','CPP6H','TPP8H','P9','P10','PPO3H','PPO1H','PPO2H','PPO4H','PO9','PO5','PO1','PO2','PO6','PO10','CBZ'};;
text(XXPLOT,YYPLOT,labels,'VerticalAlignment','bottom','HorizontalAlignment','right');
%%
