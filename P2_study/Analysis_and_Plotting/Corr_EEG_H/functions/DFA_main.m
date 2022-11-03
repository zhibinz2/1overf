function [D,Alpha1,n,F_n,FitValues]=DFA_main(DATA)
% DATA should be a time series of length(DATA) greater than 2000,and of column vector.
% But aslo works for 100 samples in my data.
% A is the alpha in the paper
% (https://aip.scitation.org/doi/pdf/10.1063/1.166141)
% D is the dimension of the time series % ?
% n can be changed to your interest

% n=100:100:1000; % win_lengths (function default)
% n=2:2:200; % all win_lengths(customized to our need)
n=5:5:length(DATA)/2; % all 32 different win_lengths (customized to our need)
Nw=length(n); % number of different win_lengths
F_n=zeros(Nw,1); % initialize RMS values in the FDA time series
 for i=1:Nw
     F_n(i)=DFA(DATA,n(i),1);
 end
n=n';
 
%  plot(n,F_n);
% plot(log10(n),log10(F_n),'bx');
% xlabel('Scale [log10(n)]') % win_lengths 
% ylabel('RMS [log10(F-n)]') % RMS values
% title('DFA')

% select scales on power scale with base of 2 [2 4 8 16 32 ...] for the fit
% scaleInds=2.^[1:floor(log2(length(n)))];
% % open nextpow2
    
% A=polyfit(log10(n(scaleInds)),log10(F_n(scaleInds)),1);
A=polyfit(log10(n),log10(F_n),1);
% A=polyfit(log10(n(1:10)),log10(F_n(1:10)),1);
Alpha1=A(1); % the slope, the first order polynomial coefficient from polyfit (Hurst Componenet >1 ?)
D=3-A(1); % ? D is the dimension of the time series % ?


% plot the fit
FitValues=polyval(A,log10(n(1:end)));
% FitValues=polyval(A,log10(n(1:10)));

% hold on;
% plot(log10(n(1:end)),FitValues,'r--');
% % plot(log10(n(1:10)),FitValues,'r--');
% legend({'Data',['Fit (DFA=' num2str(Alpha1) ')']},'Location','southeast');

% plot the fit
% FitValues=polyval(A,log10(n(11:end)));
% hold on;
% plot(log10(n(11:end)),FitValues,'b--');
% legend({'Data',['Fit (DFA=' num2str(Alpha1) ')']},'color',[0 0 1]);
return
