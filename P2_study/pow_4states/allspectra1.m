function [pow,freqs,df,corr,ampcorr,fcoef] = allspectra1(data,rate,maxfreq,goodepochs,win,varargin)
%% DESCRIPTION-
%% function to extract speatral measures from any 
%% segmented data in a SEGEEG structure .m
%% 
%% spectra = allspectra(data,maxfreq,goodepochs,win);
%% 
%% INPUTS-
%% data: data matrix, time x channel x trial 
%% maxfreq: maximum frequency to keep (Hz), default = 50 
%% goodepochs, default = all trials
%% win:sample range to analyze, default = entire epoch
%% 
%% OUPUTS-
%% pow: power using variance (variance of relative power) -relative form
%% corr: complex valued correlation cofficiences - use abs(corr).^2 to get coherence 
%% ampcorr: amplitude correlation -
%% fcoef: Fourier coefficients 
%% 
%% NOTES-
%%
%% 
%

if nargin < 3 || isempty(maxfreq)
  maxfreq = 50;
end
if nargin < 4 || isempty(goodepochs)
  goodepochs = 1:size(data,3);
end
if nargin < 5 || isempty(win)
  win = 1:size(data,1);
end

df = rate/length(win);
nbins = ceil(maxfreq/df) + 1;
freqs = (0:(nbins-1))*df;
fcoef = fft(ndetrend(data(win,:,:),1),[],1)/length(win);
eppow = abs(squeeze(mean(fcoef(1:nbins,:,goodepochs),3))).^2; % absolute power -abs - averaged across all good trials
pow = squeeze(var(fcoef(1:nbins,:,goodepochs),[],3)); % relative power - variance
if nargout > 4
for k = 1:nbins % for each frequency
   corr(k,:,:) = corrcoef(transpose(squeeze(fcoef(k,:,goodepochs)))); % coherence = correlation cofficiences
   ampcorr(k,:,:) = corrcoef(transpose(squeeze(abs(fcoef(k,:,goodepochs)))));% amplitude correlation
end
end



