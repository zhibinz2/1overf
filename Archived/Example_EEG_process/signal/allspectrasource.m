function [pow,freqs,df,eppow, corr,cprod,sourceorient,fcoefsource] = allspectrasource(data,rate,maxfreq,goodepochs,win,varargin);
%% DESCRIPTION-
%% function to extract speatral measures from any 
%% segmented data in a SEGEEG structure .m
%%
%% spectra = allspectra(SEGEEG,maxfreq,goodepochs,win);
%%
%% INPUTS-
%% data: data matrix, time x channel x trial 
%% maxfreq: maximum frequency to keep (Hz), default = 50 
%% goodepochs, default = all trials
%% win:sample range to analyze, default = entire epoch
%% 
%%
%% NOTES-
%%
%% 
%

if nargin < 3 | isempty(maxfreq),
  maxfreq = 50;
end;
if nargin < 4 | isempty(goodepochs),
  goodepochs = [1:size(data,3)];
end;
if nargin < 5 | isempty(win)
  win = [1:size(data,1)];
end;

df = rate/length(win);
nbins = ceil(maxfreq/df) + 1;
freqs = [0:(nbins-1)]*df;
fcoef = fft(ndetrend(data(win,:,:),1),[],1)/length(win);
nsources = size(fcoef,2)/3
fcoefsource = zeros(nbins,nsources,size(fcoef,3))
for k = 1:nbins 
	for l = 1:nsources
		sourcedata = squeeze(fcoef(k,(l-1)*3+1:l*3,:));
		[u,s,v] = svd(sourcedata,'econ');
		sourceorient(:,l)= u(:,1);
		fcoefsource(k,l,:) = v(:,1);
	end
end
eppow = abs(squeeze(mean(fcoefsource(1:nbins,:,goodepochs),3))).^2;
pow = squeeze(var(fcoefsource(1:nbins,:,goodepochs),[],3));
if nargout > 4
for k = 1:nbins
   sf = corrcoef(transpose(squeeze(fcoefsource(k,:,goodepochs))));
%  coh(k,:,:) = abs(sf).^2;
%  phase(k,:,:) = angle(sf);
  cprod(k,:,:) = cov(transpose(squeeze(fcoefsource(k,:,goodepochs))));
  corr(k,:,:) = sf;
end;
end;



