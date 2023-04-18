function Hd = makefilter(sr,Fpass,Fstop,Apass,Astop,doplot)
%Hd = makefilter(sr,Fpass,Fstop,Apass,Astop)
%just a simple wrapper on matlab's filter making functions to build butterworht filters for eeg
%Fpass - Passband Frequency
%Fstop - Stopband Frequency
%Apass - Passband Ripple (dB)
%Astop - Stopband Attenuation (dB)
%doplot - set to 1 for bodeplot.  Defautlts to zero
%Fpass and Fstop should be 2 element vectors for bandpass or bandstop filters.  They should have 1 element for lowpass or high pass filters


match = 'passband';  % Band to match exactly

if length(Fpass) == 1
	if Fpass < Fstop
		h  = fdesign.lowpass(Fpass, Fstop, Apass, Astop, sr);
		Hd = design(h, 'butter', 'MatchExactly', match);
	else
		h  = fdesign.highpass(Fstop, Fpass, Astop, Apass, sr);
		Hd = design(h, 'butter', 'MatchExactly', match);
	end
end
if doplot
% 	freqz(Hd,500,'half',sr)
    freqz(Hd,0.01:0.01:60,sr)
end

% [H,F] = freqz(Hd,500,sr);