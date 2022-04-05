PsychPortAudio('GetDevices')

PsychPortAudio('Stop')

PsychPortAudio('RunMode')

open BasicSoundOutputDemo

BasicSoundOutputDemo(1)

Beeper(400,0.4,1)

 
%% other means of sound output in Matlab like sound(), soundsc(), wavplay(), audioplayer()
% sound() example
load handel
sound(y,Fs)

% audioplayer() example
% only the first three seconds.
load handel;
p = audioplayer(y, Fs);
play(p, [1 (get(p, 'SampleRate') * 3)]);

