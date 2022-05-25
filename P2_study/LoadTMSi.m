function [time,samples,TRIGGERind,sr,channels_info] = LoadTMSi(Path_filename);
    d = TMSiSAGA.Poly5.read(Path_filename);

    samples=d.samples;
    sr=d.sample_rate;

    channels=d.channels;
    numbers=num2str([1:length(channels)]');
    labels=strings(length(channels),1);
    units=strings(length(channels),1);
    for i=1:length(channels)
        labels(i)=channels{i}.alternative_name;
        units(i)=channels{i}.unit_name;
    end
    channels_info=table(numbers,labels,units)

    % Create time stamps
    % num2str(d.time)
    time=[1/sr:1/sr:d.time]';

    % Plot channels of Key presses, photocells, EMG
    % look for TRIGGERS channel;
    TRIGGERind=find(labels=='TRIGGERS');
end
