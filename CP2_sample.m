load 'CP2_SoundClip.mat'

Fs = 44100; %sample rate of the sound clip

t = length(y)/Fs; %record time in seconds
plot([1:length(y)]/Fs, y)
xlabel('Time [sec]')
ylabel('Amplitude')
p8 = audioplayer(y, Fs); playblocking(p8);