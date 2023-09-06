load 'CP2_SoundClip.mat'

Fs = 44100; %sample rate of the sound clip

S = y'; % transposes in order to have consistent dimensions
w = length(y)/4; % break the spectogram up into four time windows, otherwise it
% will be too big for MATLAB/autograder to run.

S1 = S((1-1)*w+1:1*w); % this will isolate the correct window
S2 = S((2-1)*w+1:2*w); % this will isolate the correct window
S3 = S((3-1)*w+1:3*w); % this will isolate the correct window
S4 = S((4-1)*w+1:4*w); % this will isolate the correct window

L = length(S1)/Fs; % length of each window in seconds
n = length(S1);  % number of elements in each window
t = [0:1/Fs:L - 1/Fs]; % t in sec. relative to the start of the window
tau = 0:0.1:L; % discretization for the Gabor transform
k = 2*pi*(1/L/2)*[0:n/2-1 -n/2:-1]; % discretization in frequency space
ks = fftshift(k); % gotta shift them freqs.

Sgt_spec = zeros(length(ks),length(tau)); % initializing the function for the spectrogram

% Repeat this part for S1, S2, S3, and S4.
% For this part we are going to make a spectrogram, but only for the freqs.
% of interest.  So first we'll do a Gabor transform, then we'll filter
% around our peak freq. in the regime of interest, and then we'll look at
% the spectrogram of that function.

% Gabor transform each S, just like we did in the lecture
% Week4_Spectrograms.m lines 141 to 146.
% You'll have to add code between line 144 and Sgt_spec at line 145.

% After line 144 find the index of your peak frequency (in absolute value)
% within the range of interest (this range is very forgiving so you don't
% have to match the autograder exactly, just use your judgement based on
% the figure in the assignment statement).  Then make a filter centered
% about the peak frequency.  Filter the Gabor transformed function.
% this is the function you will use in line 145 from the lecture to find
% your Sgt_spec.

%% Make the Spectrogram

windows = [S1; S2; S3; S4];
a = 400;   %window width for the Gabor function g
c = 1/L;   %coeff. for the filter function f

for window = 1:4
    Sgt_spec = zeros([length(ks),length(tau)]);
    for j = 1:length(tau)
        g = exp(-a * (t - tau(j)).^2);
        Sg = windows(window, :) .* g;
        Sgt = fft(Sg);
        [maxFreq, idx] = max(Sgt(1:1800));
        f = exp(-c .* (abs(k) - abs(k(idx))).^2);   %filter function
        Sgtf = Sgt .* f;   %multiply abs of frequencies of signal to the filter
        Sgt_spec(:,j) = fftshift(abs(Sgtf));
    end
    
    if window == 1
        A1 = Sgt_spec;
    elseif window == 2
        A2 = Sgt_spec;
    elseif window == 3
        A3 = Sgt_spec;
    elseif window == 4
        A4 = Sgt_spec;
    end
    
end

% Save Sgt_spec as variable A1 after your for loop.  Repeat for S2, S3, and
% S4 and save those Sgt_spec as A2, A3, and A4.  You don't have to rewrite
% the code, just copy and paste and use the respective S's, or write a for
% loop that iterates through S1 to S4.

%% Plotting Stuff
% Plot of spectrogram for each window (for the report, not for autograder)
% just like we did in the lecture, but change your ylim to be in the range
% of interest for the sound clip.

for window = 1:4
    figure(window)
    subplot(1, 1, 1)
    pcolor(tau,abs(ks),log(specList(:, :, window) + 1))
    shading interp
    set(gca,'ylim',[0 700],'Fontsize',16)
    colormap(hot)
    xlabel('time (t)'), ylabel('frequency (k)')
    title(['window = S',num2str(window)],'Fontsize',16)
    pause
    close all
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Isolate the bassline
S = y'; % transposes in order to have consistent dimensions
L = length(S)/Fs; % total length in sec. 
n = length(S); % total number of elements in S
t = 0:1/Fs:L - 1/Fs; % time discretization
k = (1/L)*[0:n/2-1 -n/2:-1]; % freq. discretization

% Take the Fourier transform of S, and in freq. space isolate all freqs. 
% (in absolute value) that you determine should be part of the baseline 
% according to spectrogram (or also just by listening); that is, all points
% in the transformed function not within the frequency range you determined
% should be set to zero (kind of like a Shannon filter, but simpler than
% what we did in lecture).
% You may have to do this part a few times with different thresholds to get
% it right.

St = fft(S);
newSt = St;
newSt(abs(k) > 150) = 0;
newS = ifft(newSt);

% After thresholding the transformed function, take the inverse transform
% of the thresholded function and save it as A5.

A5 = newS';     %Shape:  1x1938240 double

%% Play sound (not for autograder)

p8 = audioplayer(newS, Fs); playblocking(p8);

%% Plotting Stuff
%Plot the amplitude S over time (for the report, not for the autograder)

plot(t, newS)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Isolate the guitar
%Same exact process as the baseline above, but you'll have to be more
%careful about the frequency range.
S = y'; %reinitialize the S from the previous part above.

St = fft(S);
newSt = St;
newSt(abs(k) < 250 | abs(k) > 1200) = 0;
newS = ifft(newSt);

A6 = newS';     %Shape:  1x1938240 double

%% Play Sound
%Play sound (not for autograder)

p8 = audioplayer(newS, Fs); playblocking(p8);

%% Plotting Stuff
%Plot the amplitude S over time (for the report, not for the autograder)

plot(t, newS)
