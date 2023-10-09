
%% 
% # The sound becomes more shrill. A fourier transform can't be used to obtain 
% temporal behaviour. Create the signal, plot the spectrum. Temporal behavior 
% is supposed to be needed. Deep blue (low) red (high). A line going upward as 
% a function of time. Windowing is used, to quench the spectral leakage. Chirp 
% freq keeps inc from one value to another value. 
%% 
% frequency resolution: more data samples, window length increase will improve. 
% The contours on heat map becomes more and more focused (else it's vague). 
% 
% temporal resolution: 
% 
% time resolution: tradeoff with frequency resoution. 
% 
% 2.  window length 150 instead of 100, what it does is 
% 
% 3. 

%Question 1
alpha = 2;
Fs = 100;                 % Sampling rate (samples per second)
duration = 10;            % Duration of the signal (seconds)
t = 0:1/Fs:duration;      % Time vector

% Generate chirp signal
F_start = 2 + 2*alpha  % Initial frequency (Hz)
F_end = 5 + 5*alpha    % Final frequency (Hz)
frequencies = linspace(F_start, F_end, length(t));
x = sin(2*pi*frequencies.*t);

% 1. Plot the signal as a function of time
figure;

plot(t, x);
%title('Chirp Signal x(t)');
xlabel('Time (s)');
ylabel('Amplitude');

% 2. Compute and plot the frequency spectrum using FFT
N = length(x);            % Length of the signal
frequencies = Fs*(0:(N/2))/N;
X = fft(x);
X_magnitude = abs(X(1:N/2+1));


plot(frequencies, X_magnitude);
%title('Frequency Spectrum');
xlabel('Frequency (Hz)');
ylabel('Magnitude');

% 3. Plot the spectrogram of the signal
window_length = 100;      % Hamming window length (samples)
overlap = 10;             % Overlap (samples)


spectrogram(x, blackman(window_length), overlap, Fs, 'yaxis');
%title('Spectrogram');
%% 
% This problem aims to explore and analyze a chirp signal generated with a linearly 
% increasing frequency. We plot a spectrogram of the chirp signal, providing a 
% time-frequency representation, and compare it with the FFT analysis.
% 
% Methodology: The chirp signal is generated and its fourier tranform obtained 
% using the fft() function. The plot is then used to identify the different frequency 
% components of the signal. 
% 
% The spectogram() function is used to obtain the frequency-time plot of the 
% signal, the results are analysed by varying window length and using different 
% window techniques. 
% 
% Results: The DFT representation shows the entire frequency spectrum of the 
% signal in a single plot. It can provide a clear picture of the frequency content 
% of the entire signal.The spectrogram representation provides a time-frequency 
% analysis of the signal. It can reveal how the frequency content evolves over 
% time, which is particularly useful for signals with changing characteristics 
% like chirps. The choice of window length and windowing technique affects the 
% spectrogram's resolution in both time and frequency domains.
% 
% 
% 
% Window Length:
% 
% A longer window length in the spectrogram provides better frequency resolution 
% but poorer time resolution. It can capture finer details in the frequency domain 
% but may miss rapid changes in the signal.
% 
% A shorter window length provides better time resolution but poorer frequency 
% resolution. It can capture rapid changes in the signal but might not reveal 
% fine details in the frequency domain. So there's a tradeoff between the both. 
% 
% 
% 
% Windowing Techniques:
% 
% The hanning window gives a balanced compromise between frequency resolution 
% and sidelobe suppression while  hamming window provides sidelobe suppression 
% over frequency resolution, still giving a reasonable tradeoff.  The blackman 
% window minimizes spectral leakage and reduced frequency resolution.
% 
% Preference:  For a chirp signal, where the frequency changes linearly, the 
% spectrogram is generally preferred because it can clearly show how the frequency 
% evolves over time. The choice of window length and windowing technique depends 
% on the specific analysis goals. A balance between time and frequency resolution 
% should be considered.

%Question 2
% Load the audio file
[audio, Fs] = audioread('instru2.wav');

% Plot the spectrogram
figure;
spectrogram(audio, hann(512), 256, Fs, 'yaxis');

%%
% Question 2 part 1
N = length(audio);
spectrum = fft(audio);
spectrum_magnitude = abs(spectrum(1:N/2+1));
frequencies = Fs * (0:(N/2)) / N;

% Plot the spectrum
figure; grid on;
plot(frequencies, spectrum_magnitude);
xlabel('Frequency (Hz)');
ylabel('Magnitude');

%% 
% Peak frequency: 2369.12 Hz

%Question 2 part 2 
% Load the audio file
[audio, Fs] = audioread('opera.wav');

% Plot the spectrogram
figure;
spectrogram(audio, hann(512), 256, Fs, 'yaxis');

%% 
% This problem aims to to analyze the fundamental pitch in audio recordings 
% using spectrogram and conventional spectrum analysis (FFT) and obtain the pitch 
% of a given recording. 
% 
% Methodology: For each audio file, spectrogram  is used to visualize pitch 
% and identify fundamental pitch from the spectrogram. The effectiveness of spectrogram-based 
% pitch tracking and the conventional fft based pitch extraction is evaluated. 
% The pitch is obtained by locating the lowest dark band which corresponds to 
% the fundamental frequency. it's frequency is then obtained from the pot. 
% 
% Results: 
% 
% The spectrogram provides a more detailed view of pitch variations in "opera.wav" 
% compared to conventional frequency analysis. The conventional FFT analysis, 
% gives a static representation of the frequency content, and it does not capture 
% how the pitch evolves over time. 
% 
% On the other hand, a spectrogram is a time-frequency representation of the 
% audio signal. It provides a dynamic view of how the frequency content of the 
% signal changes over time. A spectrogram breaks the audio signal into small time 
% windows and computes the FFT for each window. By doing this, it reveals how 
% the pitch (fundamental frequency) and its variations change as the music progresses.
% 
% Obtained peak frequency - 2369.12 Hz - represents the pitch of the instument 
% in audio  instru2.wav 
% 
% Pitch is primarily associated with the fundamental frequency of a sound produced 
% by a musical instrument. The fundamental frequency is the lowest frequency component 
% in a complex sound wave. So, pitch can be defined as the  perceived highness 
% or lowness of a musical note.
% 
% In musical instruments, it corresponds to the primary vibration mode of the 
% instrument's vibrating element, such as a string, air column, or membrane.  
% Different musical instruments produce a range of pitches based on their design 
% and construction. Factors like the length, tension, thickness, and material 
% of the vibrating element determine the instrument's pitch range and timbre.
% 
% Pitch is organized into a system of octaves, where each octave represents 
% a doubling of the frequency. Notes within the same octave have a harmonic relationship 
% and sound similar but differ in frequency.

%Question 3
Fs = 4000;  % Sampling rate (Hz)
recObj = audiorecorder(Fs, 16, 1);

disp('Recording...');
recordblocking(recObj, 2);
disp('Recording complete.');

% Save the recorded audio as a .wav file
audioData = getaudiodata(recObj);
audiowrite('your_name.wav', audioData, Fs);
%%
% Load the recorded audio file (replace 'your_name.wav' with the actual filename)
[audio, Fs] = audioread('your_name.wav');

% Create a spectrogram
figure;
spectrogram(audio, hamming(512), 256, [], Fs, 'yaxis');
xlabel('Time (s)');
ylabel('Frequency (Hz)');

% Analyze the spectrogram to map the location of phonemes
% You can visually identify patterns in the spectrogram corresponding to different phonemes and sounds in your name.
%% 
% Results:
% 
% The color intensity (often represented as a colormap) indicates the magnitude 
% or power of the frequencies. Darker regions indicate higher energy in the corresponding 
% frequency range.
% 
% Phenome mapping: 
% 
% Vowels: are characterized by distinct formants, appearing as dark horizontal 
% bands in the spectrogram.
% 
% Consonants: often result in transient bursts of energy at specific frequencies 
% and are seen as vertical lines or short-duration events.
% 
% Transitions: appear as sloping patterns connecting vowels and consonants.
% 
% Methodology: The getaudio() function is first used to record the audio which 
% is then written into the folder. The spectogram() function is used to plot the 
% spectogram using a hamming window for the same. 
% 
% Aim: To enhance understandin of speech signal analysis by recording and analyzing 
% one's name, focusing on interpreting spectrogram representations.
% 
% 
% 
% 
% 
%