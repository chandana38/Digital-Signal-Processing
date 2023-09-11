
%%
% 1 part 1 and 2 
alpha = 2; 
sampling_rate = 120;
duration = 2; 
frequency = 15 * alpha; 
angular_frequency = 2 * pi * frequency;

% Generate time vector
t = 0 : 1/sampling_rate : duration;

% Generate the sinusoidal signal
x_t = cos(angular_frequency * t);

% Take the first 120 samples
x_samples = x_t(1 : 120);

% Compute the DFT
dft = fft(x_samples);

% Calculate the corresponding frequencies in Hertz
freqs = (0 : length(x_samples) - 1) * sampling_rate / length(x_samples);

% Plot magnitude of the DFT against frequency
figure;
plot(freqs, abs(dft));
xlabel('Frequency (Hz)'); %no spectral leakage 
ylabel('Magnitude');

x_samples_130 = x_t(1 : 130); %spectral leakage
dft_130 = fft(x_samples_130);
freqs = (0 : length(x_samples_130) - 1) * sampling_rate / length(x_samples_130);
plot(freqs, abs(dft_130), 'r', 'DisplayName', 'DFT of 130 samples');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
%%
% part 3
alpha = 2; % Replace with your desired value of alpha
sampling_rate = 120; % Samples per second
duration = 2; % Seconds
frequency = 15 * alpha; % Hz
angular_frequency = 2 * pi * frequency; % Radians per second

% Generate time vector
t = 0 : 1/sampling_rate : duration;

x_t = cos(angular_frequency * t);


x_samples = x_t(1 : 120); %no spectral leakage would be for 240 Hz, periodicity of this signal is 
%origitnal 50, Fs = 150 then periodicity = 3
%so periodcity here would be frequency/sampling_rate, 
frequency/sampling_rate


 
%% 
% # aim: This problem aims at investigating the impact of sample size on the 
% DFT.
%% 
% results: 
%% 
% # For the first 120 samples of the signal, the DFT plot shows a peak at the 
% specified frequency (30 Hz) with a magnitude of 60, as expected for a unit amplitude 
% sinusoid. 
% # For the first 130 samples, it exhibits a peak at the same frequency. The 
% magnitude is now distributed across multiple frequency bins. This indicates 
% a reduction in frequency resolution. 
% # To match the DFT of the first 120 samples with a different sample size (N 
% ≠ 120), we find that N = 240 would produce an identical DFT. This is because, 
% the DFT is periodic with period N, and doubling the sample size aligns it with 
% the original.
%% 
% 
% 
% 

%question 2
% Constants
alpha = 2;
A = 140;
B = 146;
sampling_rate = 200; % Samples per second
duration = 10; % Seconds

% Generate time vector
t = 0 : 1/sampling_rate : duration;

% Generate the signal xa(t)
xa_t = 0.1 * sin(A * pi * t) + cos(B * pi * t);

% Define different sample sizes
sample_sizes = [215, 415, 1115, 1515, 1915];

% Calculate and plot the DFT for each sample size
for i = 1:length(sample_sizes)
    N = sample_sizes(i);
    xa_samples = xa_t(1:N);
    
    dft = fft(xa_samples, N);

    freqs = (0 : N-1) * sampling_rate / N;
   
    figure;
    plot(freqs, abs(dft));
    xlabel('Frequency (Hz)');
    ylabel('Magnitude');
    title(['N = ' num2str(N) '']);
    grid on;
end
%% 
% 2. aim: This problem aims at demonstarating the impact of sample size on frequency 
% resolution and spectral leakage. 
% 
% 
% 
% results: The DFT of the signal with 215 samples results in limited information 
% on frequency components. The main lobe of the DFT is broad which makes it harder 
% to distinguish between closely spaced frequency components. The frequency components 
% present in the signal may not be accurately resolved. 
% 
% On the other hand, the DFT with 1915 samples provides the highest frequency 
% resolution in this analysis. The main lobe is narrow and accurate detection 
% detection of closely spaced frequency components is possible and the spectral 
% leakage is minimal. 
% 
% in conclusion, Increasing the sample size in the DFT analysis results in better 
% frequency resolution and reduced spectral leakage. Smaller sample sizes result 
% in broader main lobes, limiting the ability to distinguish closely spaced frequency 
% components.The choice of sample size is critical for accurate spectral analysis, 
% and a larger N provides more detailed frequency information.
% 
% 

%question 3
% Constants
alpha = 2;
A = 140;
B = 146;
sampling_rate = 200; % Samples per second
duration = 10; % Seconds

% Generate time vector
t = 0 : 1/sampling_rate : duration;

% Generate the signal xa(t)
xa_t = 0.1 * sin(A * pi * t) + cos(B * pi * t);

% Define different sample sizes (N values)
N_values = [215, 415, 1115, 1515, 1915];

% Initialize arrays to store DFT results
dft_results = cell(1, length(N_values));

% Calculate and plot the DFT for each N value with Hanning window
for i = 1:length(N_values)
    N = N_values(i);
    
    % Apply Hanning window
    hanning_window = hanning(N)';
    xa_t_windowed = xa_t(1:N) .* hanning_window;
    
    % Compute the DFT
    dft = fft(xa_t_windowed, N);
    dft_results{i} = dft;
    
    % Calculate the corresponding frequencies in Hertz
    freqs = (0 : N-1) * sampling_rate / N;
    
    % Plot the magnitude of the DFT
    figure;
    plot(freqs, abs(dft));
    xlabel('Frequency (Hz)');
    ylabel('Magnitude');
    title(['N = ' num2str(N) '']);
    grid on;
end


%% 
% 3. results: In Problem 2, without windowing, we observe broader main lobes 
% and pronounced side lobes in the DFT indicating spectral leakage. This hinders 
% the precise identification of frequency components. However, in Problem 3, the 
% Hanning window reduces spectral leakage, resulting in narrower main lobes and 
% diminished side lobes. This enhances the frequency resolution and allows accurate 
% localization and separation of closely spaced frequency components. 
% 
% aim: This problem is aimed at  demonstrating the benefits of windowing techniques 
% in spectral analysis.

%problem 4 part 1
% Load the data from the file
data = load('Exp4Data2.txt');

% Constants
Fs = 1000; % Sampling frequency in Hz (provided in the data)

% Apply Hamming window
hamming_window = hamming(length(data))';
data_windowed = data .* hamming_window;

% Pad zeros to the signal (e.g., 10000 samples)
N = 10000;
data_padded = [data_windowed, zeros(1, N - length(data))];

% Compute the DFT
dft = fft(data_padded);

% Calculate the corresponding frequencies in Hz
freqs = (0 : N-1) * Fs / N;

% Normalize the x-axis
normalized_freqs = freqs / 10000;

% Plot the magnitude of the DFT
figure;
plot(normalized_freqs, abs(dft));
xlabel('Frequency (Fs)');
ylabel('Magnitude');
grid on;

% Find the two dominant frequency components
[sorted_magnitudes, sorted_indices] = sort(abs(dft), 'descend');
dominant_frequencies = normalized_freqs(sorted_indices(1:2));
fprintf('Estimated Dominant Frequencies (in terms of Fs): %.3f, %.3f\n', dominant_frequencies);
%%
%problem 4 part 2
% Load the data from the file
data = load('Exp4Data2.txt');

% Constants
Fs = 1000; % Sampling frequency in Hz (provided in the data)

% Pad zeros to the signal (e.g., 10000 samples)
N = 10000;
data_padded = [data, zeros(1, N - length(data))];

% Compute the DFT
dft = fft(data_padded);

% Calculate the corresponding frequencies in Hz
freqs = (0 : N-1) * Fs / N;

% Normalize the x-axis
normalized_freqs = freqs / 10000;

% Plot the magnitude of the DFT
figure;
plot(normalized_freqs, abs(dft));
xlabel('Frequency (Fs)');
ylabel('Magnitude');
grid on;

% Find the two dominant frequency components
[sorted_magnitudes, sorted_indices] = sort(abs(dft), 'descend');
dominant_frequencies = normalized_freqs(sorted_indices(1:2));
fprintf('Estimated Dominant Frequencies (in terms of Fs): %.3f, %.3f\n', dominant_frequencies);
%% 
% 4. aim: This problem aims at highlighting the impact of windowing and zero-padding 
% techniques on the accuracy of estimating dominant frequencies within a signal. 
% 
% results: The Hamming window is applied to the signal to eliminate spectral 
% leakage and improve frequency resolution. Zero-padding is used to increase the 
% number of samples for better frequency interpolation. The two dominant frequency 
% components in the signal are identified from the DFT, expressed in terms of 
% Fs given as 0.026Fs (26 Hz) and 0.074Fs (74 Hz)
% 
% The choice of windowing technique  impacts the accuracy and resolution of 
% frequency component estimation. The Hamming window, by reducing spectral leakage, 
% provides more precise frequency localization and reduces side lobes in the DFT, 
% making it effective for analyzing signals. The rectangular window, which has 
% no windowing, is simpler but can result in a broader main lobe and increased 
% spectral leakage.