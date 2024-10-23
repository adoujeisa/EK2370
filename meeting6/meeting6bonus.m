% Load SDR data from .wav file
[data, fs] = audioread('CW_1023.WAV');

% Define constants
c = 3e8;  % Speed of light (m/s)
f1 = 2.45e9;  % Frequency 1 (Hz)
f2 = f1 + 50e6;  % Frequency 2 with a 50 MHz step

% Time vector for the data
t = (0:length(data)-1) / fs;

% Process the data with two-step frequencies
% Assume data is alternating between two frequencies
data_f1 = data(1:2:end);  % Extract samples at f1
data_f2 = data(2:2:end);  % Extract samples at f2

% Calculate the phase difference between two frequency steps
phase_f1 = angle(hilbert(data_f1));
phase_f2 = angle(hilbert(data_f2));
phase_diff = unwrap(phase_f2 - phase_f1);

% Calculate the range based on the phase difference
% Phase difference is related to range by: range = (c * phase_diff) / (4 * pi * freq_step)
freq_step = f2 - f1;  % Step size in Hz (50 MHz)
range = (c * phase_diff) /( 4 * pi * freq_step);

% Plot the estimated range over time
figure;
plot(t(1:2:end), range);
xlabel('Time (s)');
ylabel('Range (m)');
title('Estimated Range from SDR Data');
grid on;
