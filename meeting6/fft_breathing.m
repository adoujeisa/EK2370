clear all;
close all;
clc;

% Read the real and imaginary parts of the signal
[I,Fs] = audioread('SDR_CW_BREATHING_MIGUEL_REAL.wav'); 
[Q,Fs] = audioread('SDR_CW_BREATHING_MIGUEL_IMAG.wav'); 

% Composite the complex value
data1 = complex(I, Q);
data2 = conj(data1);

% Parameters
c = 3e8;                % Speed of light [m/s]
f_center = 2.45e9;             % Center Frequency [Hz]
Tp = 0.1;                     % Pulse width [s]
N = Tp * Fs;                  % Number of samples per pulse
nTargets = 1;

% Compute the phase of the complex signal
phase_data = angle(data1);  % Extract the phase information (in radians)

% Compute range using the provided formula R_i = c * phi / (4 * pi * f_0)
range = (c * phase_data) / (4 * pi * f_center);

% Define time axis for plotting (based on sampling rate)
time = (0:length(phase_data)-1) / Fs;




% Compute and plot the FFT of the composite signal
dataFFT = fft(data1, size(data1, 1) * 4); % Use a longer FFT length for better resolution
dataFFT = dataFFT(1:floor(end/2)); % Take only the positive frequency components
delta_f = linspace(0, Fs/2, length(dataFFT)); % Frequency axis

% Convert to dB scale
dataFFT_dB = 20 * log10(abs(dataFFT) + eps); % Add eps to avoid log(0)

%% Plot the FFT
figure();
plot(delta_f, dataFFT_dB);
xlabel('Frequency [Hz]');
ylabel('Amplitude [dB]');
title('FFT of Composite Signal');
xlim([0, 2]); 
ylim([-10 30]);  % Set y-axis limit from -0.025 to 0 meters
grid on;
