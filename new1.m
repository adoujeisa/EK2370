clear all;
close all;
clc;

%% COMPUTE VELOCITY FROM AUDIO FILE
% Read the audiofile 
[y, fs] = audioread('Velocity_Test_File.m4a'); 

% Take care of data inversion by the sound card
data = -y(:, 1); % Intensity of the received signal

% Parameters
c = 3e8;                % Speed of light [m/s]
f_center = 2.43e9;            % Center Frequency [Hz]
Tp = 0.1;                     % Pulse width [s]
N = round(Tp * fs);           % Number of samples per pulse
N_total = length(data);       % Total number of samples in read file
T = N_total / fs;             % Total time duration of signals in read file
num_pulse = round(T / Tp);    % Total number of pulses

% Initialize time-domain matrix
mat_time = zeros(num_pulse, N); 

for i = 1:num_pulse-1
    mat_time(i, :) = data(1 + N * (i-1) : i * N); % Assign amplitude
end

%% Mean subtraction clutter removal (MS)
data_mean = mean(data); % Mean of the entire data
mat_time = mat_time - data_mean; % Remove clutter by subtracting the mean
%%


%%zero padding for fft
for i = 1:num_pulse-1
    mat_freq(i,:)= 20*log(abs((fft(mat_time(i,:),4*N))));
 end


%% FFT and Velocity Calculation
fft_data = fft(mat_time, [], 2); % Apply FFT along the second dimension (pulses)
fft_data_db = 20 * log10(abs(fft_data)); % Convert magnitude to dB

 %% Normalize 1
fft_data_norm = fft_data_db - max(max(fft_data_db));
 %% Normalize 2
% Calculate the maximum value for each column
max_vals = max(fft_data_db, [], 2);
% Subtract each column's maximum value from each element in that column
fft_data_norm_2 = fft_data_db - max_vals;


% Velocity array
f_doppler = linspace(0, fs/2, N); % Doppler frequency axis
velocity = f_doppler * c / (2 * f_center); % Convert Doppler frequency to velocity

% Time array
time_array(:,1) = Tp:Tp:Tp*num_pulse;
%time_array = (0:num_pulse-2) * Tp; % Time array for the y-axis

%% Plot using imagesc
figure (1);
imagesc(velocity, time_array, fft_data_norm); % Plot the normalized FFT data
xlabel('Velocity (m/s)');
ylabel('Time (s)');
title('Velocity-Time Map with noormalization 1');
colorbar; % Add color bar to represent intensity

xlim([0 30]); % Set the range for the x-axis (velocity) from -10 to 10 m/s
clim([-45 0]); % Set the color axis range for the colorbar, e.g., from -30 dB to 0 dB

figure (2);
imagesc(velocity, time_array, fft_data_norm_2); % Plot the normalized FFT data
xlabel('Velocity (m/s)');
ylabel('Time (s)');
title('Velocity-Time Map with noormalization 2');
colorbar; % Add color bar to represent intensity

xlim([0 30]); % Set the range for the x-axis (velocity) from -10 to 10 m/s
clim([-15 0]); % Set the color axis range for the colorbar, e.g., from -30 dB to 0 dB