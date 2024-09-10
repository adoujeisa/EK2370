clear all;
close all;
clc;

%% COMPUTE VELOCITY FROM AUDIO FILE
% Read the audio file
[y, fs] = audioread('Range_Test_File.m4a'); 

% Assuming data_1 and data_2 are the correct channels
backscatter = -y(:, 1); % Intensity of the received signal
sync = -y(:, 2);

% Parameters
c = 3e8;                % Speed of light [m/s]
f_start = 2.408e9;      % Start Frequency [Hz]
f_stop = 2.495e9;       % Stop Frequency [Hz]
bandwidth = f_stop - f_start; % Bandwidth [Hz]
delta_r = c / (2*bandwidth); % Range resolution [m]
f_center = (f_start + f_stop) / 2; % Center Frequency [Hz]
Tp = 0.02;              % Pulse width [s]
N = round(Tp * fs);     % Number of samples per pulse
N_total = length(backscatter); % Total number of samples in read file
max_range = (N * delta_r)/2; % Maximum range [m]
num_pulse = floor((N_total - N) / (N * 2)); % Total number of pulses

%% Parse up-chirp data according to the sync data
sync_pulse = zeros(length(sync),2);
sync_pulse(:,1) = (sync > 0); % Set sync square waveform between 0 and 1
sync_pulse(:,2) = 1:1:length(sync_pulse(:,1)); % Set indexes

time_temp = find(sync_pulse==1);
% Initialize time-domain matrix
mat_time = zeros(num_pulse, N);

for i = 1:num_pulse-1
    mat_time(i, :) = backscatter(1 + N * (i-1) : i * N); % Assign amplitude
end

range_data = linspace(0, max_range, 4*N);

%% Build matrix
up_data_parsed = zeros(num_pulse, N); % Pre-allocate 
down_data_parsed = zeros(num_pulse, N);
time = zeros(1, num_pulse); % Pre-allocate
k = 1;
for i = 2:(size(sync_pulse, 1) - 2*N) 
    if sync_pulse(i,1) == 1 && sync_pulse(i-1,1) == 0 % First value of a row = first up-chirp value
        up_data_parsed(k,:) = backscatter(i:i+N-1)';
        down_data_parsed(k,:) = backscatter(i+N:i+2*N-1)';
        time(k) = sync_pulse(i,2) / fs;
        k = k + 1;
    end
end

%% Mean subtraction clutter removal (MS)
data_mean = mean(up_data_parsed); % Mean of the entire data
mat_time = mat_time - data_mean; % Remove clutter by subtracting the mean

%% 2-step MTI
MTI2 = zeros(size(up_data_parsed));
for t = 2:size(up_data_parsed,1)
    MTI2(t,:) = up_data_parsed(t,:) - up_data_parsed(t-1,:);
end

MTI2_final = zeros(size(MTI2));
for t = 2:size(MTI2,1)
    MTI2_final(t,:) = MTI2(t,:) - MTI2(t-1,:);
end

% MTI2 FFT
ifft_MTI2 = 20*log10(abs(fft(MTI2_final, 4*N, 2)));  % dft using zero padding
ifft_MTI2 = ifft_MTI2(:, 1:end/2) - max(max(ifft_MTI2)); % Normalize data 
ifft_MTI2(ifft_MTI2 < -1000000) = -1000000; % Make sure no -Inf values messing it up

% Define variables for tfridge
fridgeLength = size(ifft_MTI2, 2); % Length of the frequency axis
freqs = linspace(-fs/2, fs/2, 4*N); % Frequency axis for FFT
numRidges = 1; % Number of ridges to extract

% Find fridges
sfft_fridge = fft(MTI2_final, 4*N, 2);
[fridge, ~, ~] = tfridge(rot90(sfft_fridge(:, 1:fridgeLength)), freqs, 1, 'NumRidges', numRidges);
rangeExtracted = 15*(fridge * c * Tp) / delta_r;

% Assembling the figure
figure();
sgtitle('Single target range measurement') 
imagesc(range_data, time, ifft_MTI2);
xlim([0 100]);
xlabel('Range (m)');
ylabel('Time (s)');
colorbar;
clim([-50 0]);
title("Range-spectogram");


