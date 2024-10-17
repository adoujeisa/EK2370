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

range_data = linspace(0, max_range, 6*N);

%% Build matrix for up-chirp and down-chirp data
up_data_parsed = zeros(num_pulse, N);   % Matrix to store up-chirp data
down_data_parsed = zeros(num_pulse, N); % Matrix to store down-chirp data
time = zeros(1, num_pulse);             % Vector to store time for each pulse
pulse_count = 0;                        % Counter for storing data row-by-row

% Loop through sync_pulse to find start of each pulse
for i = 2:(length(sync_pulse) - 6*N)
    % Detect the rising edge in sync signal (new pulse start)
    if sync_pulse(i,1) == 1 && sync_pulse(i-1,1) == 0 
        pulse_count = pulse_count + 1;  % Increment pulse count
        
        % Extract up-chirp and down-chirp data
        up_data_parsed(pulse_count, :) = backscatter(i:i+N-1);          % Up-chirp data
        down_data_parsed(pulse_count, :) = backscatter(i+N:i+2*N-1);    % Down-chirp data
        
        % Record the time for the current pulse
        time(pulse_count) = (i-1) / fs;  % Time in seconds
    end
end

% Remove any extra pre-allocated rows in case not all num_pulse were used
up_data_parsed = up_data_parsed(1:pulse_count, :);
down_data_parsed = down_data_parsed(1:pulse_count, :);
time = time(1:pulse_count);

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

% MTI2 FFT (转换到频域)
sfft_MTI2 = fft(MTI2_final, 6*N, 2);  

%% 3-step MTI
MTI3 = zeros(size(up_data_parsed));  % Initialize 3-step MTI matrix

for t = 3:size(up_data_parsed,1)
    MTI3(t,:) = up_data_parsed(t,:) - 3/2 * up_data_parsed(t-1,:) + 1/2 * up_data_parsed(t-2,:);  
end


% IFFT 部分用于时域恢复 (可选操作，不影响最后结果)
ifft_MTI2 = ifft(sfft_MTI2, 6*N, 2);

% 对恢复的时域信号进行进一步处理（如果需要）
% 此处我们继续使用原来的 FFT 结果，而不是 IFFT 结果
% ifft_MTI2_dB = 20*log10(abs(ifft_MTI2));  % 这个部分仅在你希望分析时域信号时使用

% 使用 FFT 结果继续分析
ifft_MTI2_dB = 20*log10(abs(sfft_MTI2));
ifft_MTI2_dB = ifft_MTI2_dB(:, 1:end/2) - max(max(ifft_MTI2_dB)); % Normalize data 
ifft_MTI2_dB(ifft_MTI2_dB < -1000000) = -1000000; % Make sure no -Inf values messing it up

% Define variables for tfridge
fridgeLength = size(ifft_MTI2_dB, 2); % Length of the frequency axis
freqs = linspace(0, fs/2, 6*N); % Frequency axis for FFT
numRidges = 1; % Number of ridges to extract

% Find fridges
[fridge, ~, ~] = tfridge(rot90(sfft_MTI2(:, 1:fridgeLength)), freqs, 1, 'NumRidges', numRidges);
rangeExtracted = 15*(fridge * c * Tp) / delta_r;

% Assembling the figure
figure();
imagesc(range_data, time, ifft_MTI2_dB);
xlim([0 100]);
xlabel('Range (m)');
ylabel('Time (s)');
colorbar;
clim([-50 0]);
title("Range-spectogram after 2 steps MTI with 2N");

%% 3-step MTI with Zero Padding (4*N)
MTI3 = zeros(size(up_data_parsed));  % Initialize 3-step MTI matrix

for t = 3:size(up_data_parsed, 1)
    MTI3(t, :) = up_data_parsed(t, :) - 3/2 * up_data_parsed(t-1, :) + 1/2 * up_data_parsed(t-2, :);  
end

% Perform FFT with zero-padding (size 4*N)
sfft_MTI3 = fft(MTI3, 6*N, 2);  % 4*N zero-padding along the second dimension

% Convert to dB scale
ifft_MTI3_dB = 20*log10(abs(sfft_MTI3));
ifft_MTI3_dB = ifft_MTI3_dB(:, 1:end/2) - max(max(ifft_MTI3_dB)); % Normalize data
ifft_MTI3_dB(ifft_MTI3_dB < -1000000) = -1000000; % Avoid -Inf values

% Define range axis
range_data = linspace(0, max_range, 6*N);

%% Plotting the 3-step MTI result
figure();
imagesc(range_data, time, ifft_MTI3_dB);
xlim([0 100]);  % Adjust range as needed
xlabel('Range (m)');
ylabel('Time (s)');
colorbar;
clim([-50 0]);  % Adjust color limits as needed
title('Range-Spectrogram after 3-step MTI with 6n');

